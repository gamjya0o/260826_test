#!/usr/bin/env python3
"""
CDR Codon Optimization Tool — MULTI-SCHEME VERSION (Kabat / Chothia / IMGT)
============================================================================
Same behavior as cdr_codon_tool.py, but the FR/CDR numbering scheme used to
draw the CDR boundaries is selectable: kabat (default), chothia, or imgt.
Everything else (WT/mutant parsing, multi-domain LC+HC support, indel-aware
codon replacement, Excel export) works exactly the same as the original tool.

Usage:
    python3 cdr_codon_tool_multischeme.py <input_file> [--scheme kabat|chothia|imgt]

    (scheme defaults to kabat if omitted, so existing workflows are unaffected)
"""

import re
import sys
from anarci import run_anarci
from codon_table import CODON_TABLE

NT_CHARS = set("ACGTUN")

# ---------------------------------------------------------------------------
# CDR boundaries per numbering scheme (standard definitions)
# IMGT uses a single unified numbering/boundary set for both light and heavy
# chains; Kabat and Chothia differ between light and heavy.
# ---------------------------------------------------------------------------
SCHEME_CDR_RANGES = {
    "kabat": {
        "light": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
        "heavy": {"CDR1": (31, 35), "CDR2": (50, 65), "CDR3": (95, 102)},
    },
    "chothia": {
        "light": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
        "heavy": {"CDR1": (26, 32), "CDR2": (52, 56), "CDR3": (95, 102)},
    },
    "imgt": {
        "light": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
        "heavy": {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)},
    },
}
VALID_SCHEMES = list(SCHEME_CDR_RANGES.keys())
DEFAULT_SCHEME = "kabat"

FRAME_ORDER = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(nt_seq):
    nt_seq = nt_seq.upper()
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        aa.append(CODON_TO_AA.get(codon, "X"))
    return "".join(aa)


# ---------------------------------------------------------------------------
# Input-file parsing
# ---------------------------------------------------------------------------
def parse_records(content):
    """Split a multi-record text file into [(header, sequence_no_whitespace), ...]."""
    records = []
    cur_header, cur_seq = None, []
    for raw in content.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">") or line.endswith(">>"):
            if cur_header is not None:
                records.append((cur_header, "".join(cur_seq)))
            cur_header = line.lstrip(">").rstrip(">").strip()
            cur_seq = []
        else:
            cur_seq.append(re.sub(r"\s+", "", line))
    if cur_header is not None:
        records.append((cur_header, "".join(cur_seq)))
    return records


def looks_like_nucleotide(seq):
    if not seq:
        return False
    letters = set(seq.upper())
    return letters.issubset(NT_CHARS)


def parse_blocks_positional(content):
    """
    Fallback parser used when the input has NO 'wild'-labeled header at all
    (i.e. header lines like 'wild fasta>>' were omitted entirely). Splits the
    file into blocks separated by one or more blank lines. Each block may
    still optionally start with a header line (">" / "...>>") — if present
    it's used as the record's name; if absent, the block is left unnamed
    (auto-named later). Multiple lines within a block (e.g. LC line + HC
    line) are concatenated with no separator, exactly like the header-based
    parser does.
    """
    raw_blocks = re.split(r"\n\s*\n+", content.strip())
    blocks = []
    for raw_block in raw_blocks:
        lines = [l.strip() for l in raw_block.splitlines() if l.strip()]
        if not lines:
            continue
        header = None
        seq_lines = lines
        if lines[0].startswith(">") or lines[0].endswith(">>"):
            header = lines[0].lstrip(">").rstrip(">").strip()
            seq_lines = lines[1:]
        seq = "".join(re.sub(r"\s+", "", l) for l in seq_lines)
        if seq:
            blocks.append((header, seq))
    return blocks


def chain_category(chain_type):
    """Collapse ANARCI chain types into 'light' or 'heavy' for WT<->mutant matching."""
    return "light" if chain_type in ("K", "L") else "heavy"


def extract_wt_pairs(content, scheme=DEFAULT_SCHEME):
    """
    Returns (wt_pairs, mutants).

    Two input styles are supported:

    1) HEADER-LABELED (original style) — any record whose header contains
       'wild' (case-insensitive) is collected as WT; nucleotide-looking and
       amino-acid-looking ones are paired by translation match. Every other
       header-labeled record is a mutant. This path is used whenever at
       least one 'wild'-labeled header is found anywhere in the file.

           wild fasta>>
           DIQMT...
           EVQLV...

           wild nucleotide>>
           GACAT...
           GAGGT...

           mutant1>>
           DIQMT...
           EVQLV...

    2) HEADER-OPTIONAL (positional style) — used automatically when NO
       header anywhere contains 'wild' (headers can be present or absent
       for individual blocks; this just means none of them say "wild").
       Blocks are separated by a blank line. By convention the FIRST block
       is the WT amino-acid (FASTA) sequence, the SECOND block is the WT
       nucleotide sequence, and every block after that is a mutant (in
       order). Within a block, put LC then HC on separate lines (no header
       needed) — they'll be concatenated and ANARCI will still detect both
       domains automatically:

           DIQMT...
           EVQLV...

           GACAT...
           GAGGT...

           DIQMT...(mutant1 LC)
           EVQLV...(mutant1 HC)

           DIQMT...(mutant2 LC)
           EVQLV...(mutant2 HC)

    mutants: list of (id, aa_seq).
    """
    records = parse_records(content)
    wild_records = [(h, s) for h, s in records if h and "wild" in h.lower()]

    if wild_records:
        return _extract_wt_pairs_headered(records, scheme)
    return _extract_wt_pairs_positional(content, scheme)


def _extract_wt_pairs_headered(records, scheme):
    wild_records = [(h, s) for h, s in records if h and "wild" in h.lower()]
    other_records = [(h, s) for h, s in records if not (h and "wild" in h.lower())]

    nt_candidates = [(h, s) for h, s in wild_records if looks_like_nucleotide(s)]
    aa_candidates = [(h, s) for h, s in wild_records if (h, s) not in nt_candidates]

    if not aa_candidates or not nt_candidates:
        raise ValueError(
            "입력 파일에서 wild type의 아미노산(FASTA) 서열과 nucleotide(codon) 서열을 "
            "모두 찾지 못했습니다. 헤더에 'wild'가 포함된 레코드가 최소 하나의 아미노산 "
            "서열과 하나의 nucleotide 서열 쌍으로 있어야 합니다. "
            "예: 'wild fasta>>' / 'wild nucleotide>>' (LC/HC처럼 여러 쌍도 가능). "
            "또는 헤더를 아예 생략하고 싶다면 'wild'라는 단어가 포함된 헤더를 전부 빼고, "
            "첫 블록=WT fasta, 둘째 블록=WT nucleotide 순서로 빈 줄로만 구분해 넣어도 됩니다."
        )

    used_nt_idx = set()
    wt_pairs = []
    for aa_h, aa_s in aa_candidates:
        aa_s_up = aa_s.upper()
        chosen = None
        # 1) exact translation match (most reliable)
        for i, (nt_h, nt_s) in enumerate(nt_candidates):
            if i in used_nt_idx:
                continue
            if translate(nt_s) == aa_s_up:
                chosen = i
                break
        # 2) fallback: matching length (3x) if no exact translation match found
        if chosen is None:
            for i, (nt_h, nt_s) in enumerate(nt_candidates):
                if i in used_nt_idx:
                    continue
                if len(nt_s) == len(aa_s_up) * 3:
                    chosen = i
                    sys.stderr.write(
                        f"[WARNING] '{aa_h}' <-> '{nt_h}' 매칭: nucleotide 번역 결과가 "
                        "아미노산 서열과 정확히 일치하지 않지만, 길이가 같아 자동으로 짝지었습니다. "
                        "직접 확인해 주세요.\n"
                    )
                    break
        if chosen is None:
            raise ValueError(f"'{aa_h}' 에 대응하는 wild nucleotide 서열을 찾지 못했습니다.")

        used_nt_idx.add(chosen)
        nt_h, nt_s = nt_candidates[chosen]
        domains = get_all_domains(aa_s_up, aa_h, scheme=scheme)
        wt_pairs.append({
            "id": aa_h,
            "aa": aa_s_up,
            "nt": nt_s.upper(),
            "chain_type": domain_summary(domains),
            "category": chain_category(domains[0]["chain_type"]),
        })

    mutants = [(h, s) for h, s in other_records if s]
    return wt_pairs, mutants


def _extract_wt_pairs_positional(content, scheme):
    blocks = parse_blocks_positional(content)
    if len(blocks) < 2:
        raise ValueError(
            "입력 파일에서 최소 2개 블록(WT fasta, WT nucleotide)을 찾지 못했습니다. "
            "헤더 없이 넣을 경우, 블록은 반드시 빈 줄로 구분되어야 하고 "
            "첫 블록=WT fasta, 둘째 블록=WT nucleotide 순서여야 합니다."
        )

    wt_aa_header, wt_aa_seq = blocks[0]
    wt_nt_header, wt_nt_seq = blocks[1]

    # be forgiving about order: if block 2 isn't nucleotide but block 1 is, swap
    if not looks_like_nucleotide(wt_nt_seq) and looks_like_nucleotide(wt_aa_seq):
        wt_aa_header, wt_aa_seq, wt_nt_header, wt_nt_seq = wt_nt_header, wt_nt_seq, wt_aa_header, wt_aa_seq

    if not looks_like_nucleotide(wt_nt_seq):
        raise ValueError(
            "헤더 없이 입력하신 경우, 첫 번째 블록은 WT의 아미노산(FASTA) 서열, "
            "두 번째 블록은 WT의 nucleotide 서열이어야 합니다. 지금은 두 번째 블록이 "
            "nucleotide(ACGT) 서열로 인식되지 않습니다 — 순서를 확인해 주시거나, "
            "헤더에 'wild'를 포함시켜 명시적으로 표시해 주세요 (예: 'wild fasta>>')."
        )

    wt_aa_up = wt_aa_seq.upper()
    wt_nt_up = wt_nt_seq.upper()
    translated = translate(wt_nt_up)
    if translated != wt_aa_up:
        sys.stderr.write(
            "[WARNING] wild nucleotide 서열을 번역한 결과가 wild FASTA 아미노산 서열과 "
            f"일치하지 않습니다.\n  FASTA     : {wt_aa_up}\n  Translated: {translated}\n"
        )

    wt_id = wt_aa_header or "wild"
    domains = get_all_domains(wt_aa_up, wt_id, scheme=scheme)
    wt_pairs = [{
        "id": wt_id,
        "aa": wt_aa_up,
        "nt": wt_nt_up,
        "chain_type": domain_summary(domains),
        "category": chain_category(domains[0]["chain_type"]),
    }]

    mutants = []
    for i, (h, s) in enumerate(blocks[2:], start=1):
        mutant_id = h if h else f"mutant{i}"
        mutants.append((mutant_id, s))

    return wt_pairs, mutants


def pick_wt_for_mutant(wt_pairs, mutant_id, mutant_aa, scheme=DEFAULT_SCHEME):
    """
    Choose which WT pair a mutant should be compared against.
      - If only one WT pair exists, use it (single-chain mode, backward compatible).
      - Otherwise, number the mutant with ANARCI, get its chain category
        (light/heavy), and match to WT pair(s) of the same category.
        If more than one WT candidate shares that category, pick the one
        with the highest amino-acid identity to the mutant.
    """
    if len(wt_pairs) == 1:
        return wt_pairs[0]

    _, mut_chain_type, _ = get_kabat_numbering(mutant_aa.strip().upper(), mutant_id, scheme=scheme)
    mut_category = chain_category(mut_chain_type)
    candidates = [wt for wt in wt_pairs if wt["category"] == mut_category]

    if not candidates:
        raise ValueError(
            f"'{mutant_id}' (chain type {mut_chain_type})에 대응하는 wild type을 찾지 못했습니다. "
            f"등록된 WT: {[(w['id'], w['chain_type']) for w in wt_pairs]}"
        )
    if len(candidates) == 1:
        return candidates[0]

    # multiple same-category WT candidates -> pick best identity match
    def identity_score(wt):
        n = min(len(wt["aa"]), len(mutant_aa))
        matches = sum(1 for a, b in zip(wt["aa"][:n], mutant_aa[:n]) if a == b)
        return matches / max(n, 1)

    best = max(candidates, key=identity_score)
    sys.stderr.write(
        f"[NOTE] '{mutant_id}'을(를) 동일 계열({mut_category}) WT 후보 중 "
        f"'{best['id']}'와(과) 가장 유사하다고 판단해 매칭했습니다.\n"
    )
    return best


# ---------------------------------------------------------------------------
# ANARCI helpers (multi-domain aware: supports LC+HC or any tandem chains
# concatenated into ONE sequence, e.g. scFv-style records)
# ---------------------------------------------------------------------------
def get_all_domains(seq, seq_id="seq", scheme=DEFAULT_SCHEME):
    """
    Run ANARCI (given numbering scheme) and return a list of domains found in
    `seq`, sorted by their start position:
        [{"numbering": [...], "chain_type": "K"/"L"/"H", "start": int, "end": int}, ...]
    A sequence with just one V-domain (e.g. a lone light or heavy chain)
    returns a list of length 1. A sequence with an LC and an HC concatenated
    together returns a list of length 2 (ANARCI detects each domain natively).
    """
    out = run_anarci([(seq_id, seq)], scheme=scheme)
    numbered, aligned = out[1], out[2]
    if numbered[0] is None or not numbered[0]:
        raise ValueError(f"ANARCI could not number sequence '{seq_id}' as an antibody variable domain.")
    domains = []
    for i, dom in enumerate(numbered[0]):
        numbering, start, end = dom
        chain_type = aligned[0][i]["chain_type"]
        domains.append({"numbering": numbering, "chain_type": chain_type, "start": start, "end": end})
    domains.sort(key=lambda d: d["start"])
    return domains


def get_kabat_numbering(seq, seq_id="seq", scheme=DEFAULT_SCHEME):
    """Backward-compatible helper: returns the FIRST domain only
    (numbering_list, chain_type, query_start)."""
    domains = get_all_domains(seq, seq_id, scheme=scheme)
    d = domains[0]
    return d["numbering"], d["chain_type"], d["start"]


def domain_summary(domains):
    return "+".join(d["chain_type"] for d in domains)


def build_kabat_to_seqidx(numbering, start):
    """Map (num, insertion_code) -> 0-based index in the ORIGINAL ungapped sequence."""
    mapping = {}
    seq_idx = start
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        mapping[(num, ins)] = seq_idx
        seq_idx += 1
    return mapping


def region_of(num, chain_type, scheme=DEFAULT_SCHEME):
    ranges = SCHEME_CDR_RANGES[scheme]["light" if chain_type in ("K", "L") else "heavy"]
    c1, c2, c3 = ranges["CDR1"], ranges["CDR2"], ranges["CDR3"]
    if num < c1[0]:
        return "FR1"
    if c1[0] <= num <= c1[1]:
        return "CDR1"
    if c1[1] < num < c2[0]:
        return "FR2"
    if c2[0] <= num <= c2[1]:
        return "CDR2"
    if c2[1] < num < c3[0]:
        return "FR3"
    if c3[0] <= num <= c3[1]:
        return "CDR3"
    return "FR4"


def which_cdr(num, chain_type, scheme=DEFAULT_SCHEME):
    ranges = SCHEME_CDR_RANGES[scheme]["light" if chain_type in ("K", "L") else "heavy"]
    L = "L" if chain_type in ("K", "L") else "H"
    for i, key in enumerate(["CDR1", "CDR2", "CDR3"], start=1):
        s, e = ranges[key]
        if s <= num <= e:
            return f"CDR-{L}{i}"
    return "FR (framework)"


# ---------------------------------------------------------------------------
# FR/CDR framing (multi-domain: e.g. LC domain followed by HC domain)
# ---------------------------------------------------------------------------
def frame_one_domain(seq, domain, scheme=DEFAULT_SCHEME):
    """FR1/CDR1/.../FR4 segments for a single domain dict, using absolute seq indices."""
    numbering, chain_type, start = domain["numbering"], domain["chain_type"], domain["start"]
    pos2idx = build_kabat_to_seqidx(numbering, start)
    ordered = sorted(pos2idx.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    segments = {name: [] for name in FRAME_ORDER}
    max_idx_used = -1
    for (num, ins), idx in ordered:
        seg = region_of(num, chain_type, scheme)
        segments[seg].append(seq[idx])
        max_idx_used = max(max_idx_used, idx)
    seg_strs = {name: "".join(chars) for name, chars in segments.items()}
    space_joined = " ".join(seg_strs[name] for name in FRAME_ORDER)
    return seg_strs, space_joined, chain_type, max_idx_used


def frame_sequence(seq, seq_id="seq", scheme=DEFAULT_SCHEME):
    """
    Number `seq` with ANARCI/Kabat (possibly multiple domains, e.g. LC+HC
    concatenated) and split EACH domain into FR1/CDR1/.../FR4.
    Returns: domain_results, seq_id
      domain_results = [
        {"chain_type": "K", "seg_strs": {...}, "space_joined": "...",
         "leading_gap": "", "start": 0, "end": 106},
        ...
      ]
    Anything AFTER a domain's official FR4 (a linker before the next domain,
    a C-terminal tag, or a constant region) is appended directly onto that
    domain's FR4 rather than being reported separately — i.e. FR4 always
    runs to the start of the next domain (or the end of the sequence for the
    last domain). Residues BEFORE the very first domain (rare — e.g. an
    N-terminal tag) are still reported separately via "leading_gap" on that
    first domain.
    """
    seq = seq.strip().upper()
    domains = get_all_domains(seq, seq_id, scheme=scheme)

    domain_results = []
    prev_end = -1
    for i, dom in enumerate(domains):
        leading_gap = seq[prev_end + 1: dom["start"]] if dom["start"] > prev_end + 1 else ""
        seg_strs, space_joined, chain_type, max_idx_used = frame_one_domain(seq, dom, scheme)

        # Absorb everything up to the next domain (or end of sequence) into FR4
        next_start = domains[i + 1]["start"] if i + 1 < len(domains) else len(seq)
        if next_start > max_idx_used + 1:
            seg_strs["FR4"] += seq[max_idx_used + 1: next_start]
            space_joined = " ".join(seg_strs[name] for name in FRAME_ORDER)
            max_idx_used = next_start - 1

        domain_results.append({
            "chain_type": chain_type,
            "seg_strs": seg_strs,
            "space_joined": space_joined,
            "leading_gap": leading_gap,
            "start": dom["start"],
            "end": max_idx_used,
        })
        prev_end = max_idx_used

    trailing_tail = seq[prev_end + 1:] if prev_end + 1 < len(seq) else ""
    return domain_results, trailing_tail


def format_framing(domain_results, trailing_tail, seq_id="seq", scheme=DEFAULT_SCHEME):
    lines = [f"[{seq_id}] domain(s) found: {domain_summary([{'chain_type': d['chain_type']} for d in domain_results])}  (scheme: {scheme})"]
    for i, d in enumerate(domain_results, start=1):
        L = "L" if d["chain_type"] in ("K", "L") else "H"
        label = "Light" if L == "L" else "Heavy"
        lines.append(f"--- Domain {i} ({label} chain, type {d['chain_type']}) ---")
        if d["leading_gap"]:
            lines.append(f"  (leading/linker, outside V-domain numbering): {d['leading_gap']}")
        seg_strs = d["seg_strs"]
        lines.append(f"  FR1      : {seg_strs['FR1']}")
        lines.append(f"  CDR-{L}1   : {seg_strs['CDR1']}")
        lines.append(f"  FR2      : {seg_strs['FR2']}")
        lines.append(f"  CDR-{L}2   : {seg_strs['CDR2']}")
        lines.append(f"  FR3      : {seg_strs['FR3']}")
        lines.append(f"  CDR-{L}3   : {seg_strs['CDR3']}")
        lines.append(f"  FR4      : {seg_strs['FR4']}")
        lines.append(f"  space-joined: {d['space_joined']}")
    if trailing_tail:
        lines.append(f"  (tail, outside V-domain numbering): {trailing_tail}")
    full_joined = " | ".join(d["space_joined"] for d in domain_results)
    if trailing_tail:
        full_joined += f" | {trailing_tail}"
    lines.append(f"  FULL space-joined (all domains): {full_joined}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Mutation detection + CDR codon replacement (multi-domain aware)
# ---------------------------------------------------------------------------
def analyze(wt_id, wt_aa, wt_nt, mutant_aa_seq, mutant_id="MUT", scheme=DEFAULT_SCHEME):
    assert len(wt_aa) * 3 == len(wt_nt), "WT amino acid / nucleotide length mismatch!"
    mutant_aa_seq = mutant_aa_seq.strip().upper()

    wt_domains = get_all_domains(wt_aa, wt_id, scheme=scheme)
    mut_domains = get_all_domains(mutant_aa_seq, mutant_id, scheme=scheme)

    domain_count_note = None
    if len(wt_domains) != len(mut_domains):
        domain_count_note = (
            f"WT에서는 도메인 {len(wt_domains)}개(체인: {domain_summary(wt_domains)}), "
            f"'{mutant_id}'에서는 도메인 {len(mut_domains)}개(체인: {domain_summary(mut_domains)})가 "
            f"발견되어 개수가 다릅니다. 앞에서부터 짝지어지는 도메인까지만 비교합니다."
        )

    n_domains = min(len(wt_domains), len(mut_domains))
    codon_list = [wt_nt[i:i + 3] for i in range(0, len(wt_nt), 3)]
    deleted_wt_idx = set()          # WT codon slots to remove entirely (CDR deletions)
    insertions_after = {}           # wt_idx -> [codon, codon, ...] to splice in right after that slot
    mutations, applied, skipped = [], [], []

    for d in range(n_domains):
        wt_dom, mut_dom = wt_domains[d], mut_domains[d]
        if wt_dom["chain_type"] != mut_dom["chain_type"]:
            sys.stderr.write(
                f"[WARNING] '{mutant_id}' 도메인 {d+1}: WT chain type "
                f"{wt_dom['chain_type']} vs mutant chain type {mut_dom['chain_type']} 불일치.\n"
            )
        wt_pos2idx = build_kabat_to_seqidx(wt_dom["numbering"], wt_dom["start"])
        mut_pos2idx = build_kabat_to_seqidx(mut_dom["numbering"], mut_dom["start"])
        wt_map = {k: wt_aa[i] for k, i in wt_pos2idx.items()}
        mut_map = {k: mutant_aa_seq[i] for k, i in mut_pos2idx.items()}
        all_kabat_positions = sorted(set(wt_map) | set(mut_map), key=lambda x: (x[0], x[1]))

        last_wt_idx = wt_dom["start"] - 1  # anchor for any insertion before the first WT residue
        for kpos in all_kabat_positions:
            wt_res = wt_map.get(kpos, "-")
            mut_res = mut_map.get(kpos, "-")
            wt_seq_idx = wt_pos2idx.get(kpos)
            mut_seq_idx = mut_pos2idx.get(kpos)

            if wt_res != mut_res:
                num, ins = kpos
                label = f"{num}{ins.strip()}"
                region = which_cdr(num, wt_dom["chain_type"], scheme)
                if wt_res != "-" and mut_res != "-":
                    mtype = "substitution"
                elif wt_res != "-" and mut_res == "-":
                    mtype = "deletion"
                else:
                    mtype = "insertion"
                mutations.append({
                    "domain": d + 1,
                    "chain_type": wt_dom["chain_type"],
                    "kabat_tuple": kpos,
                    "kabat_position": label,
                    "region": region,
                    "wt_aa": wt_res,
                    "mut_aa": mut_res,
                    "wt_seq_index": wt_seq_idx,
                    "mut_seq_index": mut_seq_idx,
                    "type": mtype,
                    "anchor_wt_idx": last_wt_idx,
                })

            if wt_seq_idx is not None:
                last_wt_idx = wt_seq_idx

    for m in mutations:
        if not m["region"].startswith("CDR"):
            continue  # FR changes never touch the nucleotide sequence (by design)
        if m["type"] == "substitution":
            if m["mut_aa"] not in CODON_TABLE:
                skipped.append(m)
                continue
            codon_list[m["wt_seq_index"]] = CODON_TABLE[m["mut_aa"]]
            applied.append(m)
        elif m["type"] == "deletion":
            deleted_wt_idx.add(m["wt_seq_index"])
            applied.append(m)
        elif m["type"] == "insertion":
            if m["mut_aa"] not in CODON_TABLE:
                skipped.append(m)
                continue
            insertions_after.setdefault(m["anchor_wt_idx"], []).append(CODON_TABLE[m["mut_aa"]])
            applied.append(m)

    parts = []
    for i, codon in enumerate(codon_list):
        if i not in deleted_wt_idx:
            parts.append(codon)
        parts.extend(insertions_after.get(i, []))
    parts.extend(insertions_after.get(-1, []))  # insertions anchored before domain start, if any
    new_nt_seq = "".join(parts)
    mut_domain_results, mut_trailing_tail = frame_sequence(mutant_aa_seq, mutant_id, scheme)

    return {
        "scheme": scheme,
        "wt_domain_summary": domain_summary(wt_domains),
        "mut_domain_summary": domain_summary(mut_domains),
        "domain_count_note": domain_count_note,
        "mutations": mutations,
        "cdr_mutations_applied": applied,
        "cdr_mutations_skipped": skipped,
        "modified_nucleotide_seq": new_nt_seq,
        "mut_domain_results": mut_domain_results,
        "mut_trailing_tail": mut_trailing_tail,
        "wt_domains": wt_domains,
        "mut_domains": mut_domains,
        "mutant_aa_seq": mutant_aa_seq,
        "codon_list": codon_list,
        "deleted_wt_idx": deleted_wt_idx,
        "insertions_after": insertions_after,
    }


def format_report(wt_nt, result, mutant_id="MUT"):
    lines = []
    lines.append(f"Numbering scheme: {result['scheme']}")
    lines.append(f"WT domains: {result['wt_domain_summary']}   MUT domains: {result['mut_domain_summary']}")
    if result["domain_count_note"]:
        lines.append(f"[NOTE] {result['domain_count_note']}")

    lines.append("")
    lines.append(f"=== FR/CDR framing ({mutant_id}) ===")
    lines.append(format_framing(result["mut_domain_results"], result["mut_trailing_tail"], mutant_id, result["scheme"]))

    lines.append("")
    lines.append(f"=== Amino acid differences ({result['scheme']} numbering) ===")
    if not result["mutations"]:
        lines.append("(no differences found)")
    for m in result["mutations"]:
        lines.append(f"  [Domain {m['domain']}/{m['chain_type']}] Kabat {m['kabat_position']:>5}  "
                     f"[{m['region']:<15}]  {m['wt_aa']} -> {m['mut_aa']}")

    lines.append("")
    lines.append("=== CDR codon changes applied (per codon optimization table) ===")
    if not result["cdr_mutations_applied"]:
        lines.append("(none)")
    for m in result["cdr_mutations_applied"]:
        if m["type"] == "substitution":
            old_codon = wt_nt[m["wt_seq_index"] * 3: m["wt_seq_index"] * 3 + 3]
            new_codon = CODON_TABLE[m["mut_aa"]]
            lines.append(f"  [Domain {m['domain']}/{m['chain_type']}] Kabat {m['kabat_position']:>5}  "
                         f"[{m['region']}]  {m['wt_aa']}({old_codon}) -> {m['mut_aa']}({new_codon})   "
                         f"[nt pos {m['wt_seq_index']*3+1}-{m['wt_seq_index']*3+3}]")
        elif m["type"] == "deletion":
            old_codon = wt_nt[m["wt_seq_index"] * 3: m["wt_seq_index"] * 3 + 3]
            lines.append(f"  [Domain {m['domain']}/{m['chain_type']}] Kabat {m['kabat_position']:>5}  "
                         f"[{m['region']}]  DELETION: {m['wt_aa']}({old_codon}) removed   "
                         f"[nt pos {m['wt_seq_index']*3+1}-{m['wt_seq_index']*3+3}]")
        elif m["type"] == "insertion":
            new_codon = CODON_TABLE[m["mut_aa"]]
            lines.append(f"  [Domain {m['domain']}/{m['chain_type']}] Kabat {m['kabat_position']:>5}  "
                         f"[{m['region']}]  INSERTION: {m['mut_aa']}({new_codon}) added after WT nt pos "
                         f"{m['anchor_wt_idx']*3+3}")

    if result["cdr_mutations_skipped"]:
        lines.append("")
        lines.append("=== CDR changes NOT applied (no WT codon slot / no table entry) ===")
        for m in result["cdr_mutations_skipped"]:
            lines.append(f"  [Domain {m['domain']}/{m['chain_type']}] Kabat {m['kabat_position']:>5}  "
                         f"{m['wt_aa']} -> {m['mut_aa']}")

    lines.append("")
    lines.append("=== Final nucleotide sequence (WT backbone, CDR codons swapped) ===")
    lines.append(result["modified_nucleotide_seq"])
    return "\n".join(lines)


def run_file(path_or_text, xlsx_output=None, scheme=DEFAULT_SCHEME):
    if scheme not in VALID_SCHEMES:
        raise ValueError(f"지원하지 않는 scheme '{scheme}' 입니다. 사용 가능: {VALID_SCHEMES}")

    try:
        with open(path_or_text) as f:
            content = f.read()
        base = re.sub(r"\.[^.]*$", "", path_or_text)
        default_xlsx = f"{base}_result_{scheme}.xlsx" if scheme != DEFAULT_SCHEME else f"{base}_result.xlsx"
    except (FileNotFoundError, OSError):
        content = path_or_text
        default_xlsx = f"cdr_codon_result_{scheme}.xlsx"

    wt_pairs, mutants = extract_wt_pairs(content, scheme=scheme)

    for wt in wt_pairs:
        print(f"{'#'*20} WILD TYPE ({wt['id']}) {'#'*20}")
        domain_results, trailing_tail = frame_sequence(wt["aa"], wt["id"], scheme)
        print(format_framing(domain_results, trailing_tail, wt["id"], scheme))
        print()

    mutant_results = []
    for mutant_id, mutant_seq in mutants:
        print(f"{'='*20} {mutant_id} {'='*20}")
        wt = pick_wt_for_mutant(wt_pairs, mutant_id, mutant_seq, scheme)
        if len(wt_pairs) > 1:
            print(f"(matched against WT: {wt['id']})")
        result = analyze(wt["id"], wt["aa"], wt["nt"], mutant_seq, mutant_id=mutant_id, scheme=scheme)
        print(format_report(wt["nt"], result, mutant_id=mutant_id))
        print()
        mutant_results.append((mutant_id, result["mutant_aa_seq"], result, wt))

    try:
        import excel_export_multischeme as excel_export
        xlsx_path = xlsx_output or default_xlsx
        excel_export.build_workbook(wt_pairs, mutant_results, xlsx_path, scheme=scheme)
        print(f"[Excel 파일 저장됨] {xlsx_path}")
    except Exception as e:
        sys.stderr.write(f"[WARNING] Excel 파일 생성에 실패했습니다: {e}\n")


def _parse_cli_args(argv):
    """Minimal CLI parsing: <input_file> [--scheme kabat|chothia|imgt]"""
    if not argv:
        print("Usage: python3 cdr_codon_tool_multischeme.py <input_file> [--scheme kabat|chothia|imgt]")
        sys.exit(1)
    input_file = argv[0]
    scheme = DEFAULT_SCHEME
    args = argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--scheme" and i + 1 < len(args):
            scheme = args[i + 1].lower()
            i += 2
        else:
            i += 1
    if scheme not in VALID_SCHEMES:
        print(f"[ERROR] --scheme 값은 {VALID_SCHEMES} 중 하나여야 합니다. (입력값: '{scheme}')")
        sys.exit(1)
    return input_file, scheme


if __name__ == "__main__":
    input_file, scheme = _parse_cli_args(sys.argv[1:])
    run_file(input_file, scheme=scheme)
