# -*- coding: utf-8 -*-
"""
CDR/Frame 구분 & Codon 최적화 툴 — 웹 버전 (Streamlit)

WSL / Ubuntu / anarci 설치 없이, 브라우저에서 파일 업로드 -> 버튼 클릭 ->
결과 확인 + 엑셀 다운로드까지 되는 웹 앱입니다.
로직(cdr_codon_tool.py / excel_export.py / codon_table.py)은 기존 것을 그대로 사용합니다.
"""

import io
import sys
import contextlib
import traceback

import streamlit as st

from cdr_codon_tool import (
    VALID_SCHEMES, DEFAULT_SCHEME,
    extract_wt_pairs, pick_wt_for_mutant, analyze, format_report,
    frame_sequence, format_framing,
)
import excel_export

st.set_page_config(page_title="CDR Codon Optimization Tool", page_icon="🧬", layout="wide")

st.title("🧬 CDR/Frame 구분 & Codon 최적화 툴")
st.caption("Kabat / Chothia / IMGT 기준 CDR-Frame 자동 구분 + WT→Mutant CDR codon 치환 + 엑셀 결과 생성")

with st.expander("📖 사용법 (처음이시면 먼저 읽어주세요)", expanded=True):
    st.markdown(
        "1. **맨 처음 블록**에 wild type(야생형)의 **아미노산(fasta) 서열**을 넣습니다. "
        "LC와 HC가 있다면 **Enter(줄바꿈)로 구분해서 두 줄**로 넣으면 자동으로 인식됩니다 "
        "(한 서열로 이어붙여도 됨 — scFv 등).\n"
        "2. **한 줄 띄우고(빈 줄)**, 그 다음 블록에 같은 wild type의 **nucleotide(codon) 서열**을 "
        "넣습니다. 마찬가지로 LC/HC는 줄바꿈으로 구분합니다.\n"
        "3. **다시 한 줄 띄우고**, 그 다음 블록부터는 mutant들을 하나씩 넣습니다 "
        "(mutant는 아미노산 서열만 필요합니다). mutant도 여러 개면 **블록마다 한 줄씩 띄워서** 구분합니다.\n"
        "4. `>wild fasta`, `>wild nucleotide`, `>Mutant-01` 같은 **제목(header) 줄은 전부 생략 가능**"
        "합니다. 제목을 안 쓰면 순서대로 '① 첫 아미노산 블록 = WT / ② 첫 nucleotide 블록 = WT / "
        "③ 그 다음부터는 전부 mutant'로 자동 인식됩니다. 제목을 쓰고 싶은 경우에는 기존처럼 "
        "`>` 로 시작하는 줄을 블록 맨 앞에 추가하면 됩니다 (제목 유무를 블록별로 섞어 써도 됩니다)."
    )
    st.code(
        "DIQMTQSPSS...   (LC, wild type 아미노산)\n"
        "EVQLVESG...     (HC, 엔터로 구분)\n"
        "\n"
        "GACATCC...      (LC, wild type nucleotide)\n"
        "GAGGTGC...      (HC)\n"
        "\n"
        "DIQMTQSPSS...   (mutant 1, 아미노산만)\n"
        "EVQLVESG...\n"
        "\n"
        "DIQMTQSPSS...   (mutant 2)\n"
        "EVQLVESG...",
        language="text",
    )

with st.sidebar:
    st.header("설정")
    scheme = st.selectbox(
        "Numbering scheme",
        options=VALID_SCHEMES,
        index=VALID_SCHEMES.index(DEFAULT_SCHEME),
        help="CDR/FR 경계를 나누는 기준입니다. 잘 모르겠으면 kabat(기본값) 그대로 두세요.",
    )
    st.markdown("---")
    st.markdown(
        "**입력 형식 요약**\n\n"
        "- 블록(빈 줄로 구분): ① WT 아미노산 → ② WT nucleotide → ③ mutant들\n"
        "- 같은 블록 안에서 LC/HC는 줄바꿈(Enter)으로 구분\n"
        "- `>wild fasta` 같은 제목 줄은 **없어도 됩니다** (있으면 그대로 사용)\n"
        "- 자세한 예시는 위쪽 **'사용법'**을 펼쳐서 확인하세요."
    )

uploaded = st.file_uploader("입력 파일 업로드 (.txt / .fasta)", type=["txt", "fasta", "fa"])
pasted = st.text_area("또는 여기에 직접 붙여넣기", height=200, placeholder=">wild fasta\n...\n>wild nucleotide\n...\n>Mutant-01\n...")

content = None
if uploaded is not None:
    content = uploaded.read().decode("utf-8")
elif pasted.strip():
    content = pasted

run = st.button("🚀 분석 실행", type="primary", use_container_width=True, disabled=content is None)

if run and content:
    log_buf = io.StringIO()
    try:
        with contextlib.redirect_stderr(log_buf):
            wt_pairs, mutants = extract_wt_pairs(content, scheme=scheme)

            if not mutants:
                st.warning("mutant 서열이 하나도 발견되지 않았습니다. 입력 파일을 확인해 주세요.")
            else:
                report_sections = []
                mutant_results = []

                st.subheader("Wild Type")
                for wt in wt_pairs:
                    domain_results, trailing_tail = frame_sequence(wt["aa"], wt["id"], scheme)
                    txt = format_framing(domain_results, trailing_tail, wt["id"], scheme)
                    with st.expander(f"WT: {wt['id']}", expanded=False):
                        st.code(txt, language="text")
                    report_sections.append(txt)

                st.subheader("Mutants")
                progress = st.progress(0.0)
                for i, (mutant_id, mutant_seq) in enumerate(mutants):
                    wt = pick_wt_for_mutant(wt_pairs, mutant_id, mutant_seq, scheme)
                    result = analyze(wt["id"], wt["aa"], wt["nt"], mutant_seq, mutant_id=mutant_id, scheme=scheme)
                    report_txt = format_report(wt["nt"], result, mutant_id=mutant_id)

                    n_applied = len(result["cdr_mutations_applied"])
                    n_total = len(result["mutations"])
                    with st.expander(f"{mutant_id}  —  변이 {n_total}개 (CDR 반영 {n_applied}개)", expanded=False):
                        st.code(report_txt, language="text")

                    report_sections.append(report_txt)
                    mutant_results.append((mutant_id, result["mutant_aa_seq"], result, wt))
                    progress.progress((i + 1) / len(mutants))

                # Build excel in-memory
                xlsx_buf = io.BytesIO()
                tmp_path = "/tmp/_cdr_result.xlsx"
                excel_export.build_workbook(wt_pairs, mutant_results, tmp_path, scheme=scheme)
                with open(tmp_path, "rb") as f:
                    xlsx_buf.write(f.read())
                xlsx_buf.seek(0)

                st.success(f"완료! WT {len(wt_pairs)}개, Mutant {len(mutants)}개 처리했습니다.")

                col1, col2 = st.columns(2)
                with col1:
                    st.download_button(
                        "📊 엑셀 결과 다운로드 (.xlsx)",
                        data=xlsx_buf,
                        file_name=f"cdr_codon_result_{scheme}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                    )
                with col2:
                    full_report = "\n\n".join(report_sections)
                    st.download_button(
                        "📄 텍스트 리포트 다운로드 (.txt)",
                        data=full_report.encode("utf-8"),
                        file_name=f"cdr_codon_report_{scheme}.txt",
                        use_container_width=True,
                    )

        warnings = log_buf.getvalue().strip()
        if warnings:
            with st.expander("⚠️ 경고/참고 메시지"):
                st.code(warnings, language="text")

    except Exception as e:
        st.error(f"오류가 발생했습니다: {e}")
        with st.expander("상세 에러 로그 (개발자용)"):
            st.code(traceback.format_exc(), language="text")
        warnings = log_buf.getvalue().strip()
        if warnings:
            with st.expander("경고 로그"):
                st.code(warnings, language="text")
