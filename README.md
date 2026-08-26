# CDR/Frame 구분 & Codon 최적화 툴 — 웹 버전

기존 파이썬 스크립트(`cdr_codon_tool.py`)와 완전히 동일한 로직을, 브라우저에서
파일 업로드 → 버튼 클릭 → 결과 확인/엑셀 다운로드로 쓸 수 있게 만든 웹 앱입니다.

## 왜 "그냥 HTML 파일 하나"가 아닌가요?

이 도구는 항체 서열을 Kabat/Chothia/IMGT 번호로 정확히 매기기 위해 **ANARCI**라는
생물정보학 프로그램을 쓰는데, 이건 HMMER라는 별도 엔진으로 서열을 정렬하는
무거운 계산을 합니다. 이런 계산은 브라우저(자바스크립트)에서 돌릴 수 없어서,
"인터넷 연결 없이 파일 하나 더블클릭"으로는 만들 수 없습니다.

대신 **서버에서 한 번만 돌아가는 웹 앱**으로 만들었습니다. 실험실 동료들은
그냥 웹 주소(URL)에 접속해서 파일 업로드 → 버튼 클릭만 하면 되고, WSL/Ubuntu/anarci/hmmer
같은 건 아무것도 설치할 필요가 없습니다 (전부 서버 안에 이미 들어있음). 처음 한 번만
"서버에 올리는" 작업이 필요합니다.

## 배포 방법 (추천: Streamlit Community Cloud — 완전 무료)

> ⚠️ 예전엔 Hugging Face Spaces(Docker)도 무료였는데, 2026년 중순부터 Hugging Face가
> Docker/Gradio SDK를 유료(PRO 이상) 계정에서만 만들 수 있게 정책을 바꿨습니다.
> 그래서 대신 **Streamlit이 직접 운영하는 무료 호스팅(Community Cloud)**을 사용합니다.
> 이쪽은 Docker 개념 없이 `requirements.txt`(파이썬 패키지)와 `packages.txt`(hmmer 같은
> 시스템 패키지)만 있으면 자동으로 인식해서 설치해줍니다.

### 1단계 — GitHub에 코드 올리기 (한 번만)

1. https://github.com 에서 무료 계정 생성 (이미 있으면 로그인)
2. 우측 상단 **+ → New repository** 클릭 → 이름 입력(예: `cdr-codon-tool`) →
   Public 또는 Private 선택 → **Create repository**
3. 생성된 저장소 페이지에서 **Add file → Upload files** 클릭
4. 이 폴더 안의 파일 6개를 전부 드래그 앤 드롭으로 업로드:
   - `app.py`
   - `cdr_codon_tool.py`
   - `excel_export.py`
   - `codon_table.py`
   - `requirements.txt`
   - `packages.txt`
5. **Commit changes** 클릭

### 2단계 — Streamlit Community Cloud에 배포하기

1. https://share.streamlit.io 접속 → **Sign in with GitHub**로 로그인 (계정 연동만
   승인하면 됨, 별도 가입 절차 없음)
2. **Create app** (또는 **New app**) 클릭
3. 방금 만든 저장소(`cdr-codon-tool`), 브랜치(`main`), 메인 파일 경로에 `app.py` 지정
4. **Deploy** 클릭 → 자동으로 빌드 시작 (`packages.txt`를 보고 hmmer를 자동 설치,
   `requirements.txt`를 보고 anarci/biopython/openpyxl 등을 자동 설치) — 3~5분 소요
5. 빌드가 끝나면 `https://사용자명-cdr-codon-tool-app-xxxx.streamlit.app` 같은 주소가
   바로 그 웹 앱입니다. 이 주소를 동료들에게 공유하면 **누구든 브라우저만 열면**
   사용할 수 있습니다 (설치 0건).

> 저장소를 Private으로 만들면 로그인한 팀원만 코드를 볼 수 있고, 앱 자체의 접근 권한도
> Streamlit Cloud 앱 설정(Settings → Sharing)에서 이메일 단위로 제한할 수 있습니다.
> 서열 정보가 민감하다면 이 옵션을 권장합니다.

## 로컬 컴퓨터(사내 서버 등)에 직접 띄우고 싶다면 (Docker 사용)

Docker를 쓸 수 있는 사내 서버가 있다면, 굳이 외부 호스팅 없이 그 서버에서
직접 돌릴 수도 있습니다. 함께 넣어둔 `Dockerfile`을 사용하세요:

```bash
docker build -t cdr-tool .
docker run -p 8501:8501 cdr-tool
```

그 후 브라우저에서 `http://localhost:8501` (같은 네트워크의 다른 PC에서는
`http://해당PC의IP:8501`) 접속.

## 사용법 (동료들이 실제로 하는 것)

1. 웹 주소 접속
2. 왼쪽에서 numbering scheme 선택 (기본값 kabat이면 됨)
3. 기존과 같은 형식의 입력 파일(.txt)을 업로드하거나, 텍스트를 직접 붙여넣기
4. **"분석 실행"** 버튼 클릭
5. 결과 확인 후 **엑셀 다운로드** 버튼으로 결과 파일(.xlsx) 저장

입력 형식:
- 블록(빈 줄로 구분) 순서: ① WT 아미노산 → ② WT nucleotide → ③ mutant들
- 같은 블록 안에서 LC/HC는 그냥 줄바꿈(Enter)으로 구분하면 됩니다.
- `>wild fasta` / `>wild nucleotide` / `>Mutant-01` 같은 제목 줄은 **완전히 생략 가능**합니다.
  생략하면 순서대로 자동 인식되고, 쓰고 싶으면 기존처럼 `>`로 시작하는 줄을 블록 맨 앞에
  추가하면 됩니다 (제목 있는 블록과 없는 블록을 섞어 써도 됩니다).

예시 (제목 없이):
```
DIQMTQSPSS...   (LC, wild type 아미노산)
EVQLVESG...     (HC)

GACATCC...      (LC, wild type nucleotide)
GAGGTGC...      (HC)

DIQMTQSPSS...   (mutant 1)
EVQLVESG...

DIQMTQSPSS...   (mutant 2)
EVQLVESG...
```

## 파일 구성

- `app.py` — 웹 UI (Streamlit)
- `cdr_codon_tool.py`, `excel_export.py`, `codon_table.py` — 기존 로직 그대로
- `requirements.txt` — 파이썬 패키지 목록 (Streamlit Cloud가 자동 설치)
- `packages.txt` — hmmer 등 시스템 패키지 목록 (Streamlit Cloud가 자동 설치)
- `Dockerfile` — 사내 서버 등에 직접 띄우고 싶을 때만 사용 (선택 사항)
