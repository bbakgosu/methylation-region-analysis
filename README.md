# Methylation Analysis Pipeline

특정 genomic region의 전체 메틸화율(overall methylation percentage)을 계산하는 자동화 파이프라인입니다.

## 목적 (Why this tool?)

Oxford Nanopore 시퀀싱 데이터에서 메틸화를 분석할 때, modkit은 개별 CpG 위치별 메틸화율만 제공합니다.
이 도구는 다음과 같은 추가 기능을 제공합니다:

1. **Region-level 통합 메틸화율 계산**: 지정된 genomic region 전체의 평균 메틸화율을 자동으로 계산
2. **가독성 높은 출력**: 시각화 바 차트와 상세한 통계를 포함한 포맷팅된 결과
3. **원스톱 자동화**: BAM 파일에서 최종 분석까지 한 번의 명령으로 처리
4. **사용자 친화적**: 큰 숫자에 쉼표 지원 (예: chr1:1,000,000-2,000,000)
5. **로그 및 요약 파일**: 재현 가능한 분석을 위한 자동 문서화

## 주요 기능

- **자동화된 파이프라인**: modkit pileup → bedMethyl 생성 → 메틸화 통계 계산
- **5mC 메틸화 분석**: 시토신 메틸화(5mC)에 특화
- **상세한 통계**:
  - Region 전체의 메틸화 백분율
  - CpG 위치 수, 총 리드 수
  - Coverage 통계 (평균, 중앙값)
  - Strand 분포 (plus/minus)
  - 시각화 바 차트

## 시스템 요구사항

### 필수 소프트웨어

- **modkit**: Oxford Nanopore의 메틸화 분석 도구
- **Python 3.7+**
- **pandas**: Python 데이터 분석 라이브러리

### 입력 파일 요구사항

- BAM 파일에 메틸화 태그(Mm, Ml)가 포함되어 있어야 함
- BAM 인덱스 파일(.bai)이 권장됨 (빠른 처리를 위해)

## 설치 방법

### 1. 저장소 다운로드

먼저 이 저장소를 다운로드합니다:

```bash
# Git으로 클론
git clone https://github.com/bbakgosu/methylation-region-analysis.git

# 디렉토리 이동
cd methylation-region-analysis
```

### 2. 옵션 1: Conda 환경 사용 (권장)

모든 의존성을 포함한 완전한 환경을 한 번에 설치:

```bash
# 환경 생성 및 활성화
conda env create -f environment.yml
conda activate methylation-analysis

# 설치 확인
modkit --version
python3 count_methylation_percent.py -h
```

### 2. 옵션 2: 개별 설치

#### modkit 설치

```bash
# Conda로 설치 (권장)
conda install -c bioconda modkit

# 또는 직접 다운로드
# https://github.com/nanoporetech/modkit/releases
```

#### Python 패키지 설치

```bash
pip install -r requirements.txt
# 또는
pip install pandas
```

### 3. 실행 권한 부여

```bash
chmod +x methylation_analysis.sh
chmod +x count_methylation_percent.py
```

## 사용 방법

### 기본 사용법

```bash
./methylation_analysis.sh -b <BAM_FILE> -r <REGION>
```

### 필수 인자

- `-b, --bam FILE`: 메틸화 태그를 포함한 BAM 파일
- `-r, --region REGION`: 분석할 genomic region (형식: `chr:start-end`)

### 선택 인자

- `-o, --output DIR`: 출력 디렉토리 (기본값: `./output`)
- `-t, --threads NUM`: modkit에 사용할 스레드 수 (기본값: 8)
- `-s, --script PATH`: Python 스크립트 경로 (기본값: `./count_methylation_percent.py`)
- `-m, --modkit PATH`: modkit 실행 파일 경로 (기본값: `modkit`)
- `--simple`: 간단한 출력 형식 사용
- `-v, --verbose`: 상세한 로그 출력
- `-h, --help`: 도움말 표시

## 사용 예제

### 예제 1: 기본 사용

```bash
./methylation_analysis.sh \
  -b sample.bam \
  -r chr1:14923-15923
```

### 예제 2: 큰 숫자에 쉼표 사용

```bash
./methylation_analysis.sh \
  -b sample.bam \
  -r "chr1:1,000,000-2,000,000"
```

### 예제 3: 사용자 정의 출력 디렉토리 및 스레드 수

```bash
./methylation_analysis.sh \
  -b sample.bam \
  -r chr1:14923-15923 \
  -o results/chr1_analysis \
  -t 16
```

### 예제 4: 상세 로그와 간단한 출력

```bash
./methylation_analysis.sh \
  -b sample.bam \
  -r chr1:14923-15923 \
  -v --simple
```

### 예제 5: Python 스크립트만 사용 (bedMethyl 파일이 이미 있는 경우)

```bash
python3 count_methylation_percent.py \
  -i sample_bedmethyl.bed \
  -r chr1:14923-15923
```

## 출력 파일

파이프라인은 다음 파일들을 생성합니다:

### 1. bedMethyl 파일 (`*_bedmethyl.bed`)

modkit이 생성한 원본 메틸화 데이터 (BED 형식)
- 각 CpG 위치별 메틸화 정보
- Coverage, 메틸화 비율, 리드 카운트 등 포함

### 2. 요약 파일 (`*_methylation_summary.txt`)

분석 결과 요약:
```
Sample: sample
BAM file: sample.bam
Region: chr1:14923-15923
Analysis Date: 2024-01-04

============================================================
Methylation Analysis for Region: chr1:14923-15923
============================================================

📊 Overall Statistics:
  • Total CpG positions analyzed: 150
  • Total reads analyzed: 25,432
  • Mean coverage per position: 169.5
  • Median coverage per position: 165.0

🧬 Methylation Counts:
  • Methylated reads: 18,574
  • Unmethylated reads: 6,858

📈 Methylation Percentage:
  • Overall methylation: 73.04%
  • Visual: [█████████████████████████████░░░░░░░░░░░] 73.0%
```

### 3. 로그 파일 (`*_methylation_analysis.log`)

실행 상세 로그:
- 실행 시간, 파라미터
- modkit 출력
- 에러 메시지 (있을 경우)

## 출력 해석

### 주요 메트릭 설명

- **Total CpG positions analyzed**: Region 내에서 분석된 CpG 위치의 수
- **Total reads analyzed**: 모든 위치에서 분석된 총 리드 수 (메틸화 + 비메틸화)
- **Mean coverage per position**: CpG 위치당 평균 coverage
- **Methylated reads**: 메틸화된 시토신을 포함한 리드 수
- **Unmethylated reads**: 메틸화되지 않은 시토신을 포함한 리드 수
- **Overall methylation**: Region 전체의 메틸화 비율 (%)
  - 계산식: (Methylated reads / Total reads) × 100

### 시각화 바 차트

```
[█████████████████████████████░░░░░░░░░░░] 73.0%
```

- `█` (채워진 블록): 메틸화 비율
- `░` (빈 블록): 비메틸화 비율
- 총 40칸으로 0-100% 표시

## 문제 해결

### modkit을 찾을 수 없음

```bash
Error: modkit is not installed or not in PATH
```

**해결 방법**:
```bash
# Conda로 설치
conda install -c bioconda modkit

# 또는 PATH에 추가
export PATH="/path/to/modkit:$PATH"

# 또는 전체 경로 지정
./methylation_analysis.sh -b sample.bam -r chr1:100-200 -m /full/path/to/modkit
```

### Python 패키지 오류

```bash
ModuleNotFoundError: No module named 'pandas'
```

**해결 방법**:
```bash
pip install pandas
# 또는
conda install pandas
```

### bedMethyl 파일이 비어있음

```bash
Warning: bedMethyl file is empty. Region might not contain methylation data.
```

**가능한 원인**:
1. BAM 파일에 메틸화 태그가 없음
2. 지정된 region에 데이터가 없음
3. Region 형식이 chromosome naming과 맞지 않음 (예: "chr1" vs "1")

**해결 방법**:
```bash
# BAM 파일의 메틸화 태그 확인
samtools view sample.bam | head -1

# Chromosome naming 확인
samtools view -H sample.bam | grep @SQ
```

### Region 형식 오류

```bash
Error: Invalid region format.
```

**올바른 형식**:
- `chr1:14923-15923` ✓
- `chr1:1,000,000-2,000,000` ✓ (쉼표 허용)
- `1:100-200` ✓ (chromosome 이름은 BAM 헤더와 일치해야 함)

**잘못된 형식**:
- `chr1-14923-15923` ✗ (콜론 대신 하이픈)
- `chr1:14923` ✗ (종료 위치 누락)

## 기술적 세부사항

### Coordinate 시스템

- **입력**: 1-based coordinates (사용자 친화적)
- **내부 처리**: 0-based coordinates (BED 표준)
- **출력**: 1-based coordinates (보고서)

### 메틸화 타입

현재 5mC (시토신 메틸화)만 분석합니다. bedMethyl 파일에서 `mod_type == 'm'`인 행만 사용.

### 계산 방식

1. **modkit pileup**: BAM 파일에서 지정된 region의 메틸화 정보 추출
2. **Filtering**: 5mC 타입만 선택
3. **Aggregation**: 모든 CpG 위치의 메틸화/비메틸화 리드 수 합산
4. **Percentage**: (총 메틸화 리드 / 총 리드) × 100
