[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/koesterlab/varpubs/main.yml?branch=main&label=tests)](https://github.com/koesterlab/varpubs/actions)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/varpubs/badges/version.svg)](https://anaconda.org/bioconda/varpubs)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/varpubs/badges/latest_release_date.svg)](https://anaconda.org/bioconda/varpubs)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/varpubs/badges/license.svg)](https://anaconda.org/bioconda/varpubs)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/varpubs/badges/downloads.svg)](https://anaconda.org/bioconda/varpubs)

# varpubs

**varpubs** is a command-line tool that automates retrieval and summarization of variant-associated biomedical literature to support faster variant interpretation in oncology.

Identifying treatment-relevant evidence for rare or poorly characterized mutations often requires reviewing large volumes of PubMed articles. *varpubs* streamlines this process by querying literature for variant-linked publications and generating concise, role-specific three-sentence summaries using a large language model (LLM). Optional LLM-based “judges” can score articles based on user given terms (e.g., therapy relevance) to prioritize relevant articles.

Variant-level summaries are written to a TSV file for direct use in clinical dashboards, molecular tumor boards, or downstream pipelines.
A caching system stores previously generated summaries to reduce runtime in subsequent analyses.

## Installation

### Pixi
```
pixi global install varpubs
```

### Bioconda
```
conda install -c bioconda varpubs
```

## Usage

### 1. Deploy the PubMed variant database
```bash
varpubs deploy-db \
  --db-path pubmed.duckdb \
  --vcf-paths variants.vcf.gz other.vcf.gz \
  --email your_email@example.com
```

### 2. Summarize variant-associated articles
```bash
varpubs summarize-variants \
  --db-path pubmed.duckdb \
  --vcf-path variants_annotated.vcf.gz \
  --llm-url https://your-llm-endpoint \
  --api-key $HF_TOKEN \
  --model medgemma-27b-it \
  --role physician \
  --cache cache.duckdb \
  --output summaries.tsv \
  --output-cache tmp_cache.duckdb
```

### 3. Merge or update caches
```bash
varpubs update-cache \
  --cache tmp_cache.duckdb \
  --output cache.duckdb
```

## Authors
- Felix Wiegand
- Felix Mölder
- İlmay Taş
- Johannes Köster
