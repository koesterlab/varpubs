# Changelog

## [0.3.0](https://github.com/koesterlab/varpubs/compare/v0.2.2...v0.3.0) (2025-09-30)


### Features

* Add article judging and ranking by relevance score based on judgements ([#18](https://github.com/koesterlab/varpubs/issues/18)) ([6ad3154](https://github.com/koesterlab/varpubs/commit/6ad31541355e3f519ab95c39cb7676d9dd1b2991))

## [0.2.2](https://github.com/koesterlab/varpubs/compare/v0.2.1...v0.2.2) (2025-09-29)


### Bug Fixes

* Add batch_size parameter for PubMed queries ([#16](https://github.com/koesterlab/varpubs/issues/16)) ([e09880e](https://github.com/koesterlab/varpubs/commit/e09880eff2a40e7c43b65e7ed054ceb794202362))

## [0.2.1](https://github.com/koesterlab/varpubs/compare/v0.2.0...v0.2.1) (2025-09-10)


### Bug Fixes

* Split variant term into gene and hgvsp in output CSV and rename and reorder columns ([#14](https://github.com/koesterlab/varpubs/issues/14)) ([0ffbfd0](https://github.com/koesterlab/varpubs/commit/0ffbfd0ba7833128abd1e284a21366e6ab83bee3))

## [0.2.0](https://github.com/koesterlab/varpubs/compare/v0.1.2...v0.2.0) (2025-09-08)


### Features

* Summarize across all evidence for one specific variant of interest ([#13](https://github.com/koesterlab/varpubs/issues/13)) ([13743fd](https://github.com/koesterlab/varpubs/commit/13743fd1f11ccdad5d37d1384d503483d0a1370e))


### Bug Fixes

* Clean up term set while still searching pubmed for synonyms ([#12](https://github.com/koesterlab/varpubs/issues/12)) ([d79595f](https://github.com/koesterlab/varpubs/commit/d79595f6e307014874995f243c4d5da3f641a13c))
* Fix DOI extraction from ELocationID in PubmedDB ([#10](https://github.com/koesterlab/varpubs/issues/10)) ([2771fbb](https://github.com/koesterlab/varpubs/commit/2771fbbad9eca36b6ff700db9bd39a13ab3121c6))

## [0.1.2](https://github.com/koesterlab/varpubs/compare/v0.1.1...v0.1.2) (2025-08-27)


### Bug Fixes

* fixed annotation queries and set api-key default ([#8](https://github.com/koesterlab/varpubs/issues/8)) ([7a2ac48](https://github.com/koesterlab/varpubs/commit/7a2ac481ef19d81f78881277398c104ed68bd801))

## [0.1.1](https://github.com/koesterlab/varpubs/compare/v0.1.0...v0.1.1) (2025-08-27)


### Documentation

* Add MIT License to the project ([#6](https://github.com/koesterlab/varpubs/issues/6)) ([23da841](https://github.com/koesterlab/varpubs/commit/23da841e10fa4b21c1c85b86c3bb5f3db0fb6f33))

## 0.1.0 (2025-08-26)


### Features

* Add very basic summarize skeleton code ([#2](https://github.com/koesterlab/varpubs/issues/2)) ([24aba9c](https://github.com/koesterlab/varpubs/commit/24aba9ce98524079196027f9f5ad131072a5acc7))
* Parse HGVS.p from VCF and store related PubMed abstracts in DuckDB ([#3](https://github.com/koesterlab/varpubs/issues/3)) ([cf6f1ac](https://github.com/koesterlab/varpubs/commit/cf6f1aca96b25e081f375259d167a6c1f5128ba6))
