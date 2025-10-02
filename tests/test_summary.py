from varpubs.pubmed_db import PubmedArticle
from varpubs.summarize import PubmedSummarizer, Settings
import os
import pytest

ARTICLE = article = PubmedArticle(
    pmid=12345678,
    title="Impact of KRASG12 mutations on survival with trifluridine/tipiracil plus bevacizumab in patients with refractory metastatic colorectal cancer: post hoc analysis of the phase III SUNLIGHT trial",
    abstract="""Background: In metastatic colorectal cancer (mCRC), KRAS mutations are often associated with poorer survival; however, the prognostic impact of specific point mutations is unclear. In the phase III SUNLIGHT trial, trifluridine/tipiracil (FTD/TPI) plus bevacizumab significantly improved overall survival (OS) versus FTD/TPI alone. We assessed the impact of KRASG12 mutational status on OS in SUNLIGHT.

    Patients and methods: In the global, open-label, randomized, phase III SUNLIGHT trial, adults with mCRC who had received no more than two prior chemotherapy regimens were randomized 1 : 1 to receive FTD/TPI alone or FTD/TPI plus bevacizumab. In this post hoc analysis, OS was assessed according to the presence or absence of a KRASG12 mutation in the overall population and in patients with RAS-mutated tumors.

    Results: Overall, 450 patients were analyzed, including 302 patients in the RAS mutation subgroup (214 with a KRASG12 mutation and 88 with a non-KRASG12RAS mutation). In the overall population, similar OS outcomes were observed in patients with and without a KRASG12 mutation [median 8.3 and 9.2 months, respectively; hazard ratio (HR) 1.09, 95% confidence interval (CI) 0.87-1.4]. Similar OS outcomes were also observed in the subgroup analysis of patients with a KRASG12 mutation versus those with a non-KRASG12RAS mutation (HR 1.03, 95% CI 0.76-1.4). FTD/TPI plus bevacizumab improved OS compared with FTD/TPI alone irrespective of KRASG12 mutational status. Among patients with a KRASG12 mutation, the median OS was 9.4 months with FTD/TPI plus bevacizumab versus 7.2 months with FTD/TPI alone (HR 0.67, 95% CI 0.48-0.93), and in patients without a KRASG12 mutation, the median OS was 11.3 versus 7.1 months, respectively (HR 0.59, 95% CI 0.43-0.81).

    Conclusions: The presence of a KRASG12 mutation had no detrimental effect on OS among patients treated in SUNLIGHT. The benefit of FTD/TPI plus bevacizumab over FTD/TPI alone was confirmed independently of KRASG12 status.""",
    authors="J Tabernero et al.",
    journal="Genetics Today",
    pub_date="2023-09-15",
    doi="10.1016/j.esmoop.2024.102945",
)

THERAPY_UNRELATED_ARTICLE = PubmedArticle(
    pmid=435346,
    title="KRAS Mutations Observed in Blue-Colored Zebrafish Larvae",
    abstract=(
        "We sequenced the genomes of blue-colored zebrafish larvae and identified several novel KRAS mutations. "
        "The study focused exclusively on the distribution of these mutations across developmental stages and pigmentation patterns. "
    ),
    authors="A. Nemo et al.",
    journal="Aquatic Genetics",
    pub_date="2023-09-15",
    doi="10.1016/j.esmoop.2024.102945",
)


def settings():
    return Settings(
        api_key=os.environ["LLM_API_KEY"],
        base_url=os.environ["LLM_API_BASE_URL"],
        role="oncologist",
    )


skip_if_no_api_key = pytest.mark.skipif(
    not os.environ.get("LLM_API_KEY"), reason="requires LLM_API_KEY"
)


@skip_if_no_api_key
def test_positive_judgment():
    summarizer = PubmedSummarizer(settings())
    score = summarizer.judge(ARTICLE, "therapy related")
    assert score > 1


@skip_if_no_api_key
def test_negative_judgment():
    summarizer = PubmedSummarizer(settings())
    score = summarizer.judge(THERAPY_UNRELATED_ARTICLE, "therapy unrelated")
    assert score < 2


@skip_if_no_api_key
def test_summarization():
    summarizer = PubmedSummarizer(settings())
    summary = summarizer.summarize_article(ARTICLE, "G12")
    assert summary


@skip_if_no_api_key
def test_negative_summary_validation():
    summarizer = PubmedSummarizer(settings())
    summary = "The paper validates the impact of rainbows per day (RPD) in western countries on cancer survival while eating ice cream."
    assert not summarizer.validate_summary(ARTICLE.abstract, summary)


@skip_if_no_api_key
def test_negative_summary_validation_hard():
    summarizer = PubmedSummarizer(settings())
    summary = "The SUNLIGHT trial showed that BRAF V600E mutations did impact OS in patients with refractory metastatic colorectal cancer treated with trifluridine/tipiracil plus bevacizumab."
    assert not summarizer.validate_summary(ARTICLE.abstract, summary)


@skip_if_no_api_key
def test_positive_summary_validation():
    summarizer = PubmedSummarizer(settings())
    summary = """The SUNLIGHT trial showed that KRASG12 mutations did not impact OS in patients with refractory metastatic colorectal cancer treated with trifluridine/tipiracil plus bevacizumab.
    The benefit of FTD/TPI plus bevacizumab over FTD/TPI alone was confirmed independently of KRASG12 mutational status."""
    assert summarizer.validate_summary(ARTICLE.abstract, summary)
