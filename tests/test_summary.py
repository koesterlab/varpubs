from varpubs.pubmed_db import PubmedArticle
from varpubs.summarize import PubmedSummarizer, Settings
import os


def test_summarization():
    article = PubmedArticle(
        pmid=12345678,
        title="Genetic variants in BRCA1 and breast cancer risk",
        abstract="BRCA1 mutations significantly increase the risk of breast and ovarian cancers. In this study, we evaluate the frequency and pathogenicity of BRCA1 variants in a European cohort consisting of 1000 participants. Out of these, 777 developed breast cancer and 123 developed ovarian cancer within 2 years.",
        authors="Smith J, Doe A",
        journal="Genetics Today",
        pub_date="2023-09-15",
        doi="10.1234/genetictoday.2023.456",
    )

    settings = Settings(
        api_key=os.environ["LLM_API_KEY"], base_url=os.environ["LLM_API_BASE_URL"]
    )

    summarizer = PubmedSummarizer(settings)
    summary = summarizer.summarize(article)
    print(summary)
