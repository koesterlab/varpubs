from dataclasses import dataclass
from typing import Optional
from huggingface_hub import InferenceClient
from varpubs.pubmed_db import PubmedArticle


@dataclass
class HFSettings:
    token: str
    model: str = "facebook/bart-large-cnn"
    max_new_tokens: int = 200
    temperature: float = 0.1


@dataclass
class PubmedSummarizer:
    settings: HFSettings
    _client: Optional[InferenceClient] = None

    @property
    def client(self) -> InferenceClient:
        if self._client is None:
            self._client = InferenceClient(
                model=self.settings.model,
                token=self.settings.token,
            )
        return self._client

    def summarize(self, article: PubmedArticle) -> str:
        input_text = f"Title: {article.title}\n\n{article.abstract}"
        result = self.client.summarization(input_text)
        return result.summary_text.strip()
