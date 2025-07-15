from dataclasses import dataclass
from typing import Optional
from huggingface_hub import InferenceClient
from varpubs.pubmed_db import PubmedArticle


@dataclass
class HFSettings:
    token: str
    model: str = "mistralai/Mistral-7B-Instruct-v0.1"
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
        prompt = (
            f"Summarize the following PubMed abstract about genetic variants:\n\n"
            f"Title: {article.title}\n\n"
            f"Abstract: {article.abstract}\n"
        )
        result = self.client.text_generation(
            prompt,
            max_new_tokens=self.settings.max_new_tokens,
            temperature=self.settings.temperature,
        )
        return result.strip()
