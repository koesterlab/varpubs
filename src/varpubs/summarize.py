from dataclasses import dataclass
from typing import Optional
from openai import OpenAI
from varpubs.pubmed_db import PubmedArticle


@dataclass
class Settings:
    api_key: str
    model: str = "teuken-7b-instruct-research"
    base_url: Optional[str] = None
    max_new_tokens: int = 500
    temperature: float = 0.1


@dataclass
class PubmedSummarizer:
    settings: Settings

    @property
    def client(self) -> OpenAI:
        return OpenAI(api_key=self.settings.api_key, base_url=self.settings.base_url)

    def summarize(self, article: PubmedArticle, term: str) -> str:
        instruction_text = "You are an oncologist."
        input_text = f"Concisely summarize the information in the following text regarding {term} in at most three sentences for your colleagues:\n\nTitle: {article.title}\n\n{article.abstract}"

        response = self.client.chat.completions.create(
            model=self.settings.model,
            messages=[
                {"role": "system", "content": instruction_text},
                {"role": "user", "content": input_text},
            ],
            temperature=self.settings.temperature,
            max_tokens=self.settings.max_new_tokens,
        )
        return str(response.choices[0].message.content)
