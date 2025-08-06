from dataclasses import dataclass
from typing import Optional
from openai import OpenAI
from varpubs.pubmed_db import PubmedArticle


@dataclass
class Settings:
    api_key: str
    role: str
    model: str = "medgemma-27b-it"
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

        if not article.abstract.strip():
            return "No abstract available."

        instruction_text = f"You are an {self.settings.role}."
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

    def validate_summary(self, abstract: str, summary: str) -> bool:
        few_shots = [
            {
                "abstract": "While BRAF inhibitor combinations with EGFR and/or MEK inhibitors have improved clinical efficacy in BRAFV600E colorectal cancer (CRC), response rates remain low and lack durability. ...",
                "summary": "In BRAFV600E colorectal cancer, combining PD-1, BRAF, and MEK inhibitors demonstrated a 24.3% overall response rate and improved durability compared to historical BRAF-targeted therapies...",
                "reason": "The study's findings support this claim.",
                "validation": True,
            },
            {
                "abstract": "Purpose: BEACON CRC evaluated encorafenib plus cetuximab with or without binimetinib versus investigators' choice of irinotecan or FOLFIRI plus cetuximab in patients with BRAFV600E-mutant metastatic colorectal cancer...",
                "summary": "In 723 (524 male, 199 female) patients with BRAF V600E-mutant metastatic colorectal cancer, the BEACON CRC phase III trial showed that encorafenib plus cetuximab significantly improved overall survival...",
                "reason": "The original abstract did not provide the number of patients nor mention radiotherapy.",
                "validation": False,
            },
        ]
        few_shot_text = [
            f"Abstract: {f['abstract']}\n\nSummary: {f['summary']}\n\nReason: {f['reason']}\n\nValidation: {f['validation']}"
            for f in few_shots
        ]
        instruction_text = (
            "You are a scientific reviewer. "
            "Your task is to check whether the summary is factually accurate based ONLY on the abstract. "
            f"Examples:\n\n{few_shot_text}\n\n"
            "If any detail is incorrect, misleading, or cannot be confirmed, respond with 'False'. "
            "Otherwise, respond with 'True'. Respond with exactly one word: 'True' or 'False'."
        )

        input_text = (
            f"Abstract:\n{abstract}\n\n"
            f"Summary:\n{summary}\n\n"
            f"Is the summary factually accurate based on the abstract?"
        )

        response = self.client.chat.completions.create(
            model=self.settings.model,
            messages=[
                {"role": "system", "content": instruction_text},
                {"role": "user", "content": input_text},
            ],
            temperature=1,
            max_tokens=5,
        )

        answer = str(response.choices[0].message.content).strip().lower()
        return "true" in answer
