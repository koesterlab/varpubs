# src/varpubs/summarize.py
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

    def instruction(self) -> str:
        return f"You are an {self.settings.role}."

    def summarize(self, texts: list[tuple[PubmedArticle, str]], term: str) -> str:
        summaries = "\n".join(
            f"{article.pmid}: {summary}" for article, summary in texts
        )
        input_text = (
            f"Take the following summaries of PubMed articles related to {term}"
            "and summarize them in a single paragraph while ensuring that the summary is concise and informative."
            "When using information from the summaries, cite the PMID of the article."
            f"{summaries}"
        )
        response = self.client.chat.completions.create(
            model=self.settings.model,
            messages=[
                {"role": "system", "content": self.instruction()},
                {"role": "user", "content": input_text},
            ],
            max_tokens=self.settings.max_new_tokens,
            temperature=self.settings.temperature,
        )
        return str(response.choices[0].message.content)

    def summarize_article(self, article: PubmedArticle, term: str) -> str:
        if not article.abstract or not article.abstract.strip():
            return "No abstract available."

        input_text = (
            f"Concisely summarize the information in the following text regarding {term} "
            f"in at most three sentences for your colleagues:\n\n"
            f"Title: {article.title}\n\n{article.abstract}"
        )

        response = self.client.chat.completions.create(
            model=self.settings.model,
            messages=[
                {"role": "system", "content": self.instruction()},
                {"role": "user", "content": input_text},
            ],
            temperature=self.settings.temperature,
            max_tokens=self.settings.max_new_tokens,
        )
        return str(response.choices[0].message.content)

    def validate_summary(self, abstract: str, summary: str) -> bool:
        few_shots = [
            {
                "abstract": (
                    "While BRAF inhibitor combinations with EGFR and/or MEK inhibitors have improved clinical "
                    "efficacy in BRAFV600E colorectal cancer (CRC), response rates remain low and lack durability. "
                    "Preclinical data suggest that BRAF/MAPK pathway inhibition may augment the tumor immune response. "
                    "We performed a proof-of-concept single-arm phase 2 clinical trial of combined PD-1, BRAF and MEK "
                    "inhibition with spartalizumab (PDR001), dabrafenib and trametinib in 37 patients with BRAFV600E CRC. "
                    "The primary end point was overall response rate, and the secondary end points were progression-free "
                    "survival, disease control rate, duration of response and overall survival. The study met its primary "
                    "end point with a confirmed response rate (24.3% in all patients; 25% in microsatellite-stable "
                    "patients) and durability that were favorable relative to historical controls of BRAF-targeted "
                    "combinations alone. Single-cell RNA sequencing of 23 paired pretreatment and day 15 on-treatment "
                    "tumor biopsies revealed greater induction of tumor cell-intrinsic immune programs and more "
                    "complete MAPK inhibition in patients with better clinical outcome. Immune program induction in "
                    "matched patient-derived organoids correlated with the degree of MAPK inhibition. These data suggest "
                    "a potential tumor cell-intrinsic mechanism of cooperativity between MAPK inhibition and immune "
                    "response, warranting further clinical evaluation of optimized targeted and immune combinations in CRC."
                ),
                "summary": (
                    "In BRAFV600E colorectal cancer, combining PD-1, BRAF, and MEK inhibitors demonstrated a 24.3% overall "
                    "response rate and improved durability compared to historical BRAF-targeted therapies. Single-cell RNA "
                    "sequencing and organoid models revealed that better responses were linked to stronger MAPK inhibition "
                    "and induction of tumor-intrinsic immune programs. These findings suggest a cooperative mechanism "
                    "between MAPK suppression and immune activation in BRAFV600E CRC, supporting further clinical evaluation."
                ),
                "reason": "The study's findings support this claim.",
                "validation": True,
            },
            {
                "abstract": (
                    "Purpose: BEACON CRC evaluated encorafenib plus cetuximab with or without binimetinib versus "
                    "investigators' choice of irinotecan or FOLFIRI plus cetuximab in patients with BRAFV600E-mutant "
                    "metastatic colorectal cancer (mCRC), after progression on 1-2 prior regimens. In the previously "
                    "reported primary analysis, encorafenib, binimetinib plus cetuximab (ENCO/BINI/CETUX; triplet) and "
                    "encorafenib plus cetuximab (ENCO/CETUX; doublet) regimens improved overall survival (OS) and objective "
                    "response rate (ORR) versus standard of care. The purpose of this analysis was to report updated "
                    "efficacy and safety data.\n\nMethods: In this open-label, phase III trial, 665 patients with "
                    "BRAF V600E-mutant mCRC were randomly assigned 1:1:1 to receive triplet, doublet, or control. Primary "
                    "end points were OS and independently reviewed ORR comparing triplet to control. OS for doublet versus "
                    "control was a key secondary end point. Updated analyses include 6 months of additional follow-up and "
                    "ORR for all randomized patients.\n\nResults: Patients received triplet (n = 224), doublet (n = 220), "
                    "or control (n = 221). Median OS was 9.3 months for triplet and 5.9 months for control. Median OS for "
                    "doublet was 9.3 months. Confirmed ORR was 26.8 % for triplet, 19.5 % for doublet, and 1.8 % for control. "
                    "Adverse events were consistent with the prior primary analysis.\n\nConclusion: In the BEACON CRC study, "
                    "encorafenib plus cetuximab improved OS, ORR, and progression-free survival compared with standard "
                    "chemotherapy. Based on the analyses, encorafenib plus cetuximab is a new standard care regimen for "
                    "previously treated patients with BRAF V600E mCRC."
                ),
                "summary": (
                    "In 723 (524 male, 199 female) patients with BRAF V600E-mutant metastatic colorectal cancer, the BEACON "
                    "CRC phase III trial showed that encorafenib plus cetuximab, with or without binimetinib, significantly "
                    "improved overall survival and objective response rate compared to standard radiotherapy. Median overall "
                    "survival was approximately 9.3 months for both regimens versus 5.9 months for control."
                ),
                "reason": (
                    "The original abstract did not provide the number of patients nor mention radiotherapy; therefore the "
                    "summary contains unsubstantiated details."
                ),
                "validation": False,
            },
        ]

        few_shot_text = "\n\n".join(
            f"Abstract: {f['abstract']}\n\nSummary: {f['summary']}\n\nReason: {f['reason']}\n\nValidation: {f['validation']}"
            for f in few_shots
        )

        instruction_text = (
            "You are a scientific reviewer. "
            "Your task is to check whether the summary is factually accurate based ONLY on the abstract. "
            "Return exactly one word: 'True' or 'False'.\n\n"
            f"Examples:\n\n{few_shot_text}\n\n"
            "If any detail is incorrect, misleading, or cannot be confirmed, respond 'False'. "
            "Otherwise, respond 'True'."
        )

        user_prompt = (
            f"Abstract:\n{abstract}\n\n"
            f"Summary:\n{summary}\n\n"
            "Is the summary factually accurate based on the abstract?"
        )

        response = self.client.chat.completions.create(
            model=self.settings.model,
            messages=[
                {"role": "system", "content": instruction_text},
                {"role": "user", "content": user_prompt},
            ],
            temperature=1,
            max_tokens=5,
        )
        content = response.choices[0].message.content or ""
        answer = content.strip().lower()
        return answer == "true"
