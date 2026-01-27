import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import sqlalchemy
from sqlmodel import Field, Integer, Session, SQLModel, select

logger = logging.getLogger(__name__)


class Summary(SQLModel, table=True):
    """LLM summary cache of PubMed article stored locally."""

    term: str = Field(nullable=False, primary_key=True)
    pmid: int = Field(Integer, nullable=False, primary_key=True)
    model: str = Field(nullable=False)
    summary: str = Field(nullable=False)
    prompt_hash: str = Field(nullable=False)


class Judge(SQLModel, table=True):
    """LLM judge score of PubMed article stored locally."""

    term: str = Field(nullable=False, primary_key=True)
    pmid: int = Field(Integer, nullable=False, primary_key=True)
    model: str = Field(nullable=False)
    judge: str = Field(nullable=False)
    score: int = Field(nullable=False)
    prompt_hash: str = Field(nullable=False)


@dataclass
class Cache:
    path: Path
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        SQLModel.metadata.create_all(self.engine)
        if not self.path.exists():
            logger.info(f"New cache created under {self.path}")
        else:
            logger.info(f"Using existing cache {self.path}")

    def merge(self, other: "Cache", overwrite: bool = False) -> None:
        """Merge another CacheDB into this one."""
        with Session(self.engine) as session, Session(other.engine) as other_session:
            for summary in other_session.exec(select(Summary)):
                if (
                    overwrite
                    or not session.exec(
                        select(Summary).filter_by(
                            term=summary.term,
                            pmid=summary.pmid,
                            model=summary.model,
                            prompt_hash=summary.prompt_hash,
                        )
                    ).first()
                ):
                    session.add(Summary(**summary.model_dump()))
            for judge in other_session.exec(select(Judge)):
                if (
                    overwrite
                    or not session.exec(
                        select(Judge).filter_by(
                            term=judge.term,
                            pmid=judge.pmid,
                            model=judge.model,
                            judge=judge.judge,
                            prompt_hash=judge.prompt_hash,
                        )
                    ).first()
                ):
                    session.add(Judge(**judge.model_dump()))
            session.commit()

    def lookup_summary(
        self, term: str, pmid: int, model: str, prompt_hash: str
    ) -> Optional[Summary]:
        """Look up a summary by term and PMID."""
        with Session(self.engine) as session:
            entry = session.exec(
                select(Summary).filter_by(
                    term=term, pmid=pmid, model=model, prompt_hash=prompt_hash
                )
            ).first()
            if entry:
                logger.info(f"Found summary cache entry for {term} (pmid: {pmid})")
                return entry
            else:
                return None

    def write_summaries(self, summaries: List[Summary]) -> None:
        with Session(self.engine) as session:
            summaries = [
                summary
                for summary in summaries
                if self.lookup_summary(
                    summary.term, summary.pmid, summary.model, summary.prompt_hash
                )
            ]
            session.add_all(summaries)
            session.commit()

    def lookup_judge(
        self, term: str, pmid: int, model: str, judge: str, prompt_hash: str
    ) -> Optional[int]:
        """Look up a judge by term, PMID used model and prompt hash."""
        with Session(self.engine) as session:
            j = session.exec(
                select(Judge).filter_by(
                    term=term,
                    pmid=pmid,
                    model=model,
                    judge=judge,
                    prompt_hash=prompt_hash,
                )
            ).first()
            if j:
                logger.info(
                    f"Found judgement cache entry for {term} with judge term {judge} and PMID {pmid}: {j.score}"
                )
                return j.score
            else:
                return None

    def write_judges(self, judges: List[Judge]) -> None:
        with Session(self.engine) as session:
            judges = [
                judge
                for judge in judges
                if self.lookup_judge(
                    judge.term, judge.pmid, judge.model, judge.judge, judge.prompt_hash
                )
            ]
            session.add_all(judges)
            session.commit()

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        """Lazy-create SQLAlchemy engine for DuckDB."""
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
