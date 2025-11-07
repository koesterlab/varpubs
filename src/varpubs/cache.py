from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
import logging
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from sqlmodel import select

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


@dataclass
class Cache:
    path: Path
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        SQLModel.metadata.create_all(self.engine)
        if not self.path.exists():
            logger.info("New cache created.")
        else:
            logger.info("Existing cache found. Tables checked or updated.")

    def merge(self, other: "Cache", overwrite: bool = False) -> None:
        """Merge another CacheDB into this one."""
        with Session(self.engine) as session, Session(other.engine) as other_session:
            for summary in other_session.exec(select(Summary)):
                if (
                    overwrite
                    or not session.exec(
                        select(Summary).filter_by(
                            term=summary.term, pmid=summary.pmid, model=summary.model
                        )
                    ).first()
                ):
                    session.add(summary)
            session.commit()

    def lookup(self, term: str, pmid: int, model: str) -> Optional[Summary]:
        """Look up a summary by term and PMID."""
        with Session(self.engine) as session:
            return session.exec(
                select(Summary).filter_by(term=term, pmid=pmid, model=model)
            ).first()

    def write(self, summaries: List[Summary]) -> None:
        with Session(self.engine) as session:
            session.add_all(summaries)
            session.commit()

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        """Lazy-create SQLAlchemy engine for DuckDB."""
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
