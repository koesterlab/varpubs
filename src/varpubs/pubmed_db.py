from dataclasses import dataclass, field
from pathlib import Path
import sqlalchemy
from sqlmodel import Field, SQLModel


class PubmedArticle(SQLModel, table=True):
    pmid: int = Field(primary_key=True)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str


@dataclass
class PubmedDB:
    path: Path
    _engine: sqlalchemy.engine.base.Engine | None = field(init=False)

    def deploy(self) -> None:
        # wipe existing database
        if self.path.exists():
            self.path.unlink()
        # create parent directories if they don't exist
        if self.path.parent != Path():
            self.path.parent.mkdir(parents=True, exist_ok=True)
        # create the database
        SQLModel.metadata.create_all(self.engine)

        # load the data and insert it into the database
        # with Session(self.engine) as session:
        #     ...

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
