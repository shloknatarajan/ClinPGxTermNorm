from pydantic import BaseModel
from typing import List, Optional
from clinpgx_lookup.search_utils import (
    general_search,
    general_search_comma_list,
)
from clinpgx_lookup.data import get_tsv_path
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class ChemicalSearchResult(BaseModel):
    raw_input: str
    id: str
    name: str
    url: str
    score: float
    source: str = "clinpgx"


class ChemicalLookup:
    def __init__(self, data_path: Optional[str] = None) -> None:
        self.data_path = data_path or str(get_tsv_path("chemicals"))
        self._df: Optional[pd.DataFrame] = None

    @property
    def df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self.data_path, sep="\t")
        return self._df

    def search(
        self, chemical: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[ChemicalSearchResult]:
        """Search for a chemical by name, generic names, or trade names."""
        results = general_search(
            self.df,
            chemical,
            "Name",
            "PharmGKB Accession Id",
            threshold=threshold,
            top_k=top_k,
        )
        if results and results[0]["score"] >= threshold:
            return self._to_results(chemical, results[:top_k])

        alt_results = general_search_comma_list(
            self.df,
            chemical,
            "Generic Names",
            "PharmGKB Accession Id",
            threshold=threshold,
            top_k=top_k,
        )
        alt_results.extend(
            general_search_comma_list(
                self.df,
                chemical,
                "Trade Names",
                "PharmGKB Accession Id",
                threshold=threshold,
                top_k=top_k,
            )
        )
        all_results = results + alt_results
        all_results.sort(key=lambda x: x["score"], reverse=True)
        return self._to_results(chemical, all_results[:top_k])

    def _to_results(self, raw_input: str, results: list) -> List[ChemicalSearchResult]:
        return [
            ChemicalSearchResult(
                raw_input=raw_input,
                id=r["PharmGKB Accession Id"],
                name=r.get("Name", ""),
                url=f"https://www.clinpgx.org/chemical/{r['PharmGKB Accession Id']}",
                score=r["score"],
            )
            for r in results
        ]
