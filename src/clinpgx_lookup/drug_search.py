from pydantic import BaseModel
from typing import List, Optional
import requests
from clinpgx_lookup.search_utils import (
    calc_similarity,
    general_search,
    general_search_comma_list,
)
from clinpgx_lookup.data import get_tsv_path
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class DrugSearchResult(BaseModel):
    raw_input: str
    id: str
    name: str
    url: str
    score: float
    source: str = "clinpgx"


def get_first_rxnorm_candidate(data):
    """Get the first candidate object with RXNORM as the source."""
    candidates = data.get("approximateGroup", {}).get("candidate", [])
    for candidate in candidates:
        if candidate.get("source") == "RXNORM":
            return candidate
    return None


def rxnorm_search(drug_name: str) -> Optional[DrugSearchResult]:
    url = "https://rxnav.nlm.nih.gov/REST/approximateTerm.json"
    params = {"term": drug_name, "maxEntries": 1}
    try:
        response = requests.get(url, params=params, timeout=10)
    except requests.RequestException:
        logger.warning("RxNorm API request failed for '%s'", drug_name)
        return None
    if response.status_code == 200:
        data = response.json()
        candidate = get_first_rxnorm_candidate(data)
        if candidate:
            rxcui = candidate["rxcui"]
            result_url = f"https://ndclist.com/rxnorm/rxcui/{rxcui}"
            name = candidate["name"]
            score = calc_similarity(drug_name, name)
            return DrugSearchResult(
                raw_input=drug_name,
                id=f"RXN{rxcui}",
                name=name,
                url=result_url,
                score=score,
                source="rxnorm",
            )
    return None


class DrugLookup:
    def __init__(self, data_path: Optional[str] = None) -> None:
        self.data_path = data_path or str(get_tsv_path("drugs"))
        self._df: Optional[pd.DataFrame] = None

    @property
    def df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self.data_path, sep="\t")
        return self._df

    def _clinpgx_drug_name_search(
        self, drug_name: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[DrugSearchResult]:
        results = general_search(
            self.df,
            drug_name,
            "Name",
            "PharmGKB Accession Id",
            threshold=threshold,
            top_k=top_k,
        )
        if results:
            return [
                DrugSearchResult(
                    raw_input=drug_name,
                    id=result["PharmGKB Accession Id"],
                    name=result["Name"],
                    url=f"https://www.clinpgx.org/chemical/{result['PharmGKB Accession Id']}",
                    score=result["score"],
                )
                for result in results
            ]
        return []

    def _clinpgx_drug_alternatives_search(
        self, drug_name: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[DrugSearchResult]:
        """Checks generic names and trade names for the drug."""
        results = general_search_comma_list(
            self.df,
            drug_name,
            "Generic Names",
            "PharmGKB Accession Id",
            threshold=threshold,
            top_k=top_k,
        )
        results.extend(
            general_search_comma_list(
                self.df,
                drug_name,
                "Trade Names",
                "PharmGKB Accession Id",
                threshold=threshold,
                top_k=top_k,
            )
        )
        if results:
            return [
                DrugSearchResult(
                    raw_input=drug_name,
                    id=result["PharmGKB Accession Id"],
                    name=result["Name"],
                    url=f"https://www.clinpgx.org/chemical/{result['PharmGKB Accession Id']}",
                    score=result["score"],
                )
                for result in results
            ]
        return []

    def clinpgx_lookup(
        self, drug_name: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[DrugSearchResult]:
        """Try name search first, then alternatives if scores are too low."""
        name_results = self._clinpgx_drug_name_search(
            drug_name, threshold=threshold, top_k=top_k
        )
        if name_results and any(r.score >= threshold for r in name_results):
            return name_results

        alternatives_results = self._clinpgx_drug_alternatives_search(
            drug_name, threshold=threshold, top_k=top_k
        )
        all_results = (name_results or []) + (alternatives_results or [])
        if all_results:
            all_results.sort(key=lambda x: x.score, reverse=True)
            return all_results[:top_k]
        return []

    def rxnorm_lookup(self, drug_name: str) -> List[DrugSearchResult]:
        """Search using RxNorm and convert results back to PharmGKB format."""
        rxnorm_result = rxnorm_search(drug_name)
        if not rxnorm_result or not rxnorm_result.id:
            return []

        rxcui = rxnorm_result.id.removeprefix("RXN")
        # Try to map RxCUI back to PharmGKB
        rxcui_results = general_search(
            self.df,
            rxcui,
            "RxNorm Identifiers",
            "PharmGKB Accession Id",
            threshold=0.8,
            top_k=1,
        )
        if rxcui_results:
            return [
                DrugSearchResult(
                    raw_input=drug_name,
                    id=r["PharmGKB Accession Id"],
                    name=r["Name"],
                    url=f"https://www.clinpgx.org/chemical/{r['PharmGKB Accession Id']}",
                    score=r["score"],
                )
                for r in rxcui_results
            ]
        # If no PharmGKB mapping, return the RxNorm result directly
        return [rxnorm_result]

    def search(
        self, drug_name: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[DrugSearchResult]:
        """Search for a drug. Tries ClinPGx local data first, then RxNorm as fallback."""
        results = self.clinpgx_lookup(drug_name, threshold=threshold, top_k=top_k)
        if results:
            return results
        logger.info("No ClinPGx results for '%s', trying RxNorm", drug_name)
        return self.rxnorm_lookup(drug_name)
