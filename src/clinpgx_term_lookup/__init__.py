"""clinpgx-term-lookup: fuzzy term lookup for ClinPGx/PharmGKB drugs and variants."""

from .drug_search import DrugLookup, DrugSearchResult
from .variant_search import VariantLookup, VariantSearchResult

__version__ = "0.1.0"

__all__ = [
    "DrugLookup",
    "DrugSearchResult",
    "VariantLookup",
    "VariantSearchResult",
    "__version__",
]
