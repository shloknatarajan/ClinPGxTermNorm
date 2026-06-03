"""Command-line entry point for clinpgx-term-lookup."""

import argparse
import json

from .drug_search import DrugLookup
from .variant_search import VariantLookup


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="clinpgx-term-lookup",
        description="Look up a drug or variant term in the ClinPGx/PharmGKB databases.",
    )
    parser.add_argument(
        "query", help="The term to look up (e.g. 'warfarin', 'rs1234')."
    )
    parser.add_argument(
        "-t",
        "--type",
        choices=["drug", "variant"],
        required=True,
        help="Whether the term is a drug or a variant.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.8,
        help="Fuzzy-match threshold (default: 0.8).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=1,
        help="Maximum number of results to return (default: 1).",
    )
    args = parser.parse_args()

    if args.type == "drug":
        results = DrugLookup().search(
            args.query, threshold=args.threshold, top_k=args.top_k
        )
    else:
        results = VariantLookup().search(
            args.query, threshold=args.threshold, top_k=args.top_k
        )

    print(json.dumps([result.model_dump() for result in (results or [])], indent=2))


if __name__ == "__main__":
    main()
