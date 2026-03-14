# ClinPGx Lookup

Normalize clinical pharmacogenomics terms to their canonical [ClinPGx](https://www.clinpgx.org/) database identifiers using fuzzy matching with API fallbacks.

Supports **drugs**, **variants** (rsIDs and star alleles), **genes**, **phenotypes**, and **chemicals**.

## Installation

```bash
pip install clinpgx-lookup
```

## Quick Start

```python
from clinpgx_lookup import normalize

# Drug lookup
results = normalize("aspirin", entity_type="drug")
print(results[0].name, results[0].id, results[0].score)
# aspirin PA448497 1.0

# Variant lookup (auto-detects rsID)
results = normalize("rs1045642")
print(results[0].name, results[0].id)

# Star allele lookup (auto-detected)
results = normalize("CYP2C19*2")

# Gene lookup
results = normalize("CYP2D6", entity_type="gene")

# Phenotype lookup
results = normalize("asthma", entity_type="phenotype")
```

### Using Lookup Classes Directly

For more control, use the individual lookup classes:

```python
from clinpgx_lookup import DrugLookup, VariantLookup, GeneLookup

drug_lookup = DrugLookup()
results = drug_lookup.search("Tylenol", threshold=0.7, top_k=3)

variant_lookup = VariantLookup()
results = variant_lookup.rsid_lookup("rs1000002")
results = variant_lookup.star_lookup("CYP2C19*2")

gene_lookup = GeneLookup()
results = gene_lookup.search("CYP2D6")
```

## How It Works

### Search Strategy

Each entity type uses a cascading search strategy:

**Drugs:**
1. Fuzzy match against drug Name
2. Fuzzy match against Generic Names and Trade Names
3. Fallback to [RxNorm API](https://rxnav.nlm.nih.gov/) with mapping back to PharmGKB IDs

**Variants:**
1. Query [PharmGKB API](https://api.pharmgkb.org/) (rsID or star allele endpoint)
2. Fuzzy match against local variant names and synonyms

**Genes:**
1. Fuzzy match against gene Symbol
2. Fuzzy match against Name, Alternate Names, and Alternate Symbols

**Phenotypes / Chemicals:**
1. Fuzzy match against Name
2. Fuzzy match against Alternate Names (or Generic/Trade Names for chemicals)

### Similarity Metric

Uses Python's `difflib.SequenceMatcher` ratio for fuzzy matching, with configurable threshold (default 0.8).

## API Reference

### `normalize(term, entity_type=None, threshold=0.8, top_k=1)`

The main entry point. Auto-detects entity type when not specified:
- Terms starting with `rs` + digits → variant
- Terms containing `*` → variant (star allele)
- Everything else → drug

Returns a list of result objects with fields: `raw_input`, `id`, `name`, `url`, `score`.

## Data Attribution

This package bundles data from [PharmGKB / ClinPGx](https://www.clinpgx.org/), licensed under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](http://creativecommons.org/licenses/by-sa/4.0/).

When using this package, please cite ClinPGx as described at [clinpgx.org/page/citingClinpgx](https://www.clinpgx.org/page/citingClinpgx).

The software code in this package is licensed under the MIT License. The bundled PharmGKB data retains its CC BY-SA 4.0 license.

## Development

```bash
git clone https://github.com/shloknatarajan/clinpgx-lookup.git
cd clinpgx-lookup
pip install -e ".[dev]"
pytest
```

## License

MIT (code) / CC BY-SA 4.0 (bundled PharmGKB data) — see [LICENSE](LICENSE) for details.
