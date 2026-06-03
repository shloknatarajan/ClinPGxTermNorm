# clinpgx-term-lookup

Fuzzy term lookup for [ClinPGx](https://www.clinpgx.org/) / [PharmGKB](https://www.pharmgkb.org/):
take a free-text drug or variant term and resolve it to its ClinPGx record.

It is **API-only** — no local data files to download. Lookups are served by the
public [PharmGKB](https://api.pharmgkb.org/) and
[RxNorm](https://rxnav.nlm.nih.gov/) APIs, so results stay current and the
package stays small. An internet connection is required at lookup time.

## Installation

```bash
pip install clinpgx-term-lookup
```

## Usage

### Drugs

```python
from clinpgx_term_lookup import DrugLookup

results = DrugLookup().search("warfarin")
for r in results:
    print(r.name, r.id, r.url, r.score, r.source)
```

The drug lookup:

1. Tries an exact PharmGKB chemical name match.
2. On a miss, fuzzy-matches the term with RxNorm's `approximateTerm` endpoint,
   resolves it to its ingredient (so brand names like `Tylenol` map to
   `acetaminophen`), and re-queries PharmGKB by that name.
3. On a miss, returns the RxNorm result itself as a fallback (with
   `source="rxnorm"`).

So misspellings (`warfarn`) and trade names (`tylenol`) both resolve to the
right PharmGKB chemical.

### Variants

```python
from clinpgx_term_lookup import VariantLookup

VariantLookup().search("rs1234")        # rsID  -> variant endpoint
VariantLookup().search("CYP2C19*2")     # star allele -> haplotype endpoint
```

Terms starting with `rs` are looked up as rsIDs; everything else is treated as
a star allele.

### Result objects

Both lookups return a list of Pydantic models with these fields:

| Field       | Description                                  |
| ----------- | -------------------------------------------- |
| `raw_input` | The original query string                    |
| `id`        | PharmGKB accession ID (or `RXN<rxcui>`)      |
| `name`      | The matched name                             |
| `url`       | Link to the ClinPGx (or RxNorm) record       |
| `score`     | 0–1 fuzzy similarity to the query            |
| `source`    | `"pharmgkb"` or `"rxnorm"`                    |

## Command line

```bash
clinpgx-term-lookup warfarin --type drug
clinpgx-term-lookup rs1234 --type variant
```

Output is JSON. Use `--top-k` and `--threshold` to tune results.

## License

MIT. Note that lookups query the PharmGKB and RxNorm APIs; their data is subject
to their respective terms of use.
