# pyMEAN: Metabolomic Enrichment ANalysis

The ```pyMEAN``` package is designed to facilitate semi-automated enrichment analysis for Metabolites. The process consists of input of metabolites, and a output of the results.

## Installation

To do.

## Usage

To do.

```
compound_list = [
    ...
]

ea = EnrichmentAnalysis(compound_list, organism="hsa")

results = ea.run_analysis(pvalue_cutoff=0.05)

print(results)
```

## License

Code released under the [GPLv3](https://github.com/KeironO/pymetabenrichanalysis/blob/master/LICENSE).

