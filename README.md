# pyMEAN: Metabolomic Enrichment ANalysis

The pyMEAN package is designed to facilitate semi-automated enrichment analysis for metabolomic experiments.

## Installation

pyMEAN requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

The easiest way of installing pyMEAN is using pip:

```
pip install pymean
```

Alternatively, you can use git and pip in unison to get the development branch:


```
pip install https://github.com/KeironO/pyMEAN
```

## Usage

```
# import package into python
>>> from pymean import EnrichmentAnalysis
# create a compound list of inchikeys
>>> compound_list = [...]
# create EnrichmentAnalysis object for the analysis of hsa (Homo Sapiens)
>>> ea = EnrichmentAnalysis(compound_list, organism="hsa")
# run the analysis
>>> results = ea.run_analysis(pvalue_cutoff=0.05)
```

## License

Code released under the [GPLv3](https://github.com/KeironO/pymetabenrichanalysis/blob/master/LICENSE).

