from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
import json
import pandas as pd
import os

class EnrichmentAnalysis:
    def __init__(self, compound_list, type="compounds", database="kegg", organism="hsa"):
        self.compound_list = compound_list
        self.type = type
        self.organism = organism

        self.pathway_data = self._load_data()

    def _load_data(self):
        # This needs to be refactored quite a lot.
        




    def run_analysis(self, pvalue_cutoff=0.05, p=0, alternative="two-sided", adj_method="bonferroni"):
        results = []
        for pathway in self.pathway_data:
            pathway_compounds = self.pathway_data[pathway][self.type]
            pathway_name = self.pathway_data[pathway]["name"]
            in_pathway = [x for x in self.clusters if x in pathway_compounds]
            p_value = binom_test(len(in_pathway), len(pathway_compounds), 1/p, alternative)
            results.append([pathway, pathway_name, "(%i / %i)" % (len(in_pathway), len(pathway_compounds)), p_value, "; ".join(in_pathway)])

        results = pd.DataFrame(results, columns=["Pathway ID", "Pathway Name", "Count", "p-value", "Identifiers"])

        reject, cor_p_values, _, _ = multipletests(results["p-value"].values, method=adj_method)

        results.insert(4, "adj. p-value", cor_p_values)
        results.set_index("Pathway ID", inplace=True)

        results = results[results["adj. p-value"] < pvalue_cutoff]

        results.sort_values("adj. p-value", inplace=True)

        return results
