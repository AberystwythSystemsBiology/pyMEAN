from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
import pandas as pd
import os
from utils.get_data import get_data
from itertools import chain

class EnrichmentAnalysis:
    def __init__(self, compound_list, database: str ="kegg", organism: str="hsa"):
        self.compound_list = compound_list
        self.database = database
        self.organism = organism
        self.pathway_data = self._load_data()

    def _load_data(self) -> dict:
        return get_data(self.database, self.organism)

    def run_analysis(self, pvalue_cutoff: float=0.05, alternative: str="two-sided", adj_method: str="bonferroni") -> pd.DataFrame:

        results = []

        population = self.pathway_data["population"]


        for pathway in self.pathway_data["pathways"]:

            pathway_info = self.pathway_data["pathways"][pathway]

            pathway_name = pathway_info["name"]

            pathway_compounds = list(pathway_info["compounds"].values())

            in_pathway = [x for x in self.compound_list if x in chain(*pathway_compounds)]

            p_value = binom_test(len(in_pathway), len(pathway_compounds), 1/population, alternative)

            results.append([pathway, pathway_name, "(%i / %i)" % (len(in_pathway), len(pathway_compounds)), p_value, "; ".join(in_pathway)])

        results = pd.DataFrame(results, columns=["Pathway ID", "Pathway Name", "Count", "p-value", "Identifiers"])

        reject, cor_p_values, _, _ = multipletests(results["p-value"].values, method=adj_method)

        results.insert(4, "adj. p-value", cor_p_values)
        results.set_index("Pathway ID", inplace=True)

        results = results[results["adj. p-value"] < pvalue_cutoff]

        results.sort_values("adj. p-value", inplace=True)

        return results

if __name__ == "__main__":

    compound_list = [
        "GPRLSGONYQIRFK-FTGQXOHASA-N",
        "IVOMOUWHDPKRLL-KQYNXXCUSA-N",
        "ITGRMIJUDWFPJB-YYHNHJMXSA-N",
        "MMWCIQZXVOZEGG-XJTPDSDZSA-N"
    ]

    ea = EnrichmentAnalysis(compound_list, organism="hsa")

    print(ea.run_analysis(pvalue_cutoff=0.05))
