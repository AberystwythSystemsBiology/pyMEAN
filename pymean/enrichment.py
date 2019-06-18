from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
import pandas as pd
import os
from .utils import get_data
from itertools import chain

class EnrichmentAnalysis:
    def __init__(self, compound_list, database: str ="kegg", organism: str="hsa"):
        self.compound_list = compound_list
        self.database = database
        self.organism = organism
        self.pathway_data = self._load_data()

    def _load_data(self) -> dict:
        return get_data(self.database, self.organism)


    def _check_if_in_pathway(self, pathway_compounds: dict) -> dict:
        queried_dict = { k : [] for k in pathway_compounds.keys()}

        for inchi in self.compound_list:
            for compound, inchis in pathway_compounds.items():
                if inchi in inchis:
                    queried_dict[compound].append(inchi)

        return {k: v for k,v in queried_dict.items() if v != []}


    def _generate_in_pathway_string(self, in_pathway: dict) -> str:
        return "\t".join([" -> ".join([x, ";".join(in_pathway[x])]) for x in in_pathway])

    def _calculate_importance(self, in_pathway: dict, pathway_compounds: dict) -> float:

        return len(in_pathway) / len(pathway_compounds)


    def run_analysis(self, pvalue_cutoff: float=0.05, alternative: str="two-sided", adj_method: str="bonferroni", limiter: int= 0) -> pd.DataFrame:

        results = []

        population = self.pathway_data["population"]


        for pathway in self.pathway_data["pathways"]:

            pathway_info = self.pathway_data["pathways"][pathway]

            pathway_name = pathway_info["name"]

            pathway_compounds = list(pathway_info["compounds"].values())


            if len(pathway_compounds) >= limiter:
                in_pathway = self._check_if_in_pathway(pathway_info["compounds"])

                if len(in_pathway) != 0:

                    p_value = binom_test(len(in_pathway), len(pathway_compounds), 1/population, alternative)

                    in_pathway_str = self._generate_in_pathway_string(in_pathway)

                    importance = self._calculate_importance(in_pathway, pathway_compounds)

                    results.append([pathway, pathway_name, "(%i / %i)" % (len(in_pathway), len(pathway_compounds)), p_value, importance, in_pathway_str])

        results = pd.DataFrame(results, columns=["Pathway ID", "Pathway Name", "Count", "p-value", "Importance","Identifiers"])

        reject, cor_p_values, _, _ = multipletests(results["p-value"].values, method=adj_method)

        adj_method_str = "%s adj. p-value" % (adj_method)

        results.insert(4, adj_method_str, cor_p_values)
        results.set_index("Pathway ID", inplace=True)

        results = results[results[adj_method_str] < pvalue_cutoff]

        results.sort_values(adj_method_str, inplace=True)

        return results

