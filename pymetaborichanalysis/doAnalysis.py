from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
import json
import pandas as pd

class EnrichmentAnalysis:
    def __init__(self, clusters, type="compounds", organism="hsa", p_value_cutoff=0.05):
        self.clusters = clusters
        self.type = type
        self.organism = organism
        self.p_value_cutoff = p_value_cutoff

        self.pathway_data = self._load_data()

    def _load_data(self):
        with open("/home/keo7/Desktop/kegg_pathway.json", "r") as infile:
            data = json.load(infile)

        d = {}

        if self.organism == "all":
            # TODO: Flatten the dictionary
            pass
        else:
            for pathway_id in data[self.organism]:
                pathway_info = data[self.organism][pathway_id]
                pathway_name = pathway_info["name"]
                loi = pathway_info[self.type]
                d[pathway_id] = {"name" : pathway_name, self.type : loi}
            return d

    def _calculate_q_value(self):
        pass

    def run_analysis(self, pvalue_cutoff = 0.05, p=0.05, alternative="two-sided", adj_method="bonferroni"):
        results = []
        for pathway in self.pathway_data:
            pathway_compounds = self.pathway_data[pathway][self.type]
            pathway_name = self.pathway_data[pathway]["name"]
            in_pathway = [x for x in self.clusters if x in pathway_compounds]
            p_value = binom_test(len(in_pathway), len(pathway_compounds), p, alternative)
            results.append([pathway, pathway_name, "(%i / %i)" % (len(in_pathway), len(pathway_compounds)), p_value, "; ".join(in_pathway)])


        results =  pd.DataFrame(results, columns=["Pathway ID", "Pathway Name", "Count", "p-value", "Identifiers"])

        reject, cor_p_values, _, _ = multipletests(results["p-value"].values, method=adj_method)

        results.insert(4, "adj. p-value", cor_p_values)
        results.set_index("Pathway ID", inplace=True)

        results = results[results["adj. p-value"] < pvalue_cutoff]

        results.sort_values("adj. p-value", inplace=True)


        return results

if __name__ == "__main__":

    test_compounds = ["C00033", "C00111", "C00024", "C16255", "C00118", "C00084"]

    ea = EnrichmentAnalysis(test_compounds)
    ea = ea.run_analysis()
    print(ea)
    ea.to_csv("/home/keo7/Desktop/out.csv")