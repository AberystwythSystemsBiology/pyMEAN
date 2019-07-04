from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom, fisher_exact
import numpy as np
import pandas as pd
import os
from .utils import get_data
from itertools import chain


class EnrichmentAnalysis:
    def __init__(
        self,
        compound_list,
        database: str ="kegg",
        organism: str="hsa"
    ):
        self.compound_list = compound_list
        self.database = database
        self.organism = organism
        self.pathway_data = self._load_data()

    def _load_data(self) -> dict:
        return get_data(self.database, self.organism)

    def _check_if_in_pathway(self, pth_cpds: list) -> dict:
        
        found = {k: [] for k, _ in enumerate(pth_cpds)}

        for index, cpd_inchi in enumerate(pth_cpds):
            for inchi in self.compound_list:
                if inchi in cpd_inchi:
                    found[index].append(inchi)

        return {k: v for (k, v) in found.items() if len(v) != 0}

    def _generate_string(self, pathway_hits: dict, dbids: list) -> str:
        return "\t".join(
            [" -> ".join(
                [dbids[x], ";".join(pathway_hits[x])]
                ) for x in pathway_hits])

    def _calc_cov(self, in_pathway: int, pth_cpds: int) -> float:
        try:
            return in_pathway / pth_cpds
        except ZeroDivisionError:
            return 0.0
    
    def run_analysis(
        self,
        pvalue_cutoff: float=0.05,
        method: str="hyperg",
        limiter: int= 0
    ) -> pd.DataFrame:

        results = []

        population = self.pathway_data["population"]

        num_cpds = len(self.compound_list)
        
        for pathway in self.pathway_data["pathways"]:

            pth_info = self.pathway_data["pathways"][pathway]
            pth_name = pth_info["name"]

            dbids = list(pth_info["compounds"].keys())
            pth_cpds = [pth_info["compounds"][x] for x in dbids]

            num_dbids = len(dbids)

            if num_dbids >= limiter:
                pathway_hits = self._check_if_in_pathway(pth_cpds)
                num_hits = len(pathway_hits)

                if method == "hyperg":
                    p_value = 1-hypergeom.cdf(
                        num_hits-1,
                        population-num_dbids,
                        num_cpds,
                        num_dbids)
                else:
                    _, p_value = fisher_exact(
                        [
                            [
                                num_hits,
                                num_cpds-num_hits
                            ],
                            [
                                num_dbids-num_hits,
                                ((population-num_dbids)-num_cpds)+num_hits
                            ]
                        ]
                    )

                in_pathway_str = self._generate_string(pathway_hits, dbids)
                importance = self._calc_cov(num_hits, num_dbids)
                results.append([
                    pathway,
                    pth_name,
                    num_hits,
                    num_dbids,
                    p_value,
                    importance,
                    in_pathway_str
                    ])

        results = pd.DataFrame(
            results,
            columns=[
                "Pathway ID",
                "Pathway Name",
                "Hits",
                "Pathway Compounds",
                "p",
                "Coverage",
                "Identifiers"
                ])

        results.set_index("Pathway ID", inplace=True)

        _, holm_p, _, _ = multipletests(results["p"].values, method="holm")
        results.insert(4, "Holm p",  holm_p)

        results = results[results["Holm p"] <= pvalue_cutoff]

        results.sort_values("p", inplace=True)

        self.results = results
