from bioservices.kegg import KEGG
import json
from tqdm import tqdm

class KEGGPathway:
    def __init__(self, pathway_id, _k):
        self.pathway_id = pathway_id
        self._k = _k

        self.details = self._parse()
        self.name = self._get_name()
        self.genes = self._get_genes()
        self.compounds = self._get_compounds()

    def _get_name(self):
        data = self._k.get(self.pathway_id)
        data = self._k.parse(data)
        if type(data) == dict:
            name = data["NAME"][0]
        elif type(data) == str:
            for line in data.split("\n"):
                if line.startswith("NAME"):
                    name = line.split("NAME")[1].rstrip().lstrip()
        return str(name)

    def _parse(self):
        return self._k.parse_kgml_pathway(self.pathway_id)

    def _get_genes(self):
        entries = self.details["entries"]
        genes = [x["name"].split(" ") for x in entries if x["type"] == "gene"]
        genes = [item for sublist in genes for item in sublist]
        return genes

    def _get_compounds(self):
        entries = self.details["entries"]
        compounds = [x["name"].split(":")[1] for x in entries if x["type"] == "compound"]
        return compounds

    def to_dict(self):
        return {"name" : self.name, "compounds" : self.compounds, "genes" : self.genes}

if __name__ == "__main__":
    k = KEGG()

    pathway_dict = {}

    for organism in k.organismIds:
        k.organism = organism
        for pathway_id in tqdm(k.pathwayIds):
            pathway = KEGGPathway(pathway_id, k)
            if organism not in pathway_dict:
                pathway_dict[organism] = {}
            pathway_dict[organism][pathway_id] = pathway.to_dict()

        break

    with open("/home/keo7/Desktop/kegg_pathway.json", "w") as outfile:
        json.dump(pathway_dict, outfile, indent=4)
