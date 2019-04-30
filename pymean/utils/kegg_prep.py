'''
    # pyMEAN:: Metabolomic Enrichment ANalysis

    I don't expect you to run this file at all, it's purely used to create
    the JSON files that are required to run this litte tool.

    If you do happen to want to use this, feel free to drop me an email
    at any time and I'll get back to you with instructions on how to use
    it.
'''

import os
import click
import json
import requests
import time
import xmltodict
import bioservices
from bioservices import KEGG, ChEBI
from zeep import Client


k = KEGG(verbose=False)
map_kegg_chebi = k.conv("chebi", "compound")
c = ChEBI(verbose = False)

chebi_client = Client("https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl")



def chebi(compound_id):
    def _get_chebi(chebi_id):
        result = chebi_client.service.getCompleteEntity(chebi_id.upper())
        return result

    def _recurisve_get_parent(chebi_id):
        result = _get_chebi(chebi_id)
        if result["inchiKey"] == None:
            for child in result["OntologyParents"]:
                _recurisve_get_parent(child["chebiId"])
        else:
            return result["inchiKey"]

    chebi_id = map_kegg_chebi[compound_id]
    result = _get_chebi(chebi_id)

    if result["inchiKey"] == None:
        inchikeys = []
        for child in result["OntologyChildren"]:
            child_chebi_id = child["chebiId"]
            inchikeys.append(_recurisve_get_parent(child_chebi_id))
        return inchikeys
    else:
        return [result["inchiKey"]]



def bridgedb(compound_id):
    compound_id = compound_id.split(":")[1]
    url = "https://webservice.bridgedb.org/Human/xrefs/Ck/%s" % (compound_id)
    identifiers = requests.get(url)
    identifiers = identifiers.text.split("\n")
    inchikeys = []
    for identifier in identifiers:
        if identifier != "":
            id, database = identifier.split("\t")
            if database == "InChIKey":
                inchikeys.append(id)
    return inchikeys

def parse_kgml(species_pathway_filepath):
    with open(species_pathway_filepath, "r") as infile:
        file = infile.read()


    kegg_pathway = xmltodict.parse("".join(file))["pathway"]

    pathway_name = kegg_pathway["@title"]

    compounds = [e["@name"] for e in kegg_pathway["entry"] if e["@type"] == "compound"]

    # Sometimes, CPDs are combined.
    compounds = [x.split(" ") for x in compounds]
    compounds = [item for sublist in compounds for item in sublist]

    for index, compound in enumerate(compounds):
        inchikeys = bridgedb(compound)
        if len(inchikeys) == 0:
            inchikeys = chebi(compound)
            if len(inchikeys) == 0:
                print("Cannae find %s" % compound)

    return "NotASwearWord", "Off"


@click.command()
@click.option("--dir", help="File directory containing KGML files", required=True)
@click.option("--output", help="Output Directory", required=True)
def parse(dir, output):

    filenames = os.listdir(dir)
    species = list(set([x.split(":")[1][0:3] for x in filenames]))

    pathway_info = {x : {} for x in species}

    # TODO: Individual file for each species.
    timestamp = int(time.time())

    # Taken on 16th April 2019 @ 20:26 GMT
    number_compounds = 18505

    for species in species:
        species_pathways = [x for x in filenames if species in x]
        species_pathways_dict = {
            "pathways" : {},
            "version" : timestamp,
            "population" : number_compounds
        }

        for pathway in species_pathways:
            species_pathway_filepath = os.path.join(dir, pathway)
            name, compounds = parse_kgml(species_pathway_filepath)
            pathway = pathway.split(".")[0]
            species_pathways_dict["pathways"][pathway] = {
                "name": name,
                "compounds": compounds,
            }

        with open(os.path.join(output, "kegg_%s_timestamp.json" % species), "w") as outfile:
            json.dump({"version":timestamp}, outfile)
            outfile.close()

        with open(os.path.join(output, "kegg_%s_pathways.json" % species), "w") as outfile:
            json.dump(species_pathways_dict, outfile, indent=4)
            outfile.close()

if __name__ == "__main__":
    parse()
