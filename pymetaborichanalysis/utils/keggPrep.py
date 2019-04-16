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
import time

def parse_kgml(species_pathway_filepath):
    with open(species_pathway_filepath, "r") as infile:
        file = infile.readlines()

    name = [x for x in file if x.startswith("NAME")][0].split("   ")[2].strip()

    if any("COMPOUND" in line for line in file):
        cpd_s = [i for i, x in enumerate(file) if x.startswith("COMPOUND")][0]
        try:
            cpd_e = [i for i, x in enumerate(file[cpd_s:]) if x.strip().startswith("C") != True][0]
        except IndexError():
            cpd_e = len(cpd_s - file)

        compounds = [x.strip() for x in file[cpd_s:cpd_s+cpd_e]]

        compounds_clean = []

        for compound in compounds:
            kegg_cid = compound.replace("COMPOUND", "").strip().split(" ")[0]
            compounds_clean.append(kegg_cid)
        return name, compounds_clean
    else:
        return name, []

@click.command()
@click.option("--dir", help="File directory containing KGML files", required=True)
@click.option("--output", help="Output Directory", required=True)
def parse(dir, output):

    filenames = os.listdir(dir)
    species = list(set([x.split(":")[1][0:3] for x in filenames]))

    pathway_info = {x : {} for x in species}

    for species in species:
        species_pathways = [x for x in filenames if species in x]
        for pathway in species_pathways:
            species_pathway_filepath = os.path.join(dir, pathway)
            name, compounds = parse_kgml(species_pathway_filepath)
            pathway = pathway.split(".")[0]
            pathway_info[species][pathway] = {
                "name": name,
                "compounds": compounds
            }

    timestamp = int(time.time())

    # Taken on 16th April 2019 @ 20:26 GMT
    number_compounds = 18505

    output_dict = {
        "version": timestamp,
        "population": number_compounds,
        "species" : pathway_info
    }

    with open(os.path.join(output, "timestamp.json"), "w") as outfile:
        json.dump({"version":timestamp}, outfile)
        outfile.close()

    with open(os.path.join(output, "pathways.json"), "w") as outfile:
        json.dump(output_dict, outfile, indent=4)
        outfile.close()

if __name__ == "__main__":
    parse()
