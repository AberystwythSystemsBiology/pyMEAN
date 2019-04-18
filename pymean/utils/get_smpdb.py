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

def parse_smpdb():
    pass

@click.command()
@click.option("--dir", help="File directory containing SMPDB files", required=True)
@click.option("--output", help="Output Directory", required=True)
def parse(dir, output):

    filenames = os.listdir(dir)
    species = list(set([x.split(":")[1][0:3] for x in filenames]))

    pathway_info = {x : {} for x in species}

    # TODO: Individual file for each species.
    timestamp = int(time.time())

    # SMPDB Version 2.75
    number_compounds = 55700

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
