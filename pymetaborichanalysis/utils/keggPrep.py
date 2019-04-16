import os
import click
import json

def parse_kgml(species_pathway_filepath):
    with open(species_pathway_filepath, "r") as infile:
        file = infile.readlines()

    if any("COMPOUND" in line for line in file):
        cpd_s = [i for i, x in enumerate(file) if x.startswith("COMPOUND")][0]

        try:
            cpd_e = [i for i, x in enumerate(file[cpd_s:]) if x.strip().startswith("C") != True][0]
        except IndexError():
            cpd_e = len(cpd_s - file)



        compounds = [x.strip() for x in file[cpd_s:cpd_e]]

        if len(compounds) == 0:
            print(compounds[cpd_s:])
            print(cpd_s)
            print(cpd_e)
            print(species_pathway_filepath)

        for compound in compounds:
            kegg_cid = compound.replace("COMPOUND", "").strip().split(" ")[0]
            #print(kegg_cid)
    else:
        pass

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
            parse_kgml(species_pathway_filepath)

        break


if __name__ == "__main__":
    parse()
