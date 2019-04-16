import os
import click

class KEGGPrep:
    def __init__(self):
        pass


@click.command()
@click.option("--dir", help="File directory containing KGML files", required=True)
@click.option("--output", help="Output Directory", required=True)
def parse(dir):
    print(dir)

if __name__ == "__main__":
    parse()