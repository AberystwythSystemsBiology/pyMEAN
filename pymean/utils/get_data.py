import os
from urllib import request
import json

DATA_URL = "http://users.aber.ac.uk/keo7/pymean/"
PATHWAYS_FILE = "%s_%s_pathways.json"
TIMESTAMP_FILE = "%s_%s_timestamp.json"

def _load_data(dotfile_dir, pathways_file):
    with open(os.path.join(dotfile_dir, pathways_file), "r") as infile:
        return json.load(infile)

def _check_if_old(dotfile_dir, timestamp_file):
    timestamp = json.loads(request.urlopen(DATA_URL + timestamp_file).read())

    with open(os.path.join(dotfile_dir, timestamp_file), "r") as infile:
        our_timestamp = json.load(infile)

    if timestamp["version"] > our_timestamp["version"]:
        return True
    else:
        return False


def _check_dotdir(dotfile_dir):
    if not os.path.isdir(dotfile_dir):
        os.makedirs(dotfile_dir)

def _check_if_exists(fp):
    return os.path.isfile(fp)


def _download(dotfile_dir, pathways_file, timestamp_file):
    r = request.urlopen(DATA_URL + pathways_file).read()
    r = json.loads(r)
    with open(os.path.join(dotfile_dir, pathways_file), "w") as outfile:
        json.dump(r, outfile, indent=4)

    r = request.urlopen(DATA_URL + timestamp_file).read()
    r = json.loads(r)

    with open(os.path.join(dotfile_dir, timestamp_file), "w") as outfile:
        json.dump(r, outfile, indent=4)


def get_data(database="kegg", species="hsa"):
    home_dir = os.path.expanduser("~")
    dotfile_dir = os.path.join(home_dir, ".pymean/")
    _check_dotdir(dotfile_dir)
    pathways_file = PATHWAYS_FILE % (database, species)
    timestamp_file = TIMESTAMP_FILE % (database, species)
    fp = os.path.join(dotfile_dir, pathways_file)
    exists = _check_if_exists(fp)
    if _check_if_exists(fp):
        if _check_if_old(dotfile_dir, timestamp_file):
            _download(dotfile_dir, pathways_file, timestamp_file)
    else:
        _download(dotfile_dir, pathways_file, timestamp_file)
    return _load_data(dotfile_dir, pathways_file)
