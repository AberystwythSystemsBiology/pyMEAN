import os
from urllib import request

DATA_URL = "http://users.aber.ac.uk/keo7/pymean/"

def _check_dotdir(dotfile_dir):
    if not os.path.isdir(dotfile_dir):
        os.makedirs(dotfile_dir)

def _check_if_exists(fp):
    return os.path.isfile(fp)


def _download(dotfile_dir, database):
    r = request.urlopen(DATA_URL + database+".json").read()
    print("Hello World")

def get_data(database="kegg", species="hsa"):
    home_dir = os.path.expanduser("~")
    dotfile_dir = os.path.join(home_dir, ".pymean/")
    _check_dotdir(dotfile_dir)
    fp = os.path.join(dotfile_dir, "%s_%s_pathways.json" % (database, species))
    exists = _check_if_exists(fp)
    if _check_if_exists(fp):
        pass
    else:
        _download(dotfile_dir, database)



if __name__ == "__main__":
    get_data()
