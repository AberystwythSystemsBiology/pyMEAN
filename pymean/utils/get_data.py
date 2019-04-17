import os
from urllib import request

DATA_URL = "http://users.aber.ac.uk/keo7/pymean/"

def _check_dotdir(dotfile_dir):
    if not os.path.isdir(dotfile_dir):
        os.makedirs(dotfile_dir)

def _check_if_exists(dotfile_dir, database):
    print(os.path.join(dotfile_dir, database+".json"))
    if not os.path.join(dotfile_dir, database+".json"):
        return False
    else:
        return True

def _download(dotfile_dir, database):
    r = request.urlopen(DATA_URL + database+".json").read()
    print("Hello World")

def get_data(database="kegg"):
    home_dir = os.path.expanduser("~")
    dotfile_dir = os.path.join(home_dir, ".pymean/")
    _check_dotdir(dotfile_dir)
    exists = _check_if_exists(dotfile_dir, database)
    print(exists)
    if exists:
        pass
    else:
        _download(dotfile_dir, database)



if __name__ == "__main__":
    get_data()
