from cluster import Cluster
from analysis import aver_bond_len
from utils import traverse_file

def f(file_name):
    c = Cluster()
    return aver_bond_len(c, 2.5)

traverse_file('.', 'car', f)
