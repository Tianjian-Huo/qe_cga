import sys
sys.path.append('../..')
from cmds.src.abinitio import Abinitio
from cmds.search.genetic_algorithm import CGA
import sys
import os

if len(sys.argv) > 1:
    os.chdir(sys.argv[1])
if len(sys.argv) > 2:
    Abinitio.set_path(sys.argv[2])

ga = CGA()
ga.init_pop()
ga.iterator()
