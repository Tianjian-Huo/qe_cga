import sys
sys.path.append(r'..\..')
from cmds.src.io import read
from cmds.src.optimize import *

set_optimize_method('dmol3')
name = r'opt.xyz'

c = read(name)
print(name,end=', ')
optimize(c, folder_name=os.path.splitext(name)[0], 
         task_name=os.path.splitext(os.path.basename(name))[0])
print('energy is',c.get_energy())
