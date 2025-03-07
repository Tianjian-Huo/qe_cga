import sys
sys.path.append('../..')
from cmds import Cluster
from structure import *
from cmds.src.atom import bond
from cmds.src.io import write_xyz

D = bond('C','C')


def cluster_from_ring(r_list):
    m = len(r_list)
    c = ring(m,'C')
    for i in range(m):
        c.add_atom(c.atoms[i]*(1+D/c.atoms[i].norm()), 'C')
    for i,r in enumerate(r_list):
        a1 = c.atoms[m+i]
        a2 = c.atoms[m+(i+1)%m]
        t = (a2-a1).unit()
        n = Point(t.y, -t.x)
        if r == 5:
            c.add_atom((a1+a2)/2+n*D*0.7, 'C')
        elif r == 6:
            c.add_atom((a1+3*a2)/4+n*D*sqrt(3)/2, 'C')
            c.add_atom((3*a1+a2)/4+n*D*sqrt(3)/2, 'C')
        elif r == 7:
            c.add_atom((0.1*a1+0.9*a2)+n*D*0.97, 'C')
            c.add_atom((0.9*a1+0.1*a2)+n*D*0.97, 'C')
            c.add_atom((a1+a2)/2+n*D*1.6, 'C')
    
    write_xyz(c, 'C_odd\\C%d_%s.xyz'%(c.get_size(), ''.join([str(_) for _ in r_list])))
 

cluster_from_ring([5]*5)
cluster_from_ring([5]*6)
cluster_from_ring([5]*7)
cluster_from_ring([6]*5)
cluster_from_ring([6]*6)
cluster_from_ring([6]*7)
cluster_from_ring([7]*5)
cluster_from_ring([7]*6)
cluster_from_ring([7]*7)

for core in [5,6,7]:
    r_list = [5]*core
    while True:
        for i in range(core):
            if r_list[i] != 7:
                break
        else:
            break
        if i == core-1:
            if r_list[-1]==5 and 7 not in r_list or r_list[-1] == 6:
                break
        r_list[i] += 1
        for j in range(i):
            r_list[j] = 5
        cluster_from_ring(r_list)
        print(r_list)
            
#c = cluster_from_ring([6,5,6,7,6,6,5])            
#write_xyz(c, 'test.xyz')