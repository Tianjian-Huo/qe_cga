from cluster import Cluster
from structure import *


def poly(prototype):
    adj = prototype.adjacent()
    c = Cluster()
    d = None
    for i in range(prototype.get_size()):
        ngh = []
        for j in range(prototype.get_size()):
            if adj[i,j] == 1.:
                ngh.append(prototype.coordinate[j])
        p = (prototype.coordinate[i]+ngh[0])/2
        q = (prototype.coordinate[i]+ngh[1])/2
        r = (prototype.coordinate[i]+ngh[2])/2        
        c.add_atom(p, 'B')       
        c.add_atom(q, 'B')       
        c.add_atom(r, 'B')
        c.add_atom((p+q)/2, 'B')
        c.add_atom((r+q)/2, 'B')
        c.add_atom((p+r)/2, 'B')
        d = p.dist(q)

    unique(c)
    print c.get_size()
    c.zoom(3.2/d)
    c.write_xyz('Bx.xyz')

p=Cluster()
p.read('8e.car')
poly(p)
raise
    

def cubic4():
    c2 = Cluster()
    c = cube('')
    for i in range(c.get_size()):
        ngh = []
        for j in range(c.get_size()):
            if 0.99 < c.coordinate[i].dist(c.coordinate[j]) < 1.01:
                ngh.append(c.coordinate[j])
        p = (c.coordinate[i]+ngh[0])/2
        q = (c.coordinate[i]+ngh[1])/2
        r = (c.coordinate[i]+ngh[2])/2        
        c2.add_atom(p, 'B')       
        c2.add_atom(q, 'B')       
        c2.add_atom(r, 'B')
        c2.add_atom(p/3+q*2/3, 'B')
        c2.add_atom(p*2/3+q/3, 'B')
        c2.add_atom(p/3+r*2/3, 'B')
        c2.add_atom(p*2/3+r/3, 'B')
        c2.add_atom(q/3+r*2/3, 'B')
        c2.add_atom(q*2/3+r/3, 'B')
        c2.add_atom(p/3+q/3+r/3, 'B')

    unique(c2)
    c2.zoom(7)
    c2.write_xyz('B72.xyz')
cubic4()
raise
    
c = Cluster()
c.read('archimedes466.xyz')
c2 = Cluster()
for i in range(c.get_size()):
    for j in range(i):
        d = c.coordinate[i].dist(c.coordinate[j])
        if 1.9<d<2.1:
            c2.add_atom((c.coordinate[i]+c.coordinate[j])/2, 'B')
for i in range(36):
    for j in range(i):
        d = c2.coordinate[i].dist(c2.coordinate[j])
        if abs(d-1.732)<0.1 or abs(d-1.414)<0.1:
            c2.add_atom((c2.coordinate[i]+c2.coordinate[j])/2, 'B')
c2.zoom(2)
c2.write_xyz('B466.xyz')
