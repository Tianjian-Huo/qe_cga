import sys
sys.path.append('../src')
from material import *
from math import *

n = 20
r = 8
s = 1.8
c = Cluster()
for i in range(n):
    O = r*Point(cos(2*pi*i/n), sin(2*pi*i/n))
    c.add_atom(O, 'B')    
    c.add_atom((O+Point(0,0,s)).rotate(2*pi*i/n, O.cross(Point(0,0,s)), O), 'B')  
    c.add_atom((O+Point(0,0,-s)).rotate(-2*pi*i/n, O.cross(Point(0,0,-s)), O), 'B')
c.zoom(0.7)
c.write_car('B60_mobius.car')