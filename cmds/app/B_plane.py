import sys
sys.path.append('../src')
from cluster import Cluster
from point import Point
import element_table
from math import sqrt

bond = element_table.bond('B','B')


def plane3(n1,n2,n3):
    c = Cluster()
    start_x = (n1-1)/2.0
    for i in range(n1):
        c.add_atom(Point(start_x+i), 'B')
    for i in range(n1-n2):
        start_x += 0.5
        y = (i+1)*sqrt(3)/2
        for j in range(n1-i-1):
            c.add_atom(Point(start_x+j, y), 'B')
    start_x = (n1-1)/2.0
    for i in range(n1-n3):
        start_x += 0.5
        y = -(i+1)*sqrt(3)/2
        for j in range(n1-i-1):
            c.add_atom(Point(start_x+j, y), 'B')
    c.zoom(bond)
    return c


def plane6(n1,n2,n3):
    '''abcabc'''
    n1,n2,n3 = sorted([n1,n2,n3],reverse=True)
    c = Cluster()
    x, y = 1./2, sqrt(3)/2
    for i in range(n3):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+i):
            c.add_atom(Point(x+j, y), 'B')
    for i in range(n2-n3):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+n3-1):
            c.add_atom(Point(x+j, y), 'B')
    for i in range(1,n3):
        x += 1./2
        y -= sqrt(3)/2
        for j in range(n1+n3-1-i):
            c.add_atom(Point(x+j, y), 'B')
    c.zoom(bond)
    c.write_xyz('B%d_plane.xyz' %c.get_size())
    return c


def plane6_(n1,n2):
    '''ababab'''
    c = Cluster()
    x, y = 1./2, sqrt(3)/2
    for i in range(n2):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+i):
            c.add_atom(Point(x+j, y), 'B')
    for i in range(1,n1):
        x += 1./2
        y -= sqrt(3)/2
        for j in range(n1+n2-1-i):
            c.add_atom(Point(x+j, y), 'B')
    c.zoom(bond)
    c.write_xyz('B%d_plane.xyz' %c.get_size())
    return c

c = plane3(7,4,4)
print c.get_size()
c.write_car('temp.car')
