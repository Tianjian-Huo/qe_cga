from cluster import Cluster
import numpy as np


def sphere_ratio(file_name):
    c = Cluster()
    c.read(file_name)
    c.center()
    coord = np.empty((c.get_size(), 3))
    for i in xrange(c.get_size()):
        coord[i,0] = c.coordinate[i].x
        coord[i,1] = c.coordinate[i].y
        coord[i,2] = c.coordinate[i].z
    coord = np.mat(coord)
    U,S,V = np.linalg.svd(coord)
    print S[2] / S[0]
    
sphere_ratio('B56_bilayer.car')    