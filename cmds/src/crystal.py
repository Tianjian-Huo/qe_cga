'''crystal.py
晶体类. 该类实现了晶体有关的基本操作
作者: 赛琳伟
创建时间: 2014-11-13
'''
from .cluster import Cluster
from .point import Point
from math import pi,sin,cos,acos,sqrt
import numpy as np
from copy import deepcopy

    
def param6_cell(a, b, c, alpha=90., beta=90., gamma=90.):
    '''将a,b,c,alpha,beta,gamma转成晶胞矩阵'''
    eps = 2 * np.spacing(90.0, dtype=np.float) #
    if abs(abs(alpha) - 90) < eps:
        cos_alpha = 0.0
    else:
        cos_alpha = cos(alpha * pi / 180.0)
    if abs(abs(beta) - 90) < eps:
        cos_beta = 0.0
    else:
        cos_beta = cos(beta * pi / 180.0)
    if abs(gamma - 90) < eps:
        cos_gamma = 0.0
        sin_gamma = 1.0
    elif abs(gamma + 90) < eps:
        cos_gamma = 0.0
        sin_gamma = -1.0
    else:
        cos_gamma = cos(gamma * pi / 180.0)
        sin_gamma = sin(gamma * pi / 180.0)
        
    cell = np.eye(3, dtype=np.float) #a方向为(1,0,0)
    cell[1] = b * np.array([cos_gamma, sin_gamma, 0]) #b与a在xy平面
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_sqr = 1. - cx * cx - cy * cy
    assert cz_sqr >= 0
    cz = sqrt(cz_sqr)
    cell[2] = c * np.array([cx, cy, cz])
    
    return cell
    
def cell_param6(cell):
    '''将晶胞矩阵转成a,b,c,alpha,beta,gamma'''
    a,b,c = [np.linalg.norm(v) for v in cell]
    if b*c > 1e-8:
        alpha = acos(np.dot(cell[1], cell[2]) / (b*c)) * 180/pi
    else:
        alpha = 90.
    if a*c > 1e-8:
        beta = acos(np.dot(cell[1], cell[2]) / (b*c)) * 180/pi
    else:
        beta = 90.
    if a*b > 1e-8:
        gamma = acos(np.dot(cell[1], cell[2]) / (b*c)) * 180/pi
    else:
        gamma = 90.
    return a,b,c,alpha,beta,gamma
    
    
class Crystal(Cluster):
    '''晶体类'''

    def __init__(self, atoms=None, cell=None, group=None):
        Cluster.__init__(self, atoms)
        if cell is None:
            self.cell = np.eye(3, dtpye='float')
        else:
            self.cell = cell
        if group is None:
            self.group = 'P1'
        else:
            self.group = group

    def deepcopy(self, other=None):
        '''使用copy模块的deepcopy有时候会造成错误'''
        if other is None:
            result = Cluster.deepcopy(self)
            result.cell = self.cell.copy()
            result.group = deepcopy(self.group)
            return result
        else:
            Cluster.deepcopy(self, other)
            self.cell = other.cell.copy()
            self.group = deepcopy(other.group)

    def __repr__(self):
        '''转字符串'''
        a,b,c,alpha,beta,gamma = cell_param6(self.cell)
        s = 'pbc: %f %f %f %f %f %f %s\n'%(a,b,c,alpha,beta,gamma,self.group)
        s += str(self._energy) + '\t' + '\t'.join(self.fingerprint) + '\n'
        for a in self.atoms:
            s += str(a) + '\n'
        return s

    def dimension(self):
        return self.a - max([p.x for p in self.atom]) + min([p.x for p in self.atom]) < 8 + \
            self.b - max([p.y for p in self.atom]) + min([p.y for p in self.atom]) < 8 + \
            self.c - max([p.z for p in self.atom]) + min([p.z for p in self.atom]) < 8

    def is_cluster(self):
        return self.dimension() == 0

    def center(self):
        self.move(self.get_center() + Point(self.a, self.b, self.c) / 2)

    def get_energy(self):
        return self._energy / self.get_size()

    def bond(self, i, j):
        return self.atoms[i].bond(self.atoms[j]) or \
            self.atoms[i].bond(self.atoms[j]+Point(self.a)) or self.atoms[i].bond(self.atoms[j]-Point(self.a)) or \
            self.atoms[i].bond(self.atoms[j]+Point(0,self.b)) or self.atoms[i].bond(self.atoms[j]-Point(0,self.b)) or \
            self.atoms[i].bond(self.atoms[j]+Point(0,0,self.c)) or self.atoms[i].bond(self.atoms[j]-Point(0,0,self.c))

    def zoom(self, c):
        Cluster.zoom(self, c)
        self.a *= c
        self.b *= c
        self.c *= c
