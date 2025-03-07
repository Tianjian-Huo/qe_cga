'''cluster.py
团簇类. 该类实现了团簇有关的基本操作
作者: 赛琳伟
创建时间: 2014-05-04
'''
from .point import Point
from .atom import get_element_id, bond, tolerance_min, Atom
from . import utils
from math import sin,cos,atan2,pi,hypot
import numpy as np
from random import random
from copy import deepcopy
from functools import reduce


class Cluster:
    '''团簇类'''
    iso_tol = 0.02 #同构容忍度

    def __init__(self, stru=None, elem=None):
        '''初始化
        stru: Atom列表或团簇或n*3的numpy数组（此时elem为元素列表）'''
        if stru is None:
            self.atoms = []
        elif type(stru) == list:
            self.atoms = stru
        elif isinstance(stru, Cluster):
            self.atoms = stru.atoms
        elif isinstance(stru, np.ndarray):
            self.atoms = []
            for i,e in enumerate(elem):
                self.atoms.append(Atom(stru[i,0],stru[i,1],stru[i,2],e))
        self._energy = 0.
        self.fingerprint = np.zeros(8) #指纹，长度为8的向量

    def __repr__(self):
        '''转字符串'''
        s = '%f' %self._energy
        for f in self.fingerprint:
            s += '\t%.1f'%f
        s += '\n'
        for a in self.atoms:
            s += str(a) + '\n'
        return s

    ###########################################################################
    #                      以下为属性获取
    ###########################################################################
    def get_size(self):
        '''获取原子数'''
        return len(self.atoms)

    def get_energy(self):
        '''获取能量'''
        return self._energy

    def get_fingerprint(self):
        '''获取指纹'''
        return self.fingerprint

    def get_elements_count(self):
        '''获取元素类型及每种元素的个数
        返回字典'''
        elem_count = {}
        for a in self.atoms:
            if a.elem in elem_count:
                elem_count[a.elem] += 1
            else:
                elem_count[a.elem] = 1
        return elem_count

    def get_spin(self):
        '''闭壳层0 unrestricted，开壳层1'''
        return sum([get_element_id(a.elem) for a in self.atoms]) %2

    def get_degree(self, p):
        '''计算某原子的配位数
        p为下标或Atom类型
        配位数不包含p自身'''
        if isinstance(p, Atom):
            degree = 0
            for a in self.atoms:
                if p.bond(a) and not p is a:
                    degree += 1
            return degree
        elif isinstance(p, int):
            if p < 0 or p >= self.get_size():
                raise ValueError('parameter invalid in Cluster.get_degree.')
            return len([a for a in self.atoms if a.bond(self.atoms[p])]) - 1
        else:
            raise ValueError('parameter invalid in Cluster.get_degree.')

    def get_center(self):
        return reduce(Point.__add__, self.atoms) / self.get_size()

    def isomorphic(self, other, tolerance=iso_tol):
        '''判断两个团簇是否同构
        参数other: 团簇类型,要比较的团簇
        返回: True或False'''
        if np.sum(self.fingerprint) == 0:
            self.calc_fingerprint()
        if np.sum(other.fingerprint) == 0:
            other.calc_fingerprint()
        for f1,f2 in zip(self.fingerprint, other.fingerprint):
            if abs(f1-f2)/(f1+f2+0.2) > tolerance:
                return False
        return True

    def branchs(self):
        '''连通分支数'''
        union_find = utils.UnionFind(self.get_size())
        for i in range(self.get_size()):
            for j in range(i):
                if self.atoms[i].bond(self.atoms[j]):
                    union_find.union(i, j)
        return union_find.num()

    def is_legal(self, a=None):
        '''a: Atom类型'''
        if a is None: #判断是否有原子离的太近
            for i in range(self.get_size()):
                for j in range(i):
                    if not self.atoms[i].is_legal(self.atoms[j]):
                        return False
            return True
        else: #判断一个原子a是否合法
            for atom in self.atoms:
                if not a.is_legal(atom) and a is not atom:
                    return False
            return True

    def isIsolate(self, n):
        '''判断某个点是否孤立'''
        if n < 0 or n >= self.get_size():
            return False
        for i in range(self.get_size()):
            if self.atoms[n].bond(self.atoms[i]) and i != n:
                return False
        else:
            return True

    def isOut(self, n):
        '''判断某个原子是否在表面'''
        assert n >= 0 and n < self.get_size()
        for a in self.atoms:
            if self.atoms[n].dot(a-self.atoms[n]) > 0:
                return False
        return True

    def distance(self, p):
        '''计算p跟团簇的距离（到团簇最近点的距离），p可以为点或团簇'''
        if isinstance(p, Point):
            return min([p.dist(q) for q in self.atoms])
        elif isinstance(p, Cluster):
            return min([self.distance(q) for q in p])
        else:
            raise ValueError('parameter invalid in Cluster.distance.')

    ###########################################################################
    #                    以下为属性计算和设置
    ###########################################################################
    def coord_matrix(self):
        '''返回坐标矩阵。大小N*3'''
        coord = np.empty((self.get_size(),3), dtype='float')
        for i,a in enumerate(self.atoms):
            coord[i,0] = a.x
            coord[i,1] = a.y
            coord[i,2] = a.z
        return coord

    def adjacent(self):
        '''计算邻接矩阵
        返回: 矩阵类型,邻接矩阵'''
        adj = np.zeros((self.get_size(), self.get_size()))
        for i in range(self.get_size()):
            for j in range(i):
                if self.atoms[i].bond(self.atoms[j]):
                    adj[i,j] = 1
                    adj[j,i] = 1
        return adj

    def get_bonds(self):
        '''计算每个原子的成键原子'''
        bonds = []
        for i in range(self.get_size()):
            for j in range(i):
                if self.atoms[i].bond(self.atoms[j]):
                    bonds.append((i,j))
        return bonds

    def get_connections(self):
        '''计算每个原子的成键原子'''
        conn = []
        for i in range(self.get_size()):
            conn.append([])
            for j in range(i):
                if i in conn[j]:
                    conn[-1].append(j)
            for j in range(i+1, self.get_size()):
                if self.atoms[i].bond(self.atoms[j]):
                    conn[i].append(j)
        return conn

    def dist_matrix(self):
        coord = self.coord_matrix()
        return np.linalg.norm(coord[:,np.newaxis,:]-coord[np.newaxis,:,:], axis=-1)

    def calc_fingerprint(self):
        d = self.dist_matrix()
        fp = np.linspace(1.1, 6, 8)
        self.fingerprint = np.exp(-5.65*(d[:,:,np.newaxis]-fp[np.newaxis,np.newaxis,:])**2).sum(axis=1).mean(axis=0)
        return self.fingerprint

    def get_ring(self, first, second, third, conn=None):
        if not conn:
            conn = self.get_connections()
        ring = [first,second,third]
        while True:
            p12 = self.atoms[second]-self.atoms[first]
            p23 = self.atoms[third]-self.atoms[second]
            angle = p23.angle(p12)
            index = []
            for forth in conn[third]:
                if forth == second:
                    index.append(7.)
                    continue
                for i in range(len(ring)):
                    if ring[i] == forth:
                        return ring[i:]
                p34 = self.atoms[forth]-self.atoms[third]
                angle2 = (p34).angle(p23)
                angle_diff = abs(angle2-angle)
                biangle = (p23.cross(p12)).angle(p34.cross(p23))
                index.append(angle_diff + biangle)
            forth = conn[third][np.argmin(index)]
            ring.append(forth)
            first,second,third = second,third,forth
        return ring

    def diameter(self):
        '''直径,即最远的两原子距离'''
        d = 0.
        for i in range(self.get_size()):
            for j in range(i):
                d = max(d, self.atoms[i].dist(self.atoms[j]))
        return d

    def best_transpose(self, other):
        '''找到self变为other的最好的旋转和平移操作'''
        if self.get_size() != other.get_size():
            return None
        o1 = self.get_center()
        o2 = other.get_center()
        X = np.matrix(np.zeros((3,self.get_size())))
        for i,a in enumerate(self.atoms):
            q = a - o1
            X[0,i],X[1,i],X[2,i] = q.x, q.y, q.z
        Y = np.matrix(np.zeros((3,self.get_size())))
        for i,a in enumerate(other.atoms):
            q = a - o2
            Y[0,i],Y[1,i],Y[2,i] = q.x, q.y, q.z
        w = X*Y.T
        U,S,VT = np.linalg.svd(w)
        R = (U*VT).T
        t = np.array([o2.x,o2.y,o2.z]) - np.dot(R, np.array([o1.x,o1.y,o1.z]))
        return R,np.array([t[0,0],t[0,1],t[0,2]]) #十分恼火

    def symmetry(self, tolerance=0.1):
        '''计算对称性
        参数 tolerance: 浮点类型, 对称性的容忍度
        返回: 字符串类型, 团簇的对称性'''
        from .optimize import Symmetry
        sym = Symmetry(self, tolerance=tolerance)
        sym()
        return sym.sym

    def set_energy(self, energy):
        '''设置能量'''
        self._energy = energy

    def set_element(self, elem):
        '''设置元素名
        被调用: structure.py中的函数'''
        for a in self.atoms:
            a.elem = elem

    ###########################################################################
    #                      以下为结构生成
    ###########################################################################
    def deepcopy(self, other=None):
        '''从other深度拷贝
		2种方式：dest=deepcopy(src)或dest.deepcopy(src)'''
        if other is None:
            result = self.__class__()
            result.__dict__ = deepcopy(self.__dict__)
            return result
        else:
            self.__dict__ = deepcopy(other.__dict__)

    def adj2coo(self, adj):
        '''该函数有问题，元素如何解决'''
        self.atoms = []
        current = [0]
        remain = list(range(1, self.order))
        z = 0
        while len(remain) != 0:
            new_current = []
            for p in current:
                for i in range(self.order):
                    if adj[p,i] == 1 and i in remain:
                        new_current.append(i)
                        remain.remove(i)
            z += 2
        current = new_current
        r = sin(pi/len(current))
        for i in range(len(current)):
            self.atoms.append(Point(r*cos(2*pi/len(current)), r*sin(2*pi/len(current)), z))
        self.clean()

    def from_numpy(self, a):
        '''原来有原子，只是改坐标'''
        for i in range(self.get_size()):
            self.atoms[i].x = a[3*i]
            self.atoms[i].y = a[3*i+1]
            self.atoms[i].z = a[3*i+2]

    ###########################################################################
    #                      以下为结构修改
    ###########################################################################
    def zoom(self, c):
        '''缩放
        c为缩放倍数'''
        for a in self.atoms:
            a *= c

    def move(self, p):
        '''平移
        p为Point'''
        for a in self.atoms:
            a -= p

    def center(self):
        '''将重心放到原点'''
        self.move(self.get_center())

    def rotate(self, alpha=0, beta=0, gamma=0, array=None, fm=None,to=None):
        '''旋转
        参数alpha, beta, gamma: 浮点型,为三个轴的旋转角度,弧度制.beta=gamma=0表示绕x轴旋转,其余类同
        参数array: 3*3矩阵,若前三个参数不指定,直接提供旋转矩阵array也可'''
        if array is not None: #根据旋转矩阵旋转
            for a in self.atoms:
                a.x, a.y, a.z = array[0,0]*a.x + array[0,1]*a.y + array[0,2]*a.z, \
                    array[1,0]*a.x + array[1,1]*a.y + array[1,2]*a.z, \
                    array[2,0]*a.x + array[2,1]*a.y + array[2,2]*a.z
        elif fm is not None and to is not None: #将点fm旋转到点to
            f = fm.unit()
            t = to.unit()
            n = f.cross(t).unit()
            c = f.dot(t)
            s = f.cross(t).norm()
            array = np.array([[n.x*n.x+(1-n.x*n.x)*c, n.x*n.y*(1-c)-n.z*s, n.x*n.z*(1-c)+n.y*s],
                              [n.x*n.y*(1-c)+n.z*s, n.y*n.y+(1-n.y*n.y)*c, n.y*n.z*(1-c)-n.x*s],
                              [n.x*n.z*(1-c)-n.y*s, n.y*n.z*(1-c)+n.x*s, n.z*n.z+(1-n.z*n.z)*c]])
            print(array)
            return self.rotate(array=array)
        else: #绕xyz轴分别旋转alpha,beta,gamma
            array = np.zeros((3,3))
            array[0,0] = cos(gamma)*cos(beta)
            array[0,1] = -sin(gamma)*cos(alpha) - cos(gamma)*sin(beta)*sin(alpha)
            array[0,2] = sin(gamma)*sin(alpha)- cos(gamma)*sin(beta)*cos(alpha)

            array[1,0] = sin(gamma)*cos(beta)
            array[1,1] = cos(gamma)*cos(alpha) - sin(gamma)*sin(beta)*sin(alpha)
            array[1,2] = -cos(gamma)*sin(alpha) - sin(gamma)*sin(beta)*cos(alpha)

            array[2,0] = sin(beta)
            array[2,1] = cos(beta)*sin(alpha)
            array[2,2] = cos(beta)*cos(alpha)
            return self.rotate(array=array)

    def rigid_transform(self,P,Q,R):
        '''将团簇旋转平移，使得P.x=P.y=P.z=0; Q.y=Q.z=0; R.z=0'''
        PP = deepcopy(P)
        QQ = deepcopy(Q)
        RR = deepcopy(R)
        self.atoms.extend([Atom(PP),Atom(QQ),Atom(RR)])
        self.rotate(0, 0, atan2(QQ.y-PP.y, PP.x-QQ.x))
        self.rotate(0, atan2(QQ.z-PP.z, PP.x-QQ.x), 0)
        self.rotate(atan2(RR.z-PP.z, PP.y-RR.y), 0, 0)
        self.move(deepcopy(PP))
        self.atoms = self.atoms[:-3]

    def add_atom(self, coord, elem=''):
        '''增加一个原子
        参数elem: 字符串类型,要增加的元素名称
        参数coord: 要增加的原子坐标。Atom类型 当elem为''时候为Atom类型，否则为Point类型
        返回: 无'''
        if elem == '':
            self.atoms.append(deepcopy(coord))
        else:
            self.atoms.append(Atom(coord, elem))

    def __add__(self, other):
        c = Cluster()
        c.atoms = self.atoms+other.atoms
        return c

    def __iadd__(self, other):
        '''+=
        参数other: Cluster类型,要+=的团簇
        返回: 团簇自己'''
        self.atoms.extend(other.atoms)
        return self

    def _merge(self, c1, c2, merge=True, xy_adjust=True, z_rotate=True):
        '''将c1和c2两个list拼成起来作为self.atoms
        c1或c2可以为空
        c2在上不动，c1在下朝上移动
        merge为False，则只计算offset，不拼合'''
        #c1或c2为空
        if not c1 or not c1: 
            if merge:
                self.atoms = deepcopy(c1+c2)
            return Point(0.)

        #初始offset
        o1 = Cluster(c1).get_center()
        o2 = Cluster(c2).get_center()
        if o1.z > o2.z: #确保c2在上
            c1, c2 = c2, c1
            o1, o2 = o2, o1
        if xy_adjust:
            offset = o2-o1
        else:
            offset = Point(0,0,o2.z-o1)
        for a1 in c1:
            for a2 in c2:
                d = a2.z - a1.z - tolerance_min*bond(a1.elem, a2.elem)
                if d < offset.z:
                    offset.z = d
        c = Cluster(deepcopy(c1))
        for a in c.atoms:
            a += offset
        
        #旋转靠近
        angle = best_angle = 0.
        if z_rotate and hypot(o2.x,o2.y) < 0.5: #c2中心不在z轴不可旋转
            while True:
                for _ in range(20):
                    gamma = 2*pi*random()
                    c.rotate(0, 0, gamma)
                    angle += gamma
                    for a2 in c2:   
                        if not c.is_legal(a2):
                            break
                    else:
                        best_angle = angle
                        offset.z += 0.1
                        c.move(Point(0,0,-0.1))
                        break
                else:
                    break
                if (o1+offset).z >= o2.z:
                    break

        offset.z -= 0.1
        #拼接
        if merge:
            self.atoms = deepcopy(c1)
            self.move(offset)
            self.rotate(0,0,best_angle)
            self.atoms.extend(c2)
            
        return offset
        
    def gather(self):
        '''将分散的团簇聚合'''
        def merge_radius(c1, c2):
            '''c1不动，c2沿中心连线运动'''
            o1 = c1.get_center()
            o2 = c2.get_center()
            lower, upper = 0, (o1-o2).norm()+5
            u = (o1-o2).unit()
            while upper-lower > 0.2:
                d = (upper-lower)/2
                c2.move(-d*u)
                for a in c2.atoms:
                    if not c1.is_legal(a):
                        upper = d
                        c2.move(d*u)
                        break
                else:
                    lower,upper = 0,d
            
        n = self.get_size()
        union_find = utils.UnionFind(self.get_size())
        for i in range(n):
            for j in range(i):
                if self.atoms[i].bond(self.atoms[j]):
                    union_find.union(i, j)
        kind1 = union_find.kind(0)
        remain = set(range(n)) - set(kind1)
        c1 = Cluster([self.atoms[i] for i in kind1])
        while remain:
            kind2 = union_find.kind(list(remain)[0])
            c2 = Cluster([self.atoms[i] for i in kind2])
            merge_radius(c1, c2)
            c1.atoms.extend(c2.atoms)
            remain -= set(kind2)
        self.atoms = c1.atoms

    def place(self):
        '''将结构旋转到x轴为最长方向,z轴为最短方向,并按x坐标从小到大排序
        返回:无'''
        self.center() #Cluster.center(self)
        coord = self.coord_matrix()
        U,S,VT = np.linalg.svd(coord)
        coord = np.dot(coord, VT.T)
        for i in range(self.get_size()):
            self.atoms[i] = Atom(coord[i,0],coord[i,1],coord[i,2], self.atoms[i].elem)

    def clean(self, bonds=None):
        '''实现MS界面中的clean按钮功能'''
        from .optimize import Clean
        cl = Clean(self, bonds=bonds)
        cl()