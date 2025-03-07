'''fix_cluster.py
用于搜索的团簇类，原子数和元素类型固定
作者: 赛琳伟
创建时间: 2014-05-04
'''
from .point import Point
from . import utils
from .atom import bond as BOND
from .atom import get_element_id,tolerance_min,tolerance_max,Atom
from .cluster import Cluster
from .io import read
from .optimize import optimize, get_optimize_method
from math import pi,sqrt,sin,cos,acos,hypot
from copy import deepcopy
from random import random,randrange,choice,uniform,sample,shuffle,randint
from configparser import ConfigParser
import os
import numpy as np
from functools import reduce

__all__ = ['FixCluster','MonoCluster', 'SymCsCluster', 'SymC2Cluster', 'SymC3Cluster',
           'SymC5Cluster', 'SymCiCluster', 'SymTCluster', 'Bilayer', 'Planar',
           'BiCluster', 'SymCsBiCluster', 'SymC2BiCluster', 'SymC3BiCluster',
           'Alternate', 'MultiCluster', 'SymMultiCluster', 'MoleCluster',
           'Hydrate', 'LigandCluster', 'config']


###############################################################################
#                        抽象团簇
###############################################################################
class FixCluster(Cluster):
    '''抽象类，提取公共部分，供后面所有类继承'''
    _cls_size = 0 #原子个数
    amplitude = 0.3 #微扰幅度

    @classmethod
    def config(cls):
        cfg = ConfigParser()
        if os.path.isfile('config.ini'):
            cfg.read('config.ini')
        else:
            cfg.read(os.path.join(os.path.split(__file__)[0],'config.ini'))
        section = cfg['cluster']
        if 'iso_tol' in section:
            Cluster.iso_tol = section.getfloat('iso_tol')
        if 'perturbation_amplitude' in section:
            FixCluster.amplitude = section.getfloat('perturbation_amplitude')
        return cfg

    def check(self):
        '''检查结构是否合法'''
        pass
        
    def read(self, file_name):
        self.atoms = read(file_name).atoms
        self.after_read()
        
    def after_read(self):
        self.center()
        self.calc_fingerprint()

    def random(self):
        pass

    def mating(self, father, mother):
        pass

    def perturbation(self, step=amplitude):
        '''给结构微小扰动
        1个原子扰动10次的结果是偏离3*step,100次就是9*step,1000次是30*step'''
        for _ in range(self._cls_size*100):
            atom = choice(self.atoms)
            changes = Point.random(step)
            atom += changes
            degree = 0
            for a in self.atoms:
                if a is not atom:
                    if not atom.is_legal(a):
                        atom -= changes
                        break
                    if atom.bond(a):
                        degree += 1
            else:
                if degree == 0:
                    atom -= changes


    def unsuspend(self):
        #原子配位数为1的，修改为2
        for M in self.atoms:
            if self.get_degree(M) > 1:
                continue
            P, Q = None, None
            for a in self.atoms:
                if (P is None or a.dist(M) < P.dist(M)) and a is not M:
                    P = a
                elif (Q is None or a.dist(M) < Q.dist(M)) and a is not M:
                    Q = a
            for _ in range(100):
                M.set_coordinate(Q + Point.random(BOND(M.elem,Q.elem)))
                if self.is_legal(M) and self.get_degree(M) > 1:
                    break


    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for _ in range(100):
            idx = randrange(self.get_size())
            old_pos = deepcopy(self.atoms[idx])
            dest = choice(self.atoms)
            self.atoms[idx] = dest + Point.random(BOND(old_pos.elem,old_pos.elem))
            self.atoms[idx].elem = old_pos.elem
            if not old_pos.bond(self.atoms[idx]) and self.get_degree(self.atoms[idx]) > 1 \
                and self.is_legal(self.atoms[idx]) and self.branchs() == 1:
                return idx
            else:
                self.atoms[idx] = old_pos
        else:
            return -1

    def mutation_one_better(self):
        '''将一个原子移动到另外一个原子旁边
           尝试多次，选最好的'''
        best = deepcopy(self)
        best.set_energy(0.)
        best_move = None
        for _ in range(self._cls_size//10+3):
            c = deepcopy(self)
            move = c.mutation_one()
            if move == -1:
                continue
            optimize(c, fix_coord=[7]*move+[0]+[7]*(self.get_size()-1-move))
            if c.get_energy() < best.get_energy() and \
                not c.atoms[move].bond(self.atoms[move]):
                best = deepcopy(c)
                best_move = move
        if best.get_energy() < 0.:
            self.deepcopy(best)
            #如果移动后的原子配位数为1，修改为2
            if self.get_degree(best_move) == 1:
                M = self.atoms[best_move]
                P, Q = None, None
                for a in self.atoms:
                    if P is None or a.dist(M) < P.dist(M):
                        P = M
                    elif Q is None or a.dist(M) < Q.dist(M):
                        Q = M
                n = (M-P).cross(M-Q)
                if n.norm() > 0.1:
                    M.rotate(-0.1, n, P)
                else:
                    M = (P+Q)/2
        else:
            utils.write_log('mutation_one_better failure.')
            self.perturbation()

###############################################################################
#                        抽象对称团簇
###############################################################################
class FixSymCluster(FixCluster):
    sym = '' #对称性
    _order = 1 #对称群阶数

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(FixSymCluster,cls).config()
        if not get_optimize_method().check_symmetry_on():
            print('please set symmetry on in input file.')
            os._exit(1)
        return cfg


    def _set_symmetry_axis(self):
        #寻找对称轴，将对称轴改为z轴（Cs为对称面，xy面）
        if self.sym == 'Cs':
            x1=x2=y1=y2=z1=z2=0.
            for a in self.atoms:
                if a.x > 0:
                    x1 += a.x * a.x
                else:
                    x2 += a.x * a.x
                if a.y > 0:
                    y1 += a.y * a.y
                else:
                    y2 += a.y * a.y
                if a.z > 0:
                    z1 += a.z * a.z
                else:
                    z2 += a.z * a.z
            if abs(z1-z2)/self._cls_size < 1e-2:
                return
            elif abs(y1-y2)/self._cls_size < 1e-2:
                self.rotate(pi/2,0,0)
                return
            elif abs(x1-x2)/self._cls_size < 1e-2:
                self.rotate(0,pi/2,0)
                return
            raise utils.SymmetryErr('not a Cs symmetry structure.')
        else:
            self.calc_fingerprint()
            temp = deepcopy(self)
            temp.rotate(0,0,pi*2/3)
            temp.calc_fingerprint()
            if self.isomorphic(temp):
                return
            temp2 = deepcopy(self)
            temp2.rotate(0,pi*2/3,0)
            temp2.calc_fingerprint()
            if self.isomorphic(temp2):
                self.rotate(pi/2,0,0)
                return
            temp3 = deepcopy(self)
            temp3.rotate(pi*2/3,0,0)
            temp3.calc_fingerprint()
            if self.isomorphic(temp3):
                self.rotate(0,pi/2,0)
                return
            raise utils.SymmetryErr('not a %s symmetry structure.' %self.sym)

    def _arrange(self):
        pass

    def _is_axis(self, a):
        '''判断p是否在对称轴上。不同的类需要具体实现'''
        pass

    def _put_axis(self, a):
        '''将p修改为在对称轴上,Cs对称例外'''
        a.x = a.y = 0

    def symmetry_atom(self, a):
        '''得到p对称的全部原子。不同的类需要具体实现'''
        pass

    def after_read(self):
        '''读取xyz文件需在原Cluster的基础上再给_nonaxis_num赋值'''
        FixCluster.after_read(self)
        self._arrange()

    def perturbation(self, step=FixCluster.amplitude):
        '''给结构微小扰动'''
        for _ in range(self._cls_size*100):
            old_pos = choice(self.atoms)
            new_pos = old_pos + Point.random(step)
            self._del(old_pos)
            #不得改变是否在对称轴上
            if self._is_axis(new_pos) ^ self._is_axis(old_pos):
                self._add(old_pos)
                continue
            #旧点在在对称轴上,新点也要在对称轴上
            if self._is_axis(old_pos):
                self._put_axis(new_pos)
            #两原子不可太近
            if self.is_legal(new_pos):
                self._add(new_pos)
            else:
                self._add(old_pos)

    def mutation_one_better(self):
        '''将一个原子移动到另外一个原子旁边
           尝试多次，选最好的'''
        best = deepcopy(self)
        best.set_energy(0.)
        for _ in range(self._cls_size//5+2):
            c = deepcopy(self)
            relax_atoms = c.mutation_one()
            if relax_atoms == []:
                continue
            fix_coord = [7]*self.get_size()
            for a in relax_atoms:
                fix_coord[a] = 0
            cc = Cluster(c.atoms) #类别降为Cluster
            optimize(cc, fix_coord=fix_coord)
            if cc.get_energy() >= best.get_energy():
                continue
            if len(relax_atoms) == 1: #对称轴上的原子
                if not cc.atoms[relax_atoms[0]].bond(self.atoms[relax_atoms[0]]):
                    best.atoms = deepcopy(cc.atoms)
                    best._put_axis(cc.atoms[-1])
                    best.set_energy(cc.get_energy())
                    break
            elif isinstance(self, SymCluster): #非对称轴上的原子，且为单质团簇
                for a in c.atoms[:c._order]:
                    if a.bond(self.atoms[relax_atoms[0]]):
                        break
                else:
                    best.atoms = deepcopy(cc.atoms)
                    best.set_energy(cc.get_energy())
                    best._del(best.atoms[0])
                    best._add(cc.atoms[0])
                    break
            else: #非对称轴上的原子，且为二元团簇
                if relax_atoms[0] < self._atom_num1: #元素1
                    for a in c.atoms[self._atom_num1-self._order:self._atom_num1]:
                        if a.bond(self.atoms[relax_atoms[0]]):
                            break
                    else:
                        best.atoms = deepcopy(cc.atoms)
                        best.set_energy(cc.get_energy())
                        best._del(best.atoms[self._atom_num1-1])
                        best._add(cc.atoms[self._atom_num1-1])
                        break
                else: #元素2
                    for a in c.atoms[self._atom_num1:self._atom_num1+self._order]:
                        if a.bond(self.atoms[relax_atoms[0]]):
                            break
                    else:
                        best.atoms = deepcopy(cc.atoms)
                        best.set_energy(cc.get_energy())
                        best._del(best.atoms[self._atom_num1])
                        best._add(cc.atoms[self._atom_num1])
                        break
        self.deepcopy(best)

###############################################################################
#                        单质团簇
###############################################################################
class MonoCluster(FixCluster):
    '''单质团簇类.
    在该类中原子数和元素类型是固定的.使用前须先设置这两个值.'''

    @classmethod
    def set_element(cls, element):
        if get_element_id(element) == -1:
            print('input element error.')
            os._exit(1)
        cls._element = element #元素名
        cls._bond_len = BOND(cls._element, cls._element) #键长

    @classmethod
    def get_element(cls):
        return cls._element

    @classmethod
    def set_atom_num(cls, atom_num):
        cls._cls_size = atom_num

    @classmethod
    def get_atom_num(cls):
        return cls._cls_size

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(MonoCluster, cls).config()
        cls.set_element(cfg.get('cluster', 'element'))
        cls.set_atom_num(cfg.getint('cluster', 'atom_num'))

    def check(self):
        '''检查结构是否正确
        '''
        return self.get_size() == self._cls_size and \
            self.get_elements_count()[self._element] == self._cls_size

    def get_degree(self, pt):
        '''计算团簇外某个点的度
        pt: Point或Atom类型'''
        if isinstance(pt, Point):
            degree = 0
            for a in self.atoms:
                if pt.dist(a) < tolerance_max * self._bond_len:
                    degree += 1
            return degree
        else:
            return Cluster.get_degree(self, pt)

    def random(self):
        '''产生随机结构'''
        a,b,c = abs(random())+0.1, abs(random())+0.1, abs(random()+0.1)
        d = (self._cls_size/a/b/c)**(1./3) * 0.8 * self._bond_len
        while self.get_size() < self._cls_size:
            for _ in range(20):
                atom = Atom(Point(random()*a, random()*b, random()*c)*d, self._element)
                if self.is_legal(atom):
                    self.add_atom(atom)
                    break
            else:
                d += 0.1
                self.atoms = []

    def random2(self):
        '''产生随机结构,每个点跟前一个点相连且还跟别的一个点连'''
        self.atoms = []
        self.add_atom(Point(0), self._element)
        self.add_atom(Point(1), self._element)
        for i in range(2, self._cls_size):
            for _ in range(1000): #尝试1000次
                d = uniform(tolerance_min, tolerance_max) * self._bond_len
                a = self.atoms[i-1] + Point.random(d)
                if self.get_degree(a) > 1 and self.is_legal(a):
                    self.add_atom(a)
                    break
            else:
                print('generate initial structure failure.')
                raise
                    
    def perturbation(self, step=FixCluster.amplitude):
        '''基类的变异效率太低'''
        c = self.coord_matrix()
        d = np.linalg.norm(c[np.newaxis,:,:] - c[:,np.newaxis,:], axis=-1)
        d[np.diag_indices_from(d)] = 10
        if np.any(d < self._bond_len*tolerance_min): #有原子离的太近
            self.clean()
            return
        if np.any(np.sum(d<self._bond_len*tolerance_max, axis=0) == 0): #有原子飞掉
            self.gather()
            
        src = c.copy()
        h = 0.02
        n = (int(step/h))**2
        while n:
            change = np.random.uniform(-h, h, c.shape)
            c += change
            d = np.linalg.norm(c[np.newaxis,:,:] - c[:,np.newaxis,:], axis=-1)
            d[np.diag_indices_from(d)] = 10
            #检查是否有原子离的太近或飞掉
            union_find = utils.UnionFind(self.get_size())
            for i in range(self.get_size()):
                for j in range(i):
                    if d[i,j] < self._bond_len*tolerance_max:
                        union_find.union(i, j)
            if np.all(d > self._bond_len*tolerance_min) and union_find.num() == 1:
                n -= 1
            else:
                c -= change
            if n%10==0 and np.linalg.norm(c-src, axis=-1).mean() >= step:
                break
        for i in range(self.get_size()):
            self.atoms[i].x = c[i][0]
            self.atoms[i].y = c[i][1]
            self.atoms[i].z = c[i][2]

    def mating(self, father, mother):
        '''交叉.
        沿垂直于z轴的平面交叉,取father的下半球,mother的上半球'''
        f = deepcopy(father)
        f.center()
        f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
        f.atoms.sort(key=lambda a: a.z)
        m = deepcopy(mother)
        m.center()
        m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
        m.atoms.sort(key=lambda a: a.z)

        locus = randrange(self._cls_size//4+1, self._cls_size//2+2)
        self._merge(f.atoms[:locus], m.atoms[locus:])

    def mating_principal_direction(self, father, mother):
        '''先将双亲结构旋转到x轴为最长方向,z轴为最短方向;
        然后用任意垂直于x轴的平面切割父代和母代;
        然后随机取父代的A半边,母代的A半边,合并'''
        f = father.deepcopy()
        father.place()
        f.atoms.sort(key=lambda a: a.z)
        if random() > 0.5: #选用父代的哪半边
            for a in f.atoms:
                a.z = -a.z
            f.atoms.reverse()
        father.rotate(0, pi/2, 0)

        m = mother.deepcopy()
        m.place()
        if random() > 0.5:#选用母代的哪半边
            for a in m.atoms:
                a.z = -a.z
            m.atoms.reverse()
        mother.rotate(0, pi/2, 0)

        locus = randrange(FixCluster._cls_size/10+2, FixCluster._cls_size/2+3)
        self._merge(f.atoms[:locus], m.atoms[locus:])
            
    def mutation_one_best(self):        
        #计算删哪个原子最好
        e = []
        for i in range(self.get_size()):
            c = self.deepcopy()
            del c.atoms[i]
            optimize(c, n_iter=0) #input_name='energy'
            e.append(c.get_energy())
        #删除最差位置原子
        best_del = np.argmin(e)
        del_atom = self.atoms[best_del]
        del self.atoms[best_del]
        #加在最好的位置
        min_energy = 0.
        best_stru = None
        for neighbor in self.atoms:
            if neighbor.bond(del_atom):
                continue
            for _ in range(10):
                a = neighbor + Point.random(self._bond_len)
                if self.is_legal(a) and self.get_degree(a) > 1:
                    c = self.deepcopy()
                    c.add_atom(a)
                    optimize(c, n_iter=100)
                    if c.get_energy() < min_energy:
                        min_energy = c.get_energy()
                        best_stru = c
                    break
        self.deepcopy(best_stru)


class SymCluster(MonoCluster, FixSymCluster):
    '''对称团簇的抽象类,不能直接使用
    具体类中需实现_is_axis, _add
    坐标存储方式为:选存不在对称轴上的点，对称的点相邻存放,再存储在对称轴上的点'''

    def __init__(self):
        super(SymCluster, self).__init__()
        self._nonaxis_num = 0 #不在对称轴上的点的个数
        self._axis_num = 0 #在对称轴上的点的个数

    def _arrange(self):
        '''每次从别的地方(例如读recover)读取结构,都要调用该函数调整原则顺序，并给_nonaxis_num和axis_num赋值
        dmol优化成高对称性后会改变原子排列顺序！'''
        self._set_symmetry_axis()

        self._axis_num = self._nonaxis_num = 0
        i = 0
        while i < self._cls_size - self._axis_num:
            if self._is_axis(self.atoms[i]):
                self._axis_num += 1
                self.add_atom(self.atoms[i])
                del self.atoms[i]
                continue
            self._nonaxis_num += 1
            sym_atoms = self.symmetry_atom(self.atoms[i])
            for a in sym_atoms[1:]:
                for j in range(i+1, self._cls_size-self._axis_num):
                    if (self.atoms[j]-a).norm() < 0.1:
                        self.atoms[i+1],self.atoms[j] = self.atoms[j],self.atoms[i+1] #这里可能有错，参看SymBiCluster
                        break
                else:
                    print('cluster loose symmetry.')
                    raise utils.SymmetryErr('cluster loose symmetry.')
            i += self._order

    def get_axis_num(self):
        '''获取坐标轴上的点的个数'''
        return self._axis_num

    def get_nonaxis_num(self):
        '''获取不在坐标轴上的点的个数'''
        return self._nonaxis_num

    def _add(self, a):
        '''增加一个原子
        非对称轴原子加在前端，对称轴原子加在后端
        被调用:random,perturbation'''
        if self._is_axis(a):
            self.add_atom(a)
            self._axis_num += 1
        else:
            self.atoms = self.symmetry_atom(a) + self.atoms
            self._nonaxis_num += 1

    def _add_axis_point(self):
        '''增加一个对称轴上的点
        返回True或False'''
        for _ in range(100):
            new_atom = choice(self.atoms) + Point.random_plane(self._bond_len)
            self._put_axis(new_atom)
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1:
                self._add(new_atom)
                return True
        else:
            return False

    def _add_nonaxis_point(self):
        '''增加一个非对称轴上的点
        返回True或False'''
        for _ in range(100):
            new_atom = choice(self.atoms) + Point.random(self._bond_len)
            if self._is_axis(new_atom):
                continue
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1:
                self._add(new_atom)
                return True
        else:
            return False

    def _del(self, a):
        '''删除一个原子
        参数p：Point类型或数组下标
        被调用:perturbation,_increase_naxis'''
        if isinstance(a, Atom): #a为原子
            index = self.atoms.index(a)
        else: #a为原子的下标
            index = a
            a = self.atoms[index]
        if self._is_axis(a):
            del self.atoms[index]
            self._axis_num -= 1
        else:
            index = index // self._order * self._order
            del self.atoms[index : index+self._order]
            self._nonaxis_num -= 1

    def random(self):
        '''产生随机结构'''
        self.atoms = []
        self._axis_num = self._nonaxis_num = 0
        expect = 1.75 * self._bond_len * self._cls_size**(1.0/3)
        height = sum([random() for _ in range(12)]) / 6 * expect * 0.5
        radius = sqrt(4.0/3 * self._cls_size * self._bond_len**3 / height) * 0.5
        while self.get_size() < self._cls_size:
            for _ in range(100): #每个原子尝试100次
                now = Atom(Point.random_plane(random()*radius) + Point(0,0,random()*height), self._element)
                if self._is_axis(now): #离对称轴太近,就放到对称轴上
                    self._put_axis(now)
                #未产生原子个数已不足以产生一个非轴原子,则让该原子在对称轴上
                if self._cls_size - self.get_size() < self._order:
                    self._put_axis(now)
                #判断新产生的原子是否合法
                if not self.is_legal(now):
                    continue
                #将新产生的点放入坐标
                self._add(now)
                break
            else: #尝试100次没有找到第i个原子的合适位置
                radius *= 1.05
                height *= 1.05
                self.atoms = []
                self._axis_num = self._nonaxis_num = 0

    def mating(self, father, mother):
        '''交叉'''
        f = deepcopy(father)
        f.center()
        f.rotate(0, 0, 2*pi*random())
        m = deepcopy(mother)
        m.center()
        m.rotate(0, 0, 2*pi*random())
        if isinstance(self, SymCsCluster):
            f.rotate(pi/2, 0, 0)
            m.rotate(pi/2, 0, 0)
        else:
            if random() > 0.5: #选用母代的哪半边
                for i in range(self._cls_size):
                    m.atoms[i].z = -m.atoms[i].z

        for _ in range(20):
            cut_plane_f = uniform(0, max([a.z for a in f.atoms])*0.7) #杂交面
            up = [a for a in f.atoms[:f._nonaxis_num*f._order] if a.z>=cut_plane_f]
            up_axis = [a for a in f.atoms[f._nonaxis_num*f._order:] if a.z>=cut_plane_f]

            lower = min([a.z for a in m.atoms])
            upper = max([a.z for a in m.atoms])
            while lower < upper - 0.1: #二分查找上边界
                cut_plane_m = (lower+upper)/2
                down = [a for a in m.atoms[:m._nonaxis_num*m._order] if a.z<cut_plane_m]
                down_axis = [a for a in m.atoms[m._nonaxis_num*m._order:] if a.z<cut_plane_m]
                defect_num = len(up) + len(down) + len(up_axis) + len(down_axis) - self._cls_size
                if defect_num > 1 or len(up)+len(down) > self._cls_size:
                    upper = cut_plane_m
                elif defect_num < -1:
                    lower = cut_plane_m
                else:
                    break

            if len(up)+len(down) > self._cls_size or abs(defect_num) > 1:
                continue

            self._nonaxis_num = (len(up) + len(down)) // self._order
            self._axis_num = len(down_axis) + len(up_axis)

            z_rotate = not isinstance(self, SymCsCluster)
            move = self._merge(down+down_axis, up+up_axis, False, z_rotate=z_rotate)
            for a in up:
                a += move
            for a in up_axis:
                a += move
            nup = len(up) // self._order
            ndown = len(down) // self._order
            assert nup*self._order == len(up) and ndown*self._order == len(down)
            self.atoms = deepcopy(up + down + up_axis + down_axis)

            if isinstance(self, SymCsCluster):
                self.rotate(-pi/2, 0, 0)
                self.center()

            if defect_num == -1: #原子数缺1
                if not self._add_axis_point():
                    continue
            elif defect_num == 1: #原子数多1
                del self.atoms[-randrange(len(up_axis + down_axis))-1]
                self._axis_num -= 1

            break
        else:
            print('mating failure.')
            if isinstance(self, SymCsCluster):
                f.rotate(-pi/2, 0, 0)
            self.deepcopy(f)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    def mating_radial(self, father, mother):
        '''交叉。沿径向方向，父外母内，母旋转
        不能用于Cs'''
        f = deepcopy(father)
        f.center()
        rf = max([a.x*a.x+a.y*a.y for a in f.atoms]) #半径
        f.atoms.sort(key=lambda a:a.x*a.x+a.y*a.y, reverse=True)
        m = deepcopy(mother)
        m.center()

        for _ in range(100):
            m.rotate(0, 0, 2*pi*random())
            cut_radius_f = uniform(rf*0.3, rf*0.7) #杂交面
            inner = [a for a in m.atoms if a.x*a.x+a.y*a.y<cut_radius_f]
            if len(inner) == self._cls_size or len(inner)==0:
                continue
            outer = f.atoms[:self._cls_size-len(inner)]
            try: #检查键长是否合法
                for p in outer:
                    for q in inner:
                        if not p.is_legal(q):
                            raise
            except:
                continue
            else:
                self.atoms = outer + inner
                break
        else:
            print('mating failure.')
            self.deepcopy(f)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号
        调用：_add_nonaxis_point，_add_axis_point，perturbation'''
        for _ in range(200):
            c = deepcopy(self)
            del_idx = randrange(self._cls_size)
            c._del(del_idx)
            if del_idx < self._nonaxis_num*self._order:
                if c._add_nonaxis_point():
                    for i in range(c._order):
                        if c.atoms[i].bond(self.atoms[del_idx]):
                            break
                    else:
                        self.deepcopy(c)
                        return list(range(c._order))
            else:
                if c._add_axis_point():
                    if not c.atoms[-1].bond(self.atoms[del_idx]):
                        self.deepcopy(c)
                        return [self._cls_size-1]
        else:
            utils.write_log('mutation_one failure.')
            self.perturbation()
        return []

    def squeez(self):
        #删除厚度超过2的原子
        i = 0
        while i < self._nonaxis_num*self._order:
            for a in self.atoms[i+self._order:-self._order]:
                if not a.plane().bond(self.atoms[i].plane()):
                    continue
                if a.dist(self.atoms[i]) > 2*tolerance_max*self._bond_len:
                    if abs(self.atoms[i].z) >= abs(a.z):
                        self._del(self.atoms[i])
                    else:
                        self._del(a)
                    break
            else:
                i += self._order
        if self.get_size() == self._cls_size:
            self.perturbation()
            return

        self.set_energy(0.)
        #把这些原子再加上
        while self.get_size() < self._cls_size:
            while True:
                nb = choice(self.atoms)
                modify = Atom(nb.plane() + Point.random_plane(uniform(tolerance_min,1)*self._bond_len) + \
                    Point(0,0,nb.z+uniform(-0.5,0.5)), self._element)
                if self._cls_size-self.get_size() < self._order:
                    self._put_axis(modify)
                    break
                elif self._is_axis(modify):
                    continue
                for a in self.atoms:
                    if not modify.is_legal(a) or \
                        modify.plane().dist(a.plane()) < 0.5*self._bond_len:
                        break
                else:
                    self._add(modify)
                    break

    def change_naxis(self):
        '''增减对称轴上的点
调用：_add_axis_point，_add_nonaxis_point'''
        best = deepcopy(self)
        best.set_energy(0.)
        if isinstance(self, SymCsCluster):
            increase_probability = 0.5 if self._axis_num >= self._order else 1.
        else:
            increase_probability = 5 ** (-(self._axis_num//self._order))
        for _ in range(10):
            copy_ = deepcopy(self)
            if random() <= increase_probability: #增对称轴上的点
                copy_._del(randrange(copy_._nonaxis_num))
                for i in range(copy_._order):
                    if not copy_._add_axis_point():
                        break
                else:
                    #类别降为Cluster，这样能避免调用_calc_axis_num
                    cc = Cluster(copy_.atoms)
                    fix_coord = [7]*self.get_size()
                    for i in range(copy_._order):
                        fix_coord[self._cls_size-i-1] = 0
                    optimize(cc, fix_coord=fix_coord)
                    if cc.get_energy() < best.get_energy():
                        best.atoms = deepcopy(cc.atoms)
                        best.set_energy(cc.get_energy())
                        for i in range(best._order):
                            best._put_axis(best.atoms[self._cls_size-i-1])
            else: #减少对称轴上的点
                for i in range(copy_._order):
                    copy_._del(randrange(copy_.get_size()-copy_._axis_num, copy_.get_size()))
                if copy_._add_nonaxis_point():
                    #类别降为Cluster，这样能避免调用_calc_axis_num
                    cc = Cluster(copy_.atoms)
                    fix_coord = [7]*self.get_size()
                    for i in range(copy_._order):
                        fix_coord[i] = 0
                    optimize(cc, fix_coord=fix_coord)
                    if cc.get_energy() < best.get_energy():
                        best.atoms = deepcopy(cc.atoms)
                        best.set_energy(cc.get_energy())
                        added = cc.atoms[0]
                        best._del(0)
                        best._add(added)
        self.deepcopy(best)


class SymCsCluster(SymCluster):
    '''处理Cs对称的团簇
    对称面为xy平面'''
    _order = 2 #对称群阶数
    sym = 'Cs'

    def _is_axis(self, a):
        '''判断p是否在对称面上'''
        return abs(a.z) < tolerance_min/2 * self._bond_len

    def _put_axis(self, a):
        '''将p修改为在对称面上'''
        a.z = 0

    def symmetry_atom(self, a):
        return [a, Atom(a.x, a.y, -a.z, a.elem)]


class SymC2Cluster(SymCluster):
    '''处理C2对称的团簇,绕z轴旋转180°对称'''
    _order = 2 #对称群阶数
    sym = 'C2'

    def _is_axis(self, a):
        '''判断p是否在对称轴上'''
        return hypot(a.x, a.y) < tolerance_min/2 * self._bond_len

    def symmetry_atom(self, a):
        return [a, Atom(-a.x, -a.y, a.z, a.elem)]


class SymC3Cluster(SymCluster):
    '''处理C3对称的团簇
    对称轴为z轴'''
    _order = 3 #对称群阶数
    sym = 'C3'

    def _is_axis(self, a):
        '''判断p是否在对称轴上'''
        return hypot(a.x, a.y) < tolerance_min / sqrt(3) * self._bond_len

    def symmetry_atom(self, a):
        return [a, Atom(a.rotate(2*pi/3, Point(0,0,1)),a.elem), \
                Atom(a.rotate(4*pi/3, Point(0,0,1)),a.elem)]


class SymC5Cluster(SymCluster):
    '''处理C5对称的团簇
    对称轴为z轴'''
    _order = 5
    sym = 'C5'

    def _is_axis(self, a):
        '''判断p是否在对称轴上'''
        return hypot(a.x, a.y) < tolerance_min / (2*sin(0.2*pi)) * self._bond_len

    def symmetry_atom(self, a):
        return [a, Atom(a.rotate(2*pi/5, Point(0,0,1)), a.elem), \
                Atom(a.rotate(4*pi/5, Point(0,0,1)), a.elem), \
                Atom(a.rotate(6*pi/5, Point(0,0,1)), a.elem),
                Atom(a.rotate(8*pi/5, Point(0,0,1)), a.elem)]


class SymCiCluster(SymCluster):
    '''处理C3对称的团簇
    对称轴为z轴'''
    _order = 2 #对称群阶数
    sym = 'Ci'

    def _is_axis(self, a):
        '''判断p是否在对称轴上'''
        return 2*a.norm() < tolerance_min * self._bond_len

    def symmetry_atom(self, a):
        return [a, Atom(-a.x, -a.y, -a.z, a.elem)]

    def random(self):
        '''产生随机结构'''
        if self._cls_size%2: #奇原子，则中心有原子
            self.atoms = [Atom(Point(0), self._element)]
            self._axis_num = 1
        else:
            self.atoms = []
            self._axis_num = 0
        self._nonaxis_num = 0
        radius = self._cls_size**(1./3) / 2 * self._bond_len
        while self.get_size() < self._cls_size//2*2:
            for _ in range(100): #每个原子尝试100次
                now = Atom(Point.random(random()*radius), self._element)
                if self._is_axis(now): #这里不同
                    continue
                #判断新产生的原子是否合法
                if not self.is_legal(now):
                    continue
                #将新产生的点放入坐标
                self._add(now)
                break
            else: #尝试100次没有找到第i个原子的合适位置
                radius += 0.1
                if radius > self._cls_size**(1./3) * 4 * self._bond_len:
                    raise
                if self._cls_size%2: #奇原子，则中心有原子
                    self.atoms = [Atom(Point(0), self._element)]
                    self._axis_num = 1
                else:
                    self.atoms = []
                    self._axis_num = 0
                self._nonaxis_num = 0

    def mating(self, fater, mother):
        '''交叉'''
        f = deepcopy(fater)
        m = deepcopy(mother)
        if f._axis_num == 1:
            f.atoms.pop()
            m.atoms.pop()
        f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

        for _ in range(100):
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            cut_angle = random() * pi
            self.atoms = []
            for a in f.atoms:
                if a.z > 0 and acos(a.x / hypot(a.x, a.y)) <= cut_angle:
                    self._add(a)
            for a in m.atoms:
                if a.z > 0 and acos(a.x / hypot(a.x, a.y)) > cut_angle:
                    self._add(a)
            if f._axis_num == 1:
                self.add_atom(Point(0), self._element)
                self._axis_num = 1
            if self.get_size() == self._cls_size:
                break
        else:
            print('mating failure.')
            f.add_atom(Point(0), self._element)
            self.deepcopy(f)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

            self._nonaxis_num -= 1

class SymTCluster(SymCluster):
    _order = 12 #对称群阶数
    sym = 'T'

    def symmetry_atom(self, a):
        return [a, Atom(a.y, a.z, a.x, a.elem), Atom(a.z, a.x, a.y, a.elem),
                Atom(-a.x, -a.y, a.z, a.elem), Atom(-a.y, -a.z, a.x, a.elem), Atom(-a.z, -a.x, a.y, a.elem),
                Atom(-a.x, a.y, -a.z, a.elem), Atom(-a.y, a.z, -a.x, a.elem), Atom(-a.z, a.x, -a.y, a.elem),
                Atom(a.x, -a.y, -a.z, a.elem), Atom(a.y, -a.z, -a.x, a.elem), Atom(a.z, -a.x, -a.y, a.elem)]

    def _is_axis(self, a):
        '''判断p是否在对称轴上
        返回：0不在对称轴上，1在C3轴上，2在D2轴上。不允许在原点'''
        if not Atom(abs(a.x), abs(a.y), abs(a.z), a.elem).is_legal(Atom(abs(a.y), abs(a.z), abs(a.x), a.elem)):
            return 1
        if a.is_legal(Atom(-a.x, -a.y, a.z, a.elem)) and \
            a.is_legal(Atom(-a.x, a.y, -a.z, a.elem)) and \
            a.is_legal(Atom(a.x, -a.y, -a.z, a.elem)):
            return 0
        return 2

    def _put_axis(self, a):
        if not a.is_legal(Atom(-a.x, -a.y, a.z, a.elem)):
            a.x = a.y = 0
        elif not a.is_legal(Atom(-a.x, a.y, -a.z, a.elem)):
            a.x = a.z = 0
        elif not a.is_legal(Atom(a.x, -a.y, -a.z, a.elem)):
            a.y = a.z = 0
        else: #默认放到C3轴上
            a.set_coordinate(Point(1,1,1)*max(abs(a.x),abs(a.y),abs(a.z)))

    def _add(self, a):
        '''增加一个原子
        非对称轴原子加在前端，对称轴原子加在后端'''
        axis = self._is_axis(a)
        if axis == 1:
            self.add_atom(a)
            self.add_atom(Atom(-a.x, -a.y, a.z, a.elem))
            self.add_atom(Atom(-a.x, a.y, -a.z, a.elem))
            self.add_atom(Atom(a.x, -a.y, -a.z, a.elem))
            self._axis_num_C3 += 1
        elif axis == 2:
            if a.x != 0:
                c = a.x
            elif a.y != 0:
                c = a.y
            else:
                c = a.z
            self.add_atom(Atom(c, 0, 0, a.elem))
            self.add_atom(Atom(-c, 0, 0, a.elem))
            self.add_atom(Atom(0, c, 0, a.elem))
            self.add_atom(Atom(0, -c, 0, a.elem))
            self.add_atom(Atom(0, 0, c, a.elem))
            self.add_atom(Atom(0, 0, -c, a.elem))
            self._axis_num_D2 += 1
        else:
            self.atoms = self.symmetry_atom(a) + self.atoms
            self._nonaxis_num += 1

    def random(self):
        '''产生随机结构'''
        self.atoms = []
        self._axis_num_C3 = self._axis_num_D2 = self._nonaxis_num = 0
        radius = 0.4 * self._cls_size**(1./3) * self._bond_len
        while self.get_size() < self._cls_size:
            for _ in range(100): #每个原子尝试100次
                now = Atom(Point.random(uniform(0.8*radius, radius)), self._element)
                if now.norm() < 0.5 * self._bond_len:
                    continue
                if self._is_axis(now): #离对称轴太近,就放到对称轴上
                    self._put_axis(now)
                #未产生原子个数已不足以产生一个非轴原子,则让该原子在对称轴上
                if self._cls_size - self.get_size() < self._order:
                    self._put_axis(now)
                #判断新产生的原子是否合法
                if not self.is_legal(now):
                    continue
                #将新产生的点放入坐标
                self._add(now)
                break
            else: #尝试100次没有找到第i个原子的合适位置
                radius *= 1.05
                self.atoms = []
                self._axis_num_C3 = self._axis_num_D2 = self._nonaxis_num = 0


class Bilayer(SymCsCluster):

    def random(self):
        a = sqrt(self._cls_size/2) * self._bond_len * 0.1
        b = uniform(1., 3) * a
        while self.get_size() < self._cls_size:
            for _ in range(20):
                atom = Atom(Point(random()*a, random()*b, random()*self._bond_len), self._element)
                if atom.z < tolerance_min/2 * self._bond_len or self.get_size()==self._cls_size-1:
                    atom.z = 0.
                if self.is_legal(atom):
                    self._add(atom)
                    break
            else:
                a *= 1.1
                b *= 1.1
                self.atoms = []

    def perturbation(self):
        super(Bilayer, self).perturbation(0.08)
        for i in range(self._nonaxis_num):
            if abs(self.atoms[i].z) > self._bond_len:
                self.atoms[i].z = self._bond_len
                self.atoms[i+self._nonaxis_num].z = -self._bond_len

    def _add_nonaxis_point(self):
        '''增加一个非对称轴上的点
        返回True或False'''
        for _ in range(100):
            new_atom = choice(self.atoms) + Point.random(self._bond_len)
            if self._is_axis(new_atom):
                continue
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1 and \
                abs(new_atom.z) < self._bond_len: #这里与基类不同
                self._add(new_atom)
                return True
        else:
            return False


class Planar(MonoCluster):
    '''平面单质团簇'''

    def random(self):
        '''产生随机结构'''
        self.atoms = []
        w = sqrt(pi*self._cls_size/2) / 2
        while self.get_size() < self._cls_size:
            for _ in range(20):
                a = Atom(Point(random()*w, random()*w, random()*2), self._element)
                thickness = sum([1 if a.plane().dist(q.plane()) < 0.5 else 0 for q in self.atoms])
                if thickness > 1:
                    continue
                for i in range(self.get_size()):
                    if a.dist(self.atoms[i]) < tolerance_min:
                        break
                else:
                    self.add_atom(a)
                    break
            else:
                w += 0.1
                self.atoms = []
        self.zoom(self._bond_len)

    def mating(self, father, mother):
        '''平面结构用的交叉.
        先将双亲结构旋转到x轴为最长方向,z轴为最短方向;
        然后用任意平行z轴的平面切割父代和母代;
        然后取父代的半边,母代的半边,合并'''
        father.place()
        father.rotate(0, 0, 2*pi*random())
        father.rotate(0, pi/2, 0)
        mother.place()
        mother.rotate(0, 0, 2*pi*random())
        mother.rotate(0, pi/2, 0)
        locus = randrange(self._cls_size/10+2, self._cls_size/2+3)
        self._merge(father.atoms[:locus], mother.atoms[locus:])

    def perturbation(self, step=FixCluster.amplitude):
        '''给结构微小扰动
        平均每个原子扰动100次,每次扰动幅度不超过step'''
        self.place()
        for _ in range(self._cls_size*200):
            idx = randrange(self.get_size())
            new_pos = self.atoms[idx] + Point.random(step)
            thickness = 0
            for i in range(self.get_size()):
                if i == idx:
                    continue
                #两原子不可太近
                if not new_pos.is_legal(self.atoms[i]):
                    break
                #厚度不得增加
                if self.atoms[idx].plane().dist(self.atoms[i].plane()) < 0.5 * self._bond_len:
                    thickness += 1
                elif new_pos.plane().dist(self.atoms[i].plane()) < 0.5 * self._bond_len and thickness > 0:
                    break
            else:
                self.atoms[idx] = new_pos

    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for _ in range(100):
            src = choice(self.atoms)
            src_bak = deepcopy(src)
            dest = choice(self.atoms)
            src = dest + Point.random(self._bond_len)
            thickness = sum([1 if a.plane().dist(dest.plane()) < 0.5 else 0 for a in self.atoms])
            if src.dist(src_bak) > 1.5 * self._bond_len and \
                self.is_legal(src) and self.get_degree(src) > 1 and thickness < 2: #找到i点附近合适的放置位置
                return self.atoms.index(src)
            else:
                src = src_bak
        else:
            return -1

    def mutation_one_better(self):
        self.place()
        super(Planar, self).mutation_one_better()

    def squeez(self):
        self.place()
        thick = []
        for a in self.atoms:
            same_unit = [q for q in self.atoms if a.plane().dist(q.plane()) < 0.5*self._bond_len]
            thick.append(max([a.z for a in same_unit]) - min([a.z for a in same_unit]))
        thickness = max(thick)
        if thickness > 1.5 * self._bond_len:
            self.set_energy(0.)
        for a in self.atoms:
            a.z *= tolerance_max * self._bond_len / thickness
        self.clean()


###############################################################################
#                   二元团簇
###############################################################################
class BiCluster(FixCluster):
    '''二元团簇类'''

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(BiCluster,cls).config()
        e = cfg.get('cluster', 'element').split()
        if len(e) != 2:
            print('Please input two elements in config.ini.')
            os._exit(1)
        cls.set_element(e[0], e[1])
        atom_num = cfg.get('cluster', 'atom_num').split()
        if len(atom_num) != 2:
            print('Please input two atom num in config.ini.')
            os._exit(1)
        cls.set_atom_num(int(atom_num[0]), int(atom_num[1]))

    @classmethod
    def set_element(cls, element1, element2):
        if get_element_id(element1) == -1 or get_element_id(element2) == -1:
            print('input element error.')
            os._exit(1)
        cls._element1, cls._element2 = element1, element2
        cls._bond_len11 = BOND(cls._element1, cls._element1)
        cls._bond_len12 = BOND(cls._element1, cls._element2)
        cls._bond_len22 = BOND(cls._element2, cls._element2)

    @classmethod
    def set_atom_num(cls, atom_num1, atom_num2):
        cls._atom_num1, cls._atom_num2 = atom_num1, atom_num2
        cls._cls_size = atom_num1 + atom_num2

    @classmethod
    def get_element(cls):
        return cls._element1, cls._element2

    @classmethod
    def get_atom_num(cls):
        return cls._cls_size

    def check(self):
        '''检查结构是否正确
        '''
        for a in self.atoms[:self._atom_num1]:
            if a.elem != self._element1:
                return False
        for a in self.atoms[self._atom_num1:]:
            if a.elem != self._element2:
                return False
        return True

    @classmethod
    def get_size1(cls):
        '''被调用: Cluster.write_poscar'''
        return cls._atom_num1

    @classmethod
    def get_size2(cls):
        return cls._atom_num2

    def _merge(self, down1, down2, up1, up2, horizon_adjust=True, merge=True):
        '''将down1, down2, up1, up2拼成起来,并计算两个半球拼起来的合适gap
        up1, up2可以为空
        该函数将保持down1, down2不动，up1, up2的最靠下原子朝down1, down2的最靠上原子附近方向移动
        horizon_adjust: 是否进行水平调整。通常进行调整能更好的拼合，但在SymBiCluster中要保持对象，故为False
        merge为False，则只计算gap，不拼合'''
        if len(down1)+len(down2) == 0:
            if merge:
                self.atoms = deepcopy(up1+up2)
            return Point(0.)
        upper = max(self._bond_len11, self._bond_len12, self._bond_len22)
        if horizon_adjust:
            oa = reduce(Point.__add__, down1+down2) / len(down1+down2)
            ob = reduce(Point.__add__, up1+up2) / len(up1+up2)
            ort = (oa - ob).unit()
            upper = (oa - ob).norm()
        else:
            ort = Point(0,0,1)
            upper = max([a.z for a in down1+down2]) - min([a.z for a in up1+up2]) + 5
        lower = 0.
        #二分查找
        norm = (lower+upper)/2
        while lower < upper - 0.1:
            try:
                if len(down1) != 0 and len(up1) != 0:
                    for p in down1:
                        for q in up1:
                            if not p.is_legal(q + ort*norm):
                                raise
                if len(down1) != 0 and len(up2) != 0:
                    for p in down1:
                        for q in up2:
                            if not p.is_legal(q + ort*norm):
                                raise
                if len(down2) != 0 and len(up1) != 0:
                    for p in down2:
                        for q in up1:
                            if not p.is_legal(q + ort*norm):
                                raise
                if len(down2) != 0 and len(up2) != 0:
                    for p in down2:
                        for q in up2:
                            if not p.is_legal(q + ort*norm):
                                raise
            except:
                upper = norm
            else:
                lower = norm
            norm = (lower+upper)/2

        #将两个半球拼在一起
        if merge:
            self.atoms = deepcopy(down1 + [a+ort*norm for a in up1] + \
                                  down2 + [a+ort*norm for a in up2])
        return ort*norm

    def random(self, forbidden1=False, forbidden2=False):
        '''产生随机结构'''
        a,b,c = abs(random())+0.1, abs(random())+0.1, abs(random()+0.1)
        d = (self._cls_size/a/b/c)**(1./3)
        n_atom1 = n_atom2 = 0
        while n_atom1 + n_atom2 < self._cls_size:
            for _ in range(50):
                elem = self._element1 if randrange(self._cls_size-n_atom1-n_atom2)< \
                    self._atom_num1-n_atom1 else self._element2
                atom = Atom(Point(random()*a, random()*b, random()*c)*d, elem)
                try:
                    if not self.is_legal(atom):
                        raise
                    if forbidden1 and elem == self._element1:
                        for i in range(n_atom1):
                            if atom.bond(self.atoms[i]):
                                raise
                    if forbidden2 and elem == self._element2:
                        for i in range(n_atom1, n_atom1+n_atom2):
                            if atom.bond(self.atoms[i]):
                                raise
                except:
                    continue
                else:
                    if elem == self._element1:
                        self.atoms.insert(0,atom)
                        n_atom1 += 1
                    else:
                        self.add_atom(atom)
                        n_atom2 += 1
                    break
            else:
                d += 0.1
                self.atoms = []
                n_atom1 = n_atom2 = 0

    def mating(self, father, mother):
        '''杂交,沿z=face平面杂交'''
        f = deepcopy(father)
        m = deepcopy(mother)
        f.center()
        m.center()

        for _ in range(1000):
            f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

            face = uniform(-self._bond_len11/2, self._bond_len11/2) #杂交面
            up1 = [a for a in f.atoms[:self._atom_num1] if a.z>=face]
            up2 = [a for a in f.atoms[self._atom_num1:self._cls_size] if a.z>=face]
            down1 = [a for a in m.atoms[:self._atom_num1] if a.z<face]
            down2 = [a for a in m.atoms[self._atom_num1:self._cls_size] if a.z<face]
            if len(up1)+len(up2) == 0 or len(down1)+len(down2) == 0:
                continue
            if len(up1)+len(down1) == self._atom_num1 and \
                len(up2)+len(down2) == self._atom_num2:
                self._merge(down1, down2, up1, up2)
                break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    def mating_principal_direction(self, father, mother):
        '''先将双亲结构旋转到x轴为最长方向,z轴为最短方向;
        然后用任意垂直于x轴的平面切割父代和母代;
        然后随机取父代的A半边,母代的A半边,合并'''
        f = deepcopy(father)
        m = deepcopy(mother)
        f.place()
        if random() > 0.5:#选用父代的哪半边
            for i in range(self._cls_size):
                f.atoms[i].x = -f.atoms[i].x
        m.place()
        if random() > 0.5:#选用母代的哪半边
            for i in range(self._cls_size):
                m.atoms[i].x = -m.atoms[i].x

        #普通方式交叉
        for _ in range(1000):
            rot_x = 2*pi*random()
            rot_y = 2*pi*random()
            rot_z = 2*pi*random()
            f.rotate(rot_x, rot_y, rot_z)
            m.rotate(rot_x, rot_y, rot_z)

            face = uniform(-self._bond_len11/2, self._bond_len11/2) #杂交面
            up1 = [a for a in f.atoms[:self._atom_num1] if a.z>=face]
            up2 = [a for a in f.atoms[self._atom_num1:self._cls_size] if a.z>=face]
            down1 = [a for a in m.atoms[:self._atom_num1] if a.z<face]
            down2 = [a for a in m.atoms[self._atom_num1:self._cls_size] if a.z<face]
            if len(up1)+len(down1) == self._atom_num1 and \
                len(up2)+len(down2) == self._atom_num2:
                self.atoms = up1 + down1 + up2 + down2
                break
        else:
            print('mating failure')
            self.deepcopy(f)
            self.perturbation() #改用变异
            return

    def exchange(self):
        '''交换两个元素.'''
        best_energy = 0.
        ex1 = randrange(self._atom_num1)
        for i in sample(list(range(self._atom_num1)), self._atom_num1//4+1):
            c = deepcopy(self)
            c.atoms[i].elem = self._element2
            optimize(c, fix_coord=[7]*i+[0]+[7]*(self.get_size()-1-i))
            if c.get_energy() <= best_energy:
                best_energy = c.get_energy()
                ex1 = i

        best_energy = 0.
        ex2 = randrange(self._atom_num1, self._cls_size)
        for i in sample(list(range(self._atom_num1, self._cls_size)), self._atom_num2//4+1):
            c = deepcopy(self)
            c.atoms[i].elem = self._element1
            optimize(c, fix_coord=[7]*i+[0]+[7]*(self.get_size()-1-i))
            if c.get_energy() <= best_energy:
                best_energy = c.get_energy()
                ex2 = i

        self.atoms[ex1].exchange_coordinate(self.atoms[ex2])

    def exchange_atom_better(self):
        '''交换两元素，尝试多次'''
        best = deepcopy(self)
        best.set_energy(0.)
        c = deepcopy(self)
        for _ in range(randrange(max(self._atom_num1, self._atom_num2)//3+1)+2): #交换多次
            c = deepcopy(self)
            ex1 = randrange(self._atom_num1)
            ex2 = randrange(self._atom_num1, self._cls_size)
            c.atoms[ex1].exchange_coordinate(c.atoms[ex2])
            fix_coord = [7]*self.get_size()
            fix_coord[ex1] = fix_coord[ex2] = 0
            optimize(c, fix_coord=fix_coord)
            if c.get_energy() < best.get_energy():
                best = deepcopy(c)
        self.deepcopy(best)


class SymBiCluster(BiCluster, FixSymCluster):
    '''对称团簇的抽象类,不能直接使用
    具体类中需实现_is_axis, _add, _adjust
    存储方式为:对称轴上原子1，非对称轴上原子1，非称轴上原子2，对称轴上原子2'''

    def __init__(self):
        super(SymBiCluster, self).__init__(self)
        self._nonaxis_num1 = 0 #不在对称轴上的点的个数
        self._axis_num1 = 0 #在对称轴上的点的个数
        self._nonaxis_num2 = 0 #不在对称轴上的点的个数
        self._axis_num2 = 0 #在对称轴上的点的个数

    def get_natom1(self):
        return self._order * self._nonaxis_num1 + self._axis_num1

    def get_natom2(self):
        return self._order * self._nonaxis_num2 + self._axis_num2

    def _arrange(self):
        '''每次从别的地方(例如读recover)读取结构,都要调用该函数调整原则顺序，并给_nonaxis_num和axis_num赋值
        dmol优化成高对称性后会改变原子排列顺序！
        同SymCluster'''
        self._set_symmetry_axis()

        self._axis_num1 = self._nonaxis_num1 = self._axis_num2 = self._nonaxis_num2 = 0
        i = 0
        while i < self._cls_size - self._axis_num2:
            if self._is_axis(self.atoms[i]):
                if i < self._atom_num1:
                    self._axis_num1 += 1
                    self.atoms.insert(0,self.atoms[i])
                    i += 1
                    del self.atoms[i]
                else:
                    self._axis_num2 += 1
                    self.add_atom(self.atoms[i])
                    del self.atoms[i]
            else:
                if i < self._atom_num1:
                    self._nonaxis_num1 += 1
                else:
                    self._nonaxis_num2 += 1
                sym_atoms = self.symmetry_atom(self.atoms[i])
                for k in range(1, len(sym_atoms)):
                    for j in range(i+1, self._cls_size):
                        if (self.atoms[j]-sym_atoms[k]).norm() < 0.1:
                            self.atoms[i+k],self.atoms[j] = self.atoms[j],self.atoms[i+k]
                            break
                    else:
                        print('cluster loose symmetry.')
                        raise utils.SymmetryErr('cluster loose symmetry.')
                i += self._order

    def isomorphic(self, other):
        '''判断两个团簇是否同构
        参数other: 团簇类型,要比较的团簇
        返回: True或False
        说明: 使用前确保_inertia和_inertia2已经计算.'''
        if self._axis_num1 != other._axis_num1 or self._axis_num2 != other._axis_num2:
            return False
        return Cluster.isomorphic(self, other)

    def _add(self, a):
        '''增加一个原子，对该原子的合法性不做检查。
        被调用:random,perturbation'''
        if self._is_axis(a):
            if a.elem == self._element1:
                self.atoms.insert(0, a)
                self._axis_num1 += 1
            else:
                self.add_atom(a)
                self._axis_num2 += 1
        else:
            if a.elem == self._element1:
                self.atoms = self.atoms[:self._axis_num1] + self.symmetry_atom(a) + self.atoms[self._axis_num1:]
                self._nonaxis_num1 += 1
            else:
                self.atoms = self.atoms[:self.get_natom1()] + self.symmetry_atom(a) + self.atoms[self.get_natom1():]
                self._nonaxis_num2 += 1

    def _del(self, a):
        '''删除一个原子
        参数a：Atom类型或数组下标
        被调用:perturbation'''
        if isinstance(a,Atom): #a为点
            index = self.atoms.index(a)
        else: #p为点的下标
            index = a
            a = self.atoms[index]
        if self._is_axis(a): #删除对称轴上的原子
            del self.atoms[index]
            if a.elem == self._element1:
                self._axis_num1 -= 1
            else:
                self._axis_num2 -= 1
        else: #删除非对称轴上的原子
            if a.elem == self._element1:
                index = self._axis_num1 + (index-self._axis_num1) // self._order * self._order
                del self.atoms[index:index+self._order]
                self._nonaxis_num1 -= 1
            else:
                index = self.get_natom1() + (index-self.get_natom1()) // self._order * self._order
                del self.atoms[index:index+self._order]
                self._nonaxis_num2 -= 1

    def _add_axis_point(self, elem):
        '''增加一个对称轴上的点
        返回True或False'''
        z_min = min([a.z for a in self.atoms]) - max(self._bond_len11,self._bond_len12,self._bond_len22)
        z_max = max([a.z for a in self.atoms]) + max(self._bond_len11,self._bond_len12,self._bond_len22)
        diameter = self.diameter()
        for _ in range(100):
            if isinstance(self, SymCsBiCluster):
                new_atom = Atom(Point.random_plane(diameter), elem)
            else:
                new_atom = Atom(Point(0,0,uniform(z_min,z_max)), elem)
            #是否为孤立点
            for a in self.atoms:
                if a.bond(new_atom):
                    break
            else:
                continue
            #是否合法
            if self.is_legal(new_atom):
                self._add(new_atom)
                return True
        else:
            return False

    def _add_nonaxis_point(self, elem):
        '''增加一个非对称轴上的点
        返回True或False'''
        for _ in range(100):
            neighbour = choice(self.atoms)
            new_atom = neighbour+Point.random(BOND(neighbour.elem, elem))
            new_atom.elem = elem
            if self._is_axis(new_atom):
                continue
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1:
                self._add(new_atom)
                return True
        else:
            return False

    def random(self):
        '''产生随机结构, 在圆柱区域内'''
        self.atoms = []
        self._axis_num1 = self._nonaxis_num1 = self._axis_num2 = self._nonaxis_num2 = 0
        aver_bond = (self._bond_len11+self._bond_len12+self._bond_len22)/3
        height = max(uniform(0.5,1.5) * self._cls_size**(1.0/3), \
            (self._atom_num1%self._order+self._atom_num2%self._order)) * aver_bond
        radius = sqrt(self._cls_size*3/4*aver_bond**3 / height) * 0.7
        while self.get_size() < self._cls_size:
            if random() <= float(self._atom_num1-self.get_natom1()) / \
                (self._cls_size - self.get_natom1() - self.get_natom2()):
                atom_type = self._element1
            else:
                atom_type = self._element2
            for _ in range(100):
                now = Atom(Point.random_plane(random()*radius) + Point(0,0,random()*height), atom_type)
                #离对称轴太近,就放到对称轴上
                if self._is_axis(now):
                    self._put_axis(now)
                #未产生原子个数已不足以产生一个非轴原子,则让该原子在对称轴上
                if atom_type == self._element1 and self._atom_num1-self.get_natom1() < self._order or \
                    atom_type == self._element2 and self._atom_num2-self.get_natom2() < self._order:
                    self._put_axis(now)
                #判断新产生的原子是否合法
                if self.is_legal(now):
                    self._add(now)
                    break
            else: #产生第i个原子失败
                radius *= 1.05
                height *= 1.05
                self.atoms = []
                self._axis_num1 = self._nonaxis_num1 = self._axis_num2 = self._nonaxis_num2 = 0

    def mating(self, father, mother):
        f = deepcopy(father)
        m = deepcopy(mother)
        f.center()
        m.center()
        f.rotate(0, 0, 2*pi*random())
        m.rotate(0, 0, 2*pi*random())
        if isinstance(self, SymCsBiCluster):
            f.rotate(pi/2, 0, 0)
            m.rotate(pi/2, 0, 0)
        else:
            if random() > 0.5: #选用母代的哪半边
                for i in range(self._cls_size):
                    m.atoms[i].z = -m.atoms[i].z

        for _ in range(100):
            cut_plane_up = uniform(0, max([a.z for a in f.atoms])*0.7) #杂交面
            up1 = [a for a in f.atoms[f._axis_num1:self._atom_num1] if a.z>=cut_plane_up]
            up_axis1 = [a for a in f.atoms[:f._axis_num1] if a.z>=cut_plane_up]
            up2 = [a for a in f.atoms[self._atom_num1:
                                          self._atom_num1+f._nonaxis_num2*f._order] if a.z>=cut_plane_up]
            up_axis2 = [a for a in f.atoms[self._atom_num1+f._nonaxis_num2*f._order:] if a.z>=cut_plane_up]

            lower = min([a.z for a in m.atoms])
            upper = max([a.z for a in m.atoms])
            while lower < upper - 0.1: #二分查找上边界
                cut_plane_down = (lower+upper)/2
                down1 = [a for a in m.atoms[m._axis_num1:self._atom_num1] if a.z<cut_plane_down]
                down_axis1 = [a for a in m.atoms[:m._axis_num1] if a.z<cut_plane_down]
                defect_num1 = len(up1) + len(down1) + len(up_axis1) + len(down_axis1) - self._atom_num1
                down2 = [a for a in m.atoms[self._atom_num1: self._atom_num1+m._nonaxis_num2*m._order] if a.z<cut_plane_down]
                down_axis2 = [a for a in m.atoms[self._atom_num1+m._nonaxis_num2*m._order:] if a.z<cut_plane_down]
                defect_num2 = len(up2) + len(down2) + len(up_axis2) + len(down_axis2) - self._atom_num2

                if len(up1) + len(down1) > self._atom_num1 or \
                    len(up2) + len(down2) > self._atom_num2:
                    upper = cut_plane_down
                elif defect_num1 > 1:
                    if defect_num2 < -1:
                        break
                    else:
                        upper = cut_plane_down
                elif defect_num1 < -1:
                    if defect_num2 > 1:
                        break
                    else:
                        lower = cut_plane_down
                elif defect_num2 > 1:
                    upper = cut_plane_down
                elif defect_num2 < -1:
                    lower = cut_plane_down
                else:
                    break

            if len(up1) + len(down1) > self._atom_num1 or \
                len(up2) + len(down2) > self._atom_num2 or \
                abs(defect_num1) > 1 or abs(defect_num2) > 1:
                continue

            self._nonaxis_num1 = (len(up1) + len(down1)) // self._order
            self._axis_num1 = len(down_axis1) + len(up_axis1)
            self._nonaxis_num2 = (len(up2) + len(down2)) // self._order
            self._axis_num2 = len(down_axis2) + len(up_axis2)

            axis_atom1 = up_axis1 + down_axis1
            axis_atom2 = up_axis2 + down_axis2

            #将两个半球拼起来
            move = self._merge(down1+down_axis1, down2+down_axis2, up1+up_axis1, up2+up_axis2, False, False)
            for a in up_axis1:
                a += move
            axis1 = down_axis1 + up_axis1
            for a in up1:
                a += move
            naxis1 = down1 + up1
            for a in up2:
                a += move
            naxis2 = down2 + up2
            for a in up_axis2:
                a += move
            axis2 = down_axis2 + up_axis2
            self.atoms = deepcopy(axis1 + naxis1 + naxis2 + axis2)

            if isinstance(self, SymCsBiCluster):
                self.rotate(-pi/2, 0, 0)
                self.center()

            #增减对称轴原子使总原子数不变
            if defect_num1 == -1: #原子数缺1
                if not self._add_axis_point(self._element1):
                    continue
            elif defect_num1 == 1: #原子数多1
                del self.atoms[randrange(len(axis_atom1))]
                self._axis_num1 -= 1
            if defect_num2 == -1: #原子数缺1
                if not self._add_axis_point(self._element2):
                    continue
            elif defect_num2 == 1: #原子数多1
                del self.atoms[-randrange(len(axis_atom2))-1]
                self._axis_num2 -= 1

            break
        else:
            print('mating failure.')
            if isinstance(self, SymCsBiCluster):
                f.rotate(-pi/2, 0, 0)
            self.deepcopy(f)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for _ in range(100):
            del_idx = randrange(self._cls_size)
            old_atom = self.atoms[del_idx]
            self._del(del_idx)
            new_idx = None
            if self._is_axis(old_atom):
                if self._add_axis_point(old_atom.elem):
                    new_idx = 0 if old_atom.elem == self._element1 else -1
            else:
                if self._add_nonaxis_point(old_atom.elem):
                    if old_atom.elem == self._element1:
                        new_idx = self._axis_num1
                    else:
                        new_idx = self._atom_num1
            if new_idx is None:
                self._add(old_atom)
                continue
            elif self.get_degree(self.atoms[new_idx]) > 1:
                if new_idx != 0 and new_idx != -1:
                    for i in range(self._order):
                        if self.atoms[new_idx+i].bond(old_atom):
                            break
                    else:
                        return [new_idx+i*self._order for i in range(self._order)]
                elif self.atoms[new_idx].bond(old_atom):
                    return [new_idx]
            self._del(new_idx)
            self._add(old_atom)
        utils.write_log('mutation_one failure.')
        self.perturbation()
        return []

    def exchange(self):
        if self._nonaxis_num1 != 0 and self._nonaxis_num2 != 0:
            list1 = list(range(self._axis_num1, self._atom_num1))
            list2 = list(range(self._atom_num1, self._cls_size-self._axis_num2))
        elif self._axis_num1 != 0 and self._axis_num2 != 0:
            list1 = list(range(self._axis_num1))
            list2 = list(range(self._cls_size-self._axis_num2, self._cls_size))
        else:
            print('exchange failure')
            self.perturbation()
            return
        shuffle(list1)
        shuffle(list2)
        n_it = min(len(list1), len(list2)) // 4 + 1
        best_energy = 0.
        for i in list1[:n_it]:
            c = deepcopy(self)
            c.atoms[i].elem = self._element2
            cc = Cluster(c.atoms) #类别降为Cluster
            optimize(cc, fix_coord=[7]*i+[0]+[7]*(self.get_size()-1-i))
            if cc.get_energy() <= best_energy:
                best_energy = cc.get_energy()
                ex1 = i

        best_energy = 0.
        for i in list2[:n_it]:
            c = deepcopy(self)
            c.atoms[i].elem = self._element1
            cc = Cluster(c.atoms) #类别降为Cluster
            optimize(cc, fix_coord=[7]*i+[0]+[7]*(self.get_size()-1-i))
            if cc.get_energy() <= best_energy:
                best_energy = cc.get_energy()
                ex2 = i

        if ex1 < self._axis_num1:
            self.atoms[ex1].exchange_coordinate(self.atoms[ex2])
        else:
            a = self.atoms[ex1]
            self._del(a)
            a.elem = self._element2
            self._add(a)
            a = self.atoms[ex2]
            self._del(a)
            a.elem = self._element1
            self._add(a)

    def change_naxis(self):
        '''增减对称轴上的点
调用：_add_axis_point，_add_nonaxis_point'''
        if isinstance(self, SymCsBiCluster):
            increase_probability = 0.5 if max(self._axis_num1,self._axis_num2) >= self._order else 1.
        else:
            increase_probability = 5 ** (-(max(self._axis_num1,self._axis_num2)/self._order))
        for _ in range(10):
            copy_ = deepcopy(self)
            if random() <= increase_probability: #增对称轴上的点
                if copy_._nonaxis_num1 == 0:
                    atom_type = self._element2
                elif copy_._nonaxis_num2 == 0:
                    atom_type = self._element1
                else:
                    atom_type = self._element1 if random() < 0.5 else self._element2
                if atom_type == self._element1:
                    copy_._del(self._axis_num1 + randrange(copy_._nonaxis_num1))
                else:
                    copy_._del(self._atom_num1 + randrange(copy_._nonaxis_num2))
                for _ in range(copy_._order):
                    if not copy_._add_axis_point(atom_type):
                        break
                else:
                    self.deepcopy(copy_)
                    return
            else: #减少对称轴上的点
                if copy_._axis_num1 < copy_._order:
                    atom_type = self._element2
                elif copy_._axis_num2 < copy_._order:
                    atom_type = self._element1
                else:
                    atom_type = self._element1 if random() < 0.5 else self._element2
                if atom_type == self._element1:
                    for _ in range(copy_._order):
                        copy_._del(randrange(copy_._axis_num1))
                else:
                    for _ in range(copy_._order):
                        copy_._del(randrange(copy_.get_size()-copy_._axis_num2, copy_.get_size()))
                if copy_._add_nonaxis_point(atom_type):
                    self.deepcopy(copy_)
                    return


class SymCsBiCluster(SymBiCluster):
    '''处理Cs对称的团簇
    对称面为xy平面'''
    _order = 2 #对称群阶数
    sym = 'Cs'

    def _is_axis(self, a):
        '''判断a是否在对称面上'''
        return abs(a.z) < tolerance_min/2 * BOND(a.elem,a.elem)

    def _put_axis(self, a):
        '''将p修改为在对称面上'''
        a.z = 0

    def symmetry_atom(self, a):
        return [a, Atom(Point(a.x, a.y, -a.z), a.elem)]


class SymC2BiCluster(SymBiCluster):
    '''处理C2对称的团簇,绕z轴旋转180°对称'''
    _order = 2 #对称群阶数
    sym = 'C2'

    def _is_axis(self, a):
        '''判断a是否在对称轴上'''
        return hypot(a.x, a.y) < tolerance_min/2 * BOND(a.elem,a.elem)

    def symmetry_atom(self, a):
        return [a, Atom(Point(-a.x, -a.y, a.z), a.elem)]


class SymC3BiCluster(SymBiCluster):
    '''处理C3对称的团簇,绕z轴旋转120°对称'''
    _order = 3 #对称群阶数
    sym = 'C3'

    def _is_axis(self, a):
        '''判断a是否在对称轴上'''
        return hypot(a.x, a.y) < tolerance_min / sqrt(3) * BOND(a.elem,a.elem)

    def symmetry_atom(self, a):
        return [a, Atom(a.rotate(2*pi/3, Point(0,0,1)),a.elem), \
                Atom(a.rotate(4*pi/3, Point(0,0,1)),a.elem)]


class Alternate(BiCluster):
    #import atom
    #atom.tolerance_max = 1.15 #会导致其他地方使用该变量就不是1.25了

    def random(self):
        BiCluster.random(self, True, True)

    def clean(self,bonds=None):
        super(Alternate,self).clean(bonds)
        for bond in reversed(bonds):
            if self.atoms[bond[0]].elem == self.atoms[bond[1]].elem:
                bonds.remove(bond)
        self.write_xsd('temp.xsd',bonds)
        super(Alternate,self).clean(bonds)

    def mating(self, father, mother):
        f = deepcopy(father)
        m = deepcopy(mother)
        f.center()
        m.center()

        for _ in range(1000):
            f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

            face = uniform(-self._bond_len11/2, self._bond_len11/2) #杂交面
            up1 = [a for a in f.atoms[:self._atom_num1] if a.z>=face]
            up2 = [a for a in f.atoms[self._atom_num1:self._cls_size] if a.z>=face]
            down1 = [a for a in m.atoms[:self._atom_num1] if a.z<face]
            down2 = [a for a in m.atoms[self._atom_num1:self._cls_size] if a.z<face]
            if len(up1)+len(up2) != 0 and len(down1)+len(down2) != 0 and \
                len(up1)+len(down1) == self._atom_num1 and \
                len(up2)+len(down2) == self._atom_num2:
                for a in down1+down2:
                    a.z -= 5
                self.atoms = up1+down1+up2+down2
                bonds = self.get_bonds()
                for i,a in enumerate(up1):
                    for a2 in f.atoms:
                        if a2.z < face and a.bond(a2):
                            #与下半球有连接，该原子的键被割断
                            idx = np.argmin(np.array([a.dist(a3) for a3 in down2]))
                            bonds.append((i, self._atom_num1+len(up2)+idx))
                            break
                for i,a in enumerate(up2):
                    for a2 in f.atoms:
                        if a2.z < face and a.bond(a2):
                            #与下半球有连接，该原子的键被割断
                            idx = np.argmin(np.array([a.dist(a3) for a3 in down1]))
                            bonds.append((self._atom_num1+i, len(up1)+idx))
                            break
                self.clean(bonds)
                break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    @staticmethod
    def del_bond(bond, bonds):
        for i in range(len(bonds)):
            if bonds[i][0] == bond[0] and bonds[i][1] == bond[1]:
                del bonds[i]
                return
            elif bonds[i][0] == bond[1] and bonds[i][1] == bond[0]:
                del bonds[i]
                return

    def even_polyhefron_transform1(self):
        while True:
            bonds = self.get_bonds()
            new_bonds = deepcopy(bonds)
            conn = self.get_connections()
            a = randrange(self.get_size())
            b = choice(conn[a])
            conn[a].remove(b)
            e = choice(conn[a])
            ring0 = self.get_ring(b,a,e,conn)
            if len(ring0) != 4:
                continue
            f = ring0[-1]
            conn[a].remove(e)
            if not conn[a]:
                continue
            c = choice(conn[a])
            ring1 = self.get_ring(c,a,e,conn)
            if len(ring1) == 4:
                continue
            g = ring1[-1]
            conn[b].remove(a)
            conn[b].remove(f)
            if not conn[b]:
                continue
            d = choice(conn[b])
            ring2 = self.get_ring(d,b,f,conn)
            if len(ring2) == 4:
                continue
            h = ring1[-1]
            ring3 = self.get_ring(c,a,b,conn)
            if len(ring3) == 4:
                continue
            new_bonds.extend([[e,g],[f,h],[c,d]])
            self.del_bond([e,f], new_bonds)
            self.del_bond([c,g], new_bonds)
            self.del_bond([d,h], new_bonds)
            self.clean(new_bonds)
            break

    def even_polyhefron_transform2(self):
        while True:
            bonds = self.get_bonds()
            new_bonds = deepcopy(bonds)
            conn = self.get_connections()
            a = randrange(self.get_size())
            b,c,d = conn[a][:3]
            ring1 = self.get_ring(d,a,c,conn)
            if len(ring1) == 4:
                continue
            g = ring1[3]
            h = ring1[-1]
            ring2 = self.get_ring(b,a,c,conn)
            if len(ring2) == 4:
                continue
            e = ring2[-1]
            ring3 = self.get_ring(b,a,d,conn)
            if len(ring3) == 4:
                continue
            f = ring3[-1]
            self.del_bond([c,g], new_bonds)
            self.del_bond([d,h], new_bonds)
            self.del_bond([b,e], new_bonds)
            self.del_bond([b,f], new_bonds)
            new_bonds.extend([[b,g], [b,h], [c,e], [d,f]])
            self.clean(new_bonds)
            break

    def even_polyhefron_transform3(self):
        while True:
            bonds = self.get_bonds()
            new_bonds = deepcopy(bonds)
            conn = self.get_connections()
            a = randrange(self.get_size())
            b = choice(conn[a])
            conn[b].remove(a)
            c = choice(conn[b])
            ring = self.get_ring(a,b,c,conn)
            if len(ring) != 4:
                continue
            d = ring[-1]
            conn[a].remove(b)
            conn[a].remove(d)
            if not conn[a]:
                continue
            e = choice(conn[a])
            ring = self.get_ring(b,a,e,conn)
            x = ring[3]
            ring = self.get_ring(d,a,e,conn)
            if len(ring) == 4:
                continue
            y = ring[3]
            new_bonds.extend([[b,x], [d,y]])
            self.del_bond([a,d], new_bonds)
            self.del_bond([x,e], new_bonds)
            self.del_bond([y,e], new_bonds)
            conn[c].remove(b)
            conn[c].remove(d)
            if not conn[c]:
                continue
            f = choice(conn[c])
            ring = self.get_ring(b,c,f,conn)
            x = ring[-1]
            y = ring[3]
            self.del_bond([b,x], new_bonds)
            self.del_bond([f,y], new_bonds)
            new_bonds.extend([[e,x], [e,y], [a,f]])
            self.atoms[e].set_coordinate(self.atoms[x]+self.atoms[y])
            self.atoms[a].set_coordinate(self.atoms[x]+self.atoms[f])
            self.clean(new_bonds)
            break


###############################################################################
#                       多元团簇
###############################################################################
class MultiCluster(FixCluster):
    '''多元团簇类'''

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(MultiCluster,cls).config()
        elements = cfg.get('cluster', 'element').split()
        if len(elements)  <= 2:
            print('Please input more than two elements in config.ini.')
            os._exit(1)
        cls.set_element(elements)
        atom_numbers = cfg.get('cluster', 'atom_num').split()
        if len(atom_numbers) <= 2:
            print('Please input more than two atom numbers in config.ini.')
            os._exit(1)
        if len(atom_numbers) != len(elements):
            print('numbers of element and atom_num does not match.')
            os._exit(1)
        cls.set_atom_num([int(n) for n in atom_numbers])

    @classmethod
    def set_element(cls, elements):
        for e in elements:
            if get_element_id(e) == -1:
                print('element error in config.ini.')
                os._exit(1)
        cls._elem_types = deepcopy(elements)

    @classmethod
    def set_atom_num(cls, atom_numbers):
        for n in atom_numbers:
            if n <= 0:
                print('atom number error in config.ini.')
                os._exit(1)
        cls._atom_nums = deepcopy(atom_numbers)
        cls._cls_size = sum(atom_numbers)

    def check(self):
        '''检查结构是否合法
        '''
        elements = self.get_elements_count()
        for e,n in zip(self._elem_types, self._atom_nums):
            if elements[e] != n:
                return False
        return True

    def random(self):
        '''产生随机结构'''
        a,b,c = abs(random())+0.1, abs(random())+0.1, abs(random()+0.1)
        d = (self._cls_size/a/b/c)**(1./3)
        elements = reduce(list.__add__, [[e]*n for n,e in zip(self._atom_nums,self._elem_types)])

        while self.get_size() < self._cls_size:
            for _ in range(20):
                atom = Atom(Point(random()*a, random()*b, random()*c)*d, choice(elements))
                if self.is_legal(atom):
                    self.add_atom(atom)
                    elements.remove(atom.elem)
                    break
            else:
                d += 0.1
                self.atoms = []
                elements = reduce(list.__add__, [[e]*n for n,e in zip(self._atom_nums,self._elem_types)])

    def self_mating(self):
        self.center()
        self.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

        while True:
            locus = uniform(-1,1)
            up = Cluster([a for a in self.atoms if a.z>=locus])
            down = Cluster([a for a in self.atoms if a.z<locus])
            if up.get_size() > 0 and down.get_size() > 0:
                up.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
                down.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
                break

        self._merge(up.atoms, down.atoms)

    def exchange(self):
        '''交换两元素'''
        self.center()
        self.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

        for _ in range(1, self._cls_size//5+2):
            a1 = None
            for __ in range(100):
                if a1 is None:
                    a1 = choice(self.atoms)
                    if a1.z < 0:
                        a1 = None
                else:
                    a2 = choice(self.atoms)
                    if a2.z > 0 and a1.elem != a2.elem:
                        break
            a1.elem,a2.elem = a2.elem,a1.elem

class SymMultiCluster(MultiCluster):
    '''未测试'''
    _symmetry_list = ['Ci','Cs','C2','C3','C5']
    _last_symmetry = randrange(5)

    def set_symmetry(self, sym):
        if sym == 'Ci':
            self._order = 2
        elif sym == 'Cs':
            self._order = 2
        elif sym == 'C2':
            self._order = 2
        elif sym == 'C3':
            self._order = 3
        elif sym == 'C5':
            self._order = 5
        else:
            print('set_symmetry parameter error.')
            raise
        self.sym = sym

    def __repr__(self):
        '''转字符串'''
        s = str(self._energy) + '\t' + str(self._inertia) + '\t' + str(self._inertia2) + '\t' + self.sym+'\n'
        for a in self.atoms:
            s += str(a) + '\n'
        return s

    def _arrange(self):
        '''每次从别的地方(例如读recover)读取结构,都要调用该函数调整原则顺序，并给_nonaxis_num和axis_num赋值
        dmol优化成高对称性后会改变原子排列顺序！'''
        self._axis_num = self._nonaxis_num = 0
        self._set_symmetry_axis()

        i = 0
        while i < self._cls_size - self._axis_num:
            if self._is_axis(self.atoms[i]):
                self._axis_num += 1
                self.add_atom(self.atoms[i])
                del self.atoms[i]
                continue
            self._nonaxis_num += 1
            sym_atoms = self.symmetry_atom(self.atoms[i])
            for p in sym_atoms[1:]:
                for j in range(i+1, self._cls_size-self._axis_num):
                    if (self.atoms[j]-p).norm() < 0.1:
                        self.atoms[i+1],self.atoms[j] = self.atoms[j],self.atoms[i+1] #这里可能有错，参看SymBiCluster
                        break
                else:
                    print('cluster loose symmetry.')
                    raise utils.SymmetryErr('cluster loose symmetry.')
            i += self._order

    def read_cga_obj(self, file_obj):
        '''从文件中读取结构
        参数file: 文件对象
        废弃，该函数不该出现在该文件'''
        raise
        try:
            line = file_obj.readline() #读取第一行,能量,转动惯量
            content = line.split()
            if len(content) != 4:
                raise utils.FileFormatErr(line)
            try:
                self.set_energy(float(content[0]))
                self.set_inertia(float(content[1]), float(content[2]))
                self.set_symmetry(content[3]) #这里与cluster不同
            except:
                raise utils.FileFormatErr(line)
            self.atoms = []
            while True:
                line = file_obj.readline()
                content = line.split()
                if len(content) != 4 or len(content[0]) > 2:
                    break
                try:
                    self.add_atom(Point(float(content[1]), float(content[2]), float(content[3])), content[0])
                except:
                    raise utils.FileFormatErr(line)
        except:
            print('file format error! error at line: ', line)
            raise utils.FileFormatErr(line)
        self._arrange()
        return line

    def _is_axis(self, a):
        '''判断p是否在对称轴（面）上'''
        if self.sym == 'Ci':
            return hypot(a.x, a.y) < tolerance_min/2 * BOND(a.elem,a.elem)
        if self.sym == 'Cs':
            return abs(a.z) < tolerance_min/2 * BOND(a.elem,a.elem)
        elif self.sym == 'C2':
            return hypot(a.x, a.y) < tolerance_min/2 * BOND(a.elem,a.elem)
        elif self.sym == 'C3':
            return hypot(a.x, a.y) < tolerance_min / sqrt(3) * BOND(a.elem,a.elem)
        elif self.sym == 'C5':
            return hypot(a.x, a.y) < tolerance_min / (2*sin(0.2*pi)) * BOND(a.elem,a.elem)

    def _put_axis(self, a):
        '''将p修改为在对称面上'''
        if self.sym == 'Ci':
            a.set(Point(0.))
        if self.sym == 'Cs':
            a.z = 0.
        else:
            a.x = a.y = 0.

    def symmetry_atom(self, a):
        '''计算a的对称原子（包括a）'''
        if self.sym == 'Ci':
            return [a, -a]
        if self.sym == 'Cs':
            return [a, Atom(a.x, a.y, -a.z, a.elem)]
        elif self.sym == 'C2':
            return [a, Atom(-a.x, -a.y, a.z, a.elem)]
        elif self.sym == 'C3':
            return [a, Atom(a.rotate(2*pi/3, Point(0,0,1)),a.elem), \
                Atom(a.rotate(4*pi/3, Point(0,0,1)),a.elem)]
        elif self.sym == 'C5':
            return [a, Atom(a.rotate(2*pi/5, Point(0,0,1)), a.elem), \
                Atom(a.rotate(4*pi/5, Point(0,0,1)), a.elem), \
                Atom(a.rotate(6*pi/5, Point(0,0,1)), a.elem), \
                Atom(a.rotate(8*pi/5, Point(0,0,1)), a.elem)]

    def _add(self,a):
        '''增加一个原子
        非对称轴原子加在前端，对称轴原子加在后端'''
        if self._is_axis(a):
            self.add_atom(a)
        else:
            self.atoms = self.symmetry_atom(a) + self.atoms

    def _del(self, a):
        '''删除一个原子
        参数a：Point类型或数组下标
        被调用:perturbation,_increase_naxis'''
        if isinstance(a, Atom): #p为点
            index = self.atoms.index(a)
        else: #a为点的下标
            index = a
            a = self.atoms[index]
        if self._is_axis(a):
            del self.atoms[index]
            self._axis_num -= 1
        else:
            index = index // self._order * self._order
            del self.atoms[index : index+self._order]
            self._nonaxis_num -= 1

    def _add_axis_point(self, elem):
        '''增加一个对称轴上的点
        返回True或False'''
        new_atom = Atom(Point(0), elem)
        for _ in range(100):
            neighbour = choice(self.atoms)
            new_atom = neighbour + Point.random_plane(BOND(neighbour.elem, elem))
            self._put_axis(new_atom)
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1:
                self._add(new_atom)
                return True
        else:
            return False

    def _add_nonaxis_point(self, elem):
        '''增加一个非对称轴上的点
        返回True或False'''
        new_atom = Atom(Point(0), elem)
        for _ in range(100):
            neighbour = choice(self.atoms)
            new_atom = neighbour + Point.random(BOND(neighbour.elem, elem))
            if self._is_axis(new_atom):
                continue
            if self.is_legal(new_atom) and self.get_degree(new_atom) > 1:
                self._add(new_atom)
                return True
        else:
            return False

    def random(self):
        '''产生随机结构'''
        SymMultiCluster._last_symmetry += 1
        self.set_symmetry(SymMultiCluster._symmetry_list[SymMultiCluster._last_symmetry%len(SymMultiCluster._symmetry_list)])
        if self.sym == 'Ci' and len([x for x in self._atom_nums if x%2]) > 1:
            SymMultiCluster._last_symmetry += 1
            self.set_symmetry(SymMultiCluster._symmetry_list[SymMultiCluster._last_symmetry%len(SymMultiCluster._symmetry_list)])
        h = uniform(1, 1+self._cls_size**(1./3))
        r = sqrt(self._cls_size/h)/4
        self.atoms = []
        remain_elem = reduce(list.__add__, [[e]*n for e,n in zip(self._elem_types, self._atom_nums)])
        while self.get_size() < self._cls_size:
            e = choice(remain_elem)
            for _ in range(100): #每个原子尝试100次
                now = Atom(Point.random_plane(random()*r) + Point(0,0,random()*h), e)
                if self._is_axis(now) or remain_elem.count(e) < self._order: #离对称轴太近,就放到对称轴上
                    self._put_axis(now)
                if self.is_legal(now):
                    self._add(now)
                    if self._is_axis(now):
                        remain_elem.remove(e)
                    else:
                        for __ in range(self._order):
                            remain_elem.remove(e)
                    break
            else:
                self.atoms = []
                remain_elem = reduce(list.__add__, [[el]*n for el,n in zip(self._elem_types, self._atom_nums)])
                h += 0.1
                r += 0.1

    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for _ in range(100):
            c = deepcopy(self)
            atom = choice(c.atoms)
            c._del(atom)
            if c._is_axis(atom):
                if c._add_axis_point(atom.elem) and not c.atoms[-1].bond(atom):
                    self.deepcopy(c)
                    return [self._cls_size-1]
            else:
                if c._add_nonaxis_point(atom.elem):
                    for i in range(c._order):
                        if c.atoms[i].bond(atom):
                            break
                    else:
                        self.deepcopy(c)
                        return list(range(self._order))
        else:
            utils.write_log('mutation_one failure.')
            self.perturbation()
            return []

    def change_naxis(self):
        '''增减对称轴上的点'''
        if self.sym == 'Cs':
            increase_probability = 0.5 if self._axis_num >= self._order else 1.
        else:
            increase_probability = 5 ** (-(self._axis_num/self._order))
        if random() <= increase_probability: #增对称轴上的点
            del_idx = randrange(self._nonaxis_num)
            elem = self.atom[del_idx].elem
            self._del(del_idx)
            for _ in range(self._order):
                if not self._add_axis_point(elem):
                    break
            else:
                return
        else: #减少对称轴上的点
            for _ in range(self._order):
                self._del(randrange(self.get_size()-self._axis_num, self.get_size()))
            if self._add_nonaxis_point():
                return
        utils.write_log('change_naxis failure.')
        self.perturbation()


################################################################################
#               分子团簇
###############################################################################
class MoleCluster(FixCluster):
    '''分子团簇类
    由基本分子单元组成.坐标和元素存储规则为依次存储每个分子
    使用前需设置基本分子和分子数'''
    _mole_bond = 2.77 #分子键长

    @classmethod
    def set_unit(cls, mole):
        '''设置组成分子
        使用该类前必须调用该函数'''
        cls._unit = deepcopy(mole) #组成分子
        cls._mole_size = mole.get_size() #组成分子的原子数

    @classmethod
    def set_mole_num(cls, mole_num):
        '''设置分子数
        使用该类前必须调用该函数'''
        cls._cls_size = mole_num

    @classmethod
    def get_mole_num(cls, self):
        return cls._cls_size

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(MoleCluster,cls).config()
        cls.set_mole_num(cfg.getint('cluster', 'atom_num'))
        mole_file_name = cfg.get('cluster', 'element')
        assert os.path.isfile(mole_file_name)
        principle = read(mole_file_name)
        cls.set_unit(principle)

    def check(self):
        '''检查结构是否正确
        '''
        if self.get_size() != self._cls_size*self._unit.get_size():
            return False
        for i in range(self._cls_size):
            for j in range(self._unit.get_size()):
                if self.atoms[i*self._unit.get_size()+j].elem != self._unit.element[j]:
                    return False
        return True

    def _get_mole(self, i):
        '''获取第i个分子
        被调用: mating, perturbation'''
        return Cluster(self.atoms[i*self._mole_size : (i+1)*self._mole_size])

    def _del(self, i):
        '''删除第i个分子
        被调用: perturbation'''
        if i < 0:
            i += self._cls_size
        if i >= self._cls_size:
            print('self.del_ parameter error.')
            raise
        del self.atoms[i*self._mole_size : (i+1)*self._mole_size]

    def _arrange(self):
        '''可能由于优化的原因，一个分子的原子没有在一起，该函数将使它们重新排在一起'''
        for i in range(1,self._unit.get_size()):
            d = self._unit.atoms[0].dist(self._unit.atoms[i])
            for j in range(self.get_mole_num()):
                min_d = float("inf")
                for k in range(j*self._mole_size+1, self.get_size()):
                    if self.atoms[k].elem != self.atoms[j+i].elem:
                        continue
                    d_now = self.atoms[j*self._mole_size].dist(self.atoms[k])
                    if abs(d_now - d) < min_d:
                        min_d = abs(d_now -d)
                        min_pos = k
                self.atoms[j*self._mole_size+i],self.atoms[min_pos] = self.atoms[min_pos],self.atoms[j*self._mole_size+i]

    def is_legal(self, mole):
        '''判断团簇之外的某个分子是否合法
        被调用: random, perturbation'''
        for i in range(mole.get_size()):
            for j in range(self.get_size()):
                if mole.atoms[i].dist(self.ato[j]) < 1.5:
                    return False
        return True

    @classmethod
    def _change_mole(cls, mole, move, rot_x, rot_y, rot_z):
        '''改变一个分子(平移+旋转)
        rot_x, rot_y, rot_z单位为弧度
        被调用: random, perturbation'''
        O = mole.get_center()
        mole.move(O)
        mole.rotate(rot_x, rot_y, rot_z)
        mole.move(move-O)

    def random(self):
        '''产生随机结构, 在球形区域内'''
        for diameter in np.linspace(0.4*self.get_mole_num()**(1./3), 2*self.get_mole_num()**(1./3), 100):
            self.atoms = []
            for _ in range(self._cls_size):
                for __ in range(100): #尝试100次
                    mole = deepcopy(self._unit)
                    self._change_mole(mole,
                                        Point.random(diameter),
                                        2*pi*random(), 2*pi*random(), 2*pi*random())
                    if self.is_legal(mole): #i的位置合法,退出it循环
                        self.atoms.extend(mole.atoms)
                        break
                else: #没有找到i合适位置
                    break
            else: #找到合适结构
                return True
        else: #没有找到合适结构
            return False

    def mating(self, father, mother):
        '''交叉
        沿z轴.取father的下半球,mother的上半球'''
        f = deepcopy(father)
        f.center()
        m = deepcopy(mother)
        m.center()

        f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
        m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
        for _ in range(100):
            cut_plane = uniform(0, max([p.z for p in f.atoms])/2)
            mole_father = [f._get_mole(i) for i in range(self._cls_size)]
            mole_mother = [m._get_mole(i) for i in range(self._cls_size)]
            semi_father = [mole for mole in mole_father if mole.get_center().z<=cut_plane]
            semi_mother = [mole for mole in mole_mother if mole.get_center().z>cut_plane]
            count = len(semi_father) + len(semi_mother)
            if count == self._cls_size:
                break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            return

        self._merge(reduce(list.__add__, [mole.atoms for mole in semi_father]),
                    reduce(list.__add__, [mole.atoms for mole in semi_mother]),
                    self._mole_bond)

    def perturbation(self, step=FixCluster.amplitude):
        '''给结构微小扰动
        平均每个分子扰动100次,每次扰动不超过0.02'''
        angle = 40.*step/0.02*pi/180 #角度制转弧度制
        for _ in range(self._cls_size*100):
            index = randrange(self._cls_size)
            old_mole = self._get_mole(index)
            self._del(index)
            new_mole = deepcopy(old_mole)
            self._change_mole(new_mole,
                               Point.random(step),
                               angle*random(), angle*random(), angle*random())
            #判断新分子位置是否合法
            if not self.is_legal(new_mole):
                self.atoms.extend(old_mole.atoms)
                continue
            self.atoms.extend(new_mole.atoms)

    def adjust_H(self):
        '''调节氢键朝向，供水团簇使用'''
        def rotate_H():
            O = deepcopy(mole.atoms[0]) #以O原子为中心旋转
            mole.move(O)
            mole.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            mole.move(-O)

        n = self._cls_size if type(self)==MoleCluster else self._cls_size+1
        for i in range(n):
            mole = self._get_mole(i)
            rotate_H()
        relax_index = [i for i in range(self.get_size()) if self.element[i]=='H']
        fix_coord = [7]*self.get_size()
        for i in relax_index:
            fix_coord[i] = 0
        optimize(self, fix_coord=fix_coord)


class Hydrate(MoleCluster):
    '''水合物类
    由基本分子单元组成.坐标和元素存储规则为依次存储每个分子
    使用前需设置基本分子和分子数
    客体分子存放在最后面'''

    @classmethod
    def set_unit(cls, mole, object_):
        '''设置组成分子
        使用该类前必须调调用该函数'''
        cls.set_unit(mole)
        cls._object = object_ #客体分子

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(Hydrate,cls).config()
        mole_names = cfg.get('cluster', 'element').split()
        if len(mole_names) != 2:
            print('Please input two elements in config.ini.')
            raise
        assert os.path.isfile(mole_names[0])
        assert os.path.isfile(mole_names[1])

        nums = cfg.get('cluster', 'atom_num').split()
        if len(nums) != 2:
            print('Please input two atom numbers in config.ini.')
        if int(nums[0]) == 1:
            cls.set_mole_num(int(nums[1]))
            mole_names[0],mole_names[1] = mole_names[1],mole_names[0]
        elif int(nums[1]) == 1:
            cls.set_mole_num(int(nums[0]))
        else:
            print('only support one object molecular.')
            raise

        principle = read(mole_names[0])
        object_ = read(mole_names[1])
        cls.set_unit(principle, object_)

    def check(self):
        '''检查结构是否正确
        '''
        if self.get_size() != self._cls_size*self._unit.get_size() + Hydrate._object.get_size():
            return False
        for i in range(self._cls_size):
            for j in range(self._unit.get_size()):
                if self.atoms[i*self._unit.get_size()+j].elem != self._unit.atoms[j].elem:
                    return False
        for j in range(1, Hydrate._object.get_size()+1):
            if self.atoms[-j].elem != Hydrate._object.atoms[-j].elem:
                return False
        return True

    def _get_mole(self, index):
        '''获取第i个分子
        被调用: mating, perturbation'''
        if index == self._cls_size:
            mole = Cluster()
            mole.atoms = self.atoms[self._cls_size*self._mole_size : ]
            return mole
        else:
            return self._get_mole(self, index)

    def _del(self, index):
        '''删除第i个分子
        被调用: perturbation'''
        if index < 0:
            index += self._cls_size
        if index == self._cls_size:
            del self.atoms[self._cls_size*self._mole_size : ]
        else:
            self._del(self, index)

    def random(self):
        '''产生随机结构, 在球形区域内'''
        for diameter in np.linspace(0.4*self.get_size()**(1./3), 2*self.get_size()**(1./3), 100):
            #产生客体分子
            mole = deepcopy(Hydrate._object)
            self._change_mole(mole,
                                Point.random(diameter),
                                2*pi*random(), 2*pi*random(), 2*pi*random())
            self.atoms = mole.atoms
            #产生主体分子
            for _ in range(self._cls_size):
                for __ in range(100): #尝试100次
                    mole = deepcopy(self._unit)
                    self._change_mole(mole, Point.random(diameter),
                                      2*pi*random(), 2*pi*random(), 2*pi*random())
                    if self.is_legal(mole): #i的位置合法,退出it循环
                        self.atoms= mole.atoms + self.atoms
                        break
                else: #没有找到i合适位置
                    break
            else: #找到合适结构
                return True
        else: #没有找到合适结构
            return False

    def mating(self, father, mother):
        '''交叉
        沿z轴.取father的下半球,mother的上半球'''
        def check():
            for m in mole_father:
                for p in m.atoms:
                    for n in mole_mother:
                        for q in n.atoms:
                            if p.dist(q) < 1.5:
                                return False
            m = self._get_mole(self._cls_size)
            if father._get_mole(self._cls_size).center() > 0: #客体分子来自父代，则跟母代做检测
                for p in m.atoms:
                    for n in mole_mother:
                        for q in n.atoms:
                            if p.dist(q) < 1.5:
                                return False
            else: #客体分子来自母代，则跟父代做检测
                for p in m.atoms:
                    for n in mole_father:
                        for q in n.atoms:
                            if p.dist(q) < 1.5:
                                return False
            return True

        f = deepcopy(father)
        f.center()
        m = deepcopy(mother)
        m.center()
        for _ in range(100):
            f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            mole_father = [f._get_mole(i) for i in range(self._cls_size)] #不含客体分子
            mole_mother = [m._get_mole(i) for i in range(self._cls_size)]
            semi_father = [mole for mole in mole_father if mole.get_center().z>0]
            semi_mother = [mole for mole in mole_mother if mole.get_center().z<=0]
            if len(semi_father) == 0 or len(semi_mother) == 0:
                continue
            if len(semi_father) + len(semi_mother) == self._cls_size:
                break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            return

        semi_father = reduce(list.__add__, [mole.atoms for mole in semi_father])
        semi_mother = reduce(list.__add__, [mole.atoms for mole in semi_mother])
        if f._get_mole(self._cls_size).center() > 0:
            semi_father.extend(f._get_mole(self._cls_size).atoms)
        else:
            semi_mother.extend(m._get_mole(self._cls_size).atoms)

        self._merge(semi_father, semi_mother, self._mole_bond)

    def perturbation(self, step=FixCluster.amplitude):
        '''给结构微小扰动
        平均每个分子扰动100次,每次扰动不超过0.02
        ga中产生子代根据参数个数决定是单双亲操作,故这里改成1个参数step,去掉了angle'''
        angle = 40.*step/0.02*pi/180 #角度制转弧度制
        for _ in range(self._cls_size*100):
            index = randrange(self._cls_size)
            old_mole = self._get_mole(index)
            self._del(index)
            new_mole = deepcopy(old_mole)
            self._change_mole(new_mole,
                               Point.random(step),
                               angle*random(), angle*random(), angle*random())
            #判断新分子位置是否合法
            if not self.is_legal(new_mole):
                if index < self._cls_size: #主体分子
                    self.atoms = old_mole.atoms + self.atoms
                else: #客体分子
                    self.atoms.extend(old_mole.atoms)
                continue
            #合法
            if index < self._cls_size: #主体分子
                self.atoms = new_mole.atoms + self.atoms
            else: #客体分子
                self.atoms.extend(new_mole.atoms)


class BiMoleCluster(BiCluster, FixCluster):
    '''ToDo'''

    def mating(self, father, mother):
        '''杂交,沿z=face平面杂交'''
        f = deepcopy(father)
        m = deepcopy(mother)
        f.center()
        m.center()

        for _ in range(1000):
            f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

            face = uniform(-self._bond_len11/2, self._bond_len11/2) #杂交面
            up1 = [a for a in f.atoms[:self._atom_num1] if a.z>=face]
            up2 = [a for a in f.atoms[self._atom_num1:self._cls_size] if a.z>=face]
            down1 = [a for a in m.atoms[:self._atom_num1] if a.z<face]
            down2 = [a for a in m.atoms[self._atom_num1:self._cls_size] if a.z<face]
            if len(up1)+len(up2) == 0 or len(down1)+len(down2) == 0:
                continue
            if len(up1)+len(down1) == self._atom_num1 and \
                len(up2)+len(down2) == self._atom_num2:
                self._merge(down1, down2, up1, up2)
                break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            utils.write_log('mating failure')


class LigandCluster(BiCluster):

    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = super(BiCluster,cls).config()
        e = cfg.get('cluster', 'element').split()
        cls.set_element(e[0], e[1])
        atom_num = cfg.get('cluster', 'atom_num').split()
        if len(atom_num) != 2:
            print('Please input two atom num in config.ini.')
            os._exit(1)
        cls.set_atom_num(int(atom_num[0]), int(atom_num[1]))

    @classmethod
    def set_atom_num(cls, atom_num1, atom_num2):
        super(LigandCluster, cls).set_atom_num(atom_num1, atom_num2)
        cls._mole_size = cls.ligand.get_size()
        cls._cls_size = cls._atom_num1 + cls._atom_num2*cls._mole_size

    @classmethod
    def set_element(cls, e1, e2):
        cls.subject = e1
        assert os.path.isfile(e2)
        cls.ligand = read(e2)
        cls._bond_len11 = BOND(cls.subject, cls.subject)
        cls._bond_len12 = BOND(cls.subject, cls.ligand.atoms[0].elem)
        cls._bond_len22 = BOND(cls.ligand.atoms[0].elem, cls.ligand.atoms[0].elem)

    def check(self):
        return True
        
    def from_bonds(self, bonds):
        elem = ['Au']*self._atom_num1 + ['S','H']*self._atom_num2
        for i,e in enumerate(elem):
            print(i,e)
        print(bonds)
        self.atoms = [Atom(0,0,0, elem[0])]
        layer = 0
        previous_atoms = {0} #已有原子
        last_atoms = {0} #上一层原子
        idx = 1 #已添加的原子序号
        while True:
            current_atoms = set() #当前层原子
            for i in last_atoms:
                for j,k in bonds:
                    if j==i and k not in previous_atoms:
                        current_atoms.add(k)
                    if k==i and j not in previous_atoms:
                        current_atoms.add(j)
            if not current_atoms:
                break
            layer += 1
            n = len(current_atoms)
            r = 1.5/sin(pi/n)
            for i in range(len(current_atoms)):
                self.add_atom(Atom(r*cos(2*pi*i/n), r*sin(2*pi*i/n), layer*3, \
                    elem[idx]))
                idx += 1
            previous_atoms.update(current_atoms)
            last_atoms = current_atoms
        from .io import write_xsd
        write_xsd(self, 'pre.xsd', bonds)
        self.clean()
        write_xsd(self, 'after.xsd')

    def random(self):
        def random_Au_core():
            self.add_atom(Point(0), 'Au')
            self.add_atom(Point(1), 'Au')
            for i in range(2, self._atom_num1-self._atom_num2+len(staple_num)):
                for _ in range(1000): #尝试1000次
                    d = uniform(tolerance_min, tolerance_max) * self._bond_len11
                    a = self.atoms[i-1] + Point.random(d)
                    if self.get_degree(a) > 1 and self.is_legal(a):
                        self.add_atom(a)
                        break
                else:
                    print('generate initial structure failure.')
                    raise
                    
        staple_num = [] #每个订书针的配体个数
        n = self._atom_num2
        while n:
            staple_num.append(randint(1, min(5,n)))
            n -= staple_num[-1]
        core_num = self._atom_num1 - (self._atom_num2-len(staple_num)) #Au核原子数
        
        random_Au_core()
        self.center()
        bonds = []
        for i in range(self._atom_num1-self._atom_num2+len(staple_num)):
            for j in range(i):
                if self.atoms[i].bond(self.atoms[j]):
                    bonds.append((i,j))
        
        #安装配体
        Au_idx = core_num
        staple_idx = self._atom_num1
        atoms_SH = []
        for staple in staple_num:
            #产生订书针
            for i in range(staple-1):
                bonds.append((Au_idx+i, staple_idx+i*2)) #Au跟前面的S连
                bonds.append((staple_idx+i*2, staple_idx+i*2+1)) #前面S上的H
                bonds.append((Au_idx+i, staple_idx+(i+1)*2)) #下一个S跟刚才的Au连
            while True: #针脚连接的2个Au原子不可离太远
                v1,v2 = sample(list(range(core_num)), 2)
                if self.atoms[v1].dist(self.atoms[v2]) < 3*staple:
                    break
            bonds.append((v1, staple_idx))
            staple_idx += staple*2
            bonds.append((v2, staple_idx-2))
            Au_idx += staple-1
            bonds.append((staple_idx-2, staple_idx-1)) #最后一个S上的H
            #产生订书针坐标
            v1 = self.atoms[v1]*uniform(5,10)
            v2 = self.atoms[v2]*uniform(5,10)
            if staple != 1:
                r = (v2-v1) / (2*(staple-1))
                for i in range(staple-1):
                    atoms_SH.append(Atom(v1+2*i*r, 'S'))
                    atoms_SH.append(Atom(self.atoms[-1] * (1+3/self.atoms[-1].norm()), 'H'))
                    self.add_atom(v1+(2*i+1)*r, 'Au')
            atoms_SH.append(Atom(v2, 'S'))
            atoms_SH.append(Atom(self.atoms[-1] * (1+3/self.atoms[-1].norm()), 'H'))
        self.atoms.extend(atoms_SH)
                
        from .io import write_xsd
        write_xsd(self, 'pre.xsd', bonds)
        self.clean(bonds)
        self.clean(bonds)
        write_xsd(self, 'after.xsd', bonds)

    def random_walk(self):
        bonds = self.get_bonds()
        for it in range(100):
            idx1 = choice(list(range(self._atom_num1, self._cls_size, self._mole_size))) #随机选择一个S
            con_Au = [i for i in bonds[idx1] if self.atoms[i].elem==self.subject and len(bonds[i])>2] #与该S相邻的Au列表
            if not con_Au:
                continue
            idx2 = choice(con_Au) #随机选择一个相邻的Au
            idx3 = choice(list(range(self._atom_num1))) #选择新的Au连接
            if self.atoms[idx3].dist(self.atoms[idx1]) > 6 or len(bonds[idx3]) <= 2 or idx3==idx2:
                continue
            bonds[idx1].remove(idx2) #断掉该S-Au键
            bonds[idx2].remove(idx1)
            bonds[idx1].append(idx3) #连接新S-Au键
            bonds[idx3].append(idx1)
            self.clean(bonds)
            break

    def mating(self, father, mother):
        '''杂交,沿z=face平面杂交'''
        f = deepcopy(father)
        m = deepcopy(mother)
        f.center()
        m.center()

        for it in range(1000):
            f.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())
            m.rotate(2*pi*random(), 2*pi*random(), 2*pi*random())

            face = uniform(-self._bond_len11/2, self._bond_len11/2) #杂交面
            up_Au = [i for i in range(self._atom_num1) if f.atoms[i].z>=face]
            down_Au = [i for i in range(self._atom_num1) if m.atoms[i].z<face]
            if len(up_Au) + len(down_Au) != self._atom_num1 or len(up_Au)==0 or len(down_Au)==0:
                continue
            up_S = [i for i in range(self._atom_num1, self._cls_size, 2) if f.atoms[i].z>=face]
            down_S = [i for i in range(self._atom_num1, self._cls_size, 2) if m.atoms[i].z<face]
            if len(up_S) + len(down_S) != self._atom_num2 or len(up_S)==0 or len(down_S)==0:
                continue
                
            bonds = []
            #Au-Au
            for i,idx1 in enumerate(up_Au):
                for j,idx2 in enumerate(up_Au[:i]):
                    if f.atoms[idx1].bond(f.atoms[idx2]):
                        bonds.append((i,j))
            for i,idx1 in enumerate(down_Au):
                for j,idx2 in enumerate(down_Au[:i]):
                    if m.atoms[idx1].bond(m.atoms[idx2]):
                        bonds.append((len(up_Au)+i, len(up_Au)+j))
            for i,idx1 in enumerate(up_Au):
                for j,idx2 in enumerate(down_Au):
                    if f.atoms[idx1].bond(m.atoms[idx2]):
                        bonds.append((i, len(up_Au)+j))
            #Au-S
            for i,idx1 in enumerate(up_Au):
                for j,idx2 in enumerate(up_S):
                    if f.atoms[idx1].bond(f.atoms[idx2]):
                        bonds.append((i, self._atom_num1+j*2))
            for i,idx1 in enumerate(up_Au):
                for j,idx2 in enumerate(down_S):
                    if f.atoms[idx1].bond(m.atoms[idx2]):
                        bonds.append((i, self._atom_num1+(len(up_S)+j)*2))
            for i,idx1 in enumerate(down_Au):
                for j,idx2 in enumerate(up_S):
                    if m.atoms[idx1].bond(f.atoms[idx2]):
                        bonds.append((self._atom_num1+(len(up_S)+i)*2, j))
            for i,idx1 in enumerate(down_Au):
                for j,idx2 in enumerate(down_S):
                    if m.atoms[idx1].bond(m.atoms[idx2]):
                        bonds.append((self._atom_num1+(len(up_S)+i)*2, \
                            self._atom_num1+(len(up_S)+j)*2))
            #S-H
            for i in range(self._atom_num1, self._cls_size, 2):
                bonds.append((i,i+1))
            self.clean(bonds)
            break
        else:
            print('mating failure')
            self.deepcopy(father)
            self.perturbation() #改用变异
            utils.write_log('mating failure')

    def adjust_H(self):
        bonds = self.get_bonds()
        for i in range(self._atom_num1, self._cls_size, self._mole_size):
            if len(bonds[i+1]) != 1 or self.atoms[bonds[i+1][0]].elem != 'S':
                for it in range(100):
                    self.atoms[i+1] = self.atoms[i]+Point.random(1.41)
                    self.atoms[i+1].elem = 'H'
                    for a in self.atoms[:self._atom_num1]:
                        if a.dist(self.atoms[i+1]) < 2.5: #H要远离所有Au
                            break
                    else:
                        break
                else:
                    print('adjust_H failure.')


class AuCl(BiCluster):

    def get_staple(self, idx, bonds):
        '''根据Cl的下标查找所在的订书钉
        返回元素下标列表'''
        def find_chain(chain):
            next = None
            for v1,v2 in bonds:
                if v1 == chain[-1] and v2 != chain[-2]:
                    if next is None:
                        next = v2
                    else:
                        return
                elif v2 == chain[-1] and v1 != chain[-2]:
                    if next is None:
                        next = v1
                    else:
                        return
                chain.append(next)
                
        staple = [idx]
        
        ends = []
        for v1,v2 in bonds:
            if v1 == idx:
                ends.append(v2)
            elif v2 == idx:
                ends.append(v1)
        if len(ends) != 2:
            return staple
            
        chain1 = [idx, ends[0]]
        find_chain(chain1)
        chain2 = [idx, ends[1]]
        find_chain(chain2)
        if len(chain1)%2 or len(chain2)%2 or \
            len(set(chain1+chain2)) < len(chain1)+len(chain2):
            print(chain1, chain2)
            raise
        return list(reversed(chain1))[:-1] + chain2
    
    def enlarge_staple(self):
        pass
        
    def lessen_staple(self):
        pass
        
    def modify_staple(self):
        pass


def config():
    '''读取配置文件，获取当前使用的类，并配置该类'''
    cluster_type = FixCluster #当前使用的类
    config = ConfigParser()
    if os.path.isfile('config.ini'):
        config.read('config.ini')
    else:
        config.read(os.path.join(os.path.split(__file__)[0],'config.ini'))

    #要使用的团簇类
    sym = config.get('cluster', 'symmetry')
    #assert sym in ['C1','Cs','C2','C3','C5']
    elements = config.get('cluster', 'element').split()
    if len(elements) > 2:
        if sym == 'C1':
            cluster_type = MultiCluster
        else:
            cluster_type = SymMultiCluster
    elif len(elements) == 2:
        if sym == 'C1':
            if get_element_id(elements[0]) == -1:
                cluster_type = Hydrate
            elif get_element_id(elements[1]) == -1:
                cluster_type = LigandCluster
            else:
                cluster_type = BiCluster
        elif sym == 'Cs':
            cluster_type = SymCsBiCluster
        elif sym == 'C2':
            cluster_type = SymC2BiCluster
        elif sym == 'C3':
            cluster_type = SymC3BiCluster
        elif sym == 'P1':
            from .nanowire import NanoWire2
            cluster_type = NanoWire2
        elif sym == 'A':
            cluster_type = Alternate
    else:
        if sym == 'C1':
            if get_element_id(elements[0]) != -1:
                cluster_type = MonoCluster
            else:
                cluster_type = MoleCluster
        elif sym == 'Cs':
            cluster_type = SymCsCluster
        elif sym == 'C2':
            cluster_type = SymC2Cluster
        elif sym == 'C3':
            cluster_type = SymC3Cluster
        elif sym == 'C5': #C5
            cluster_type = SymC5Cluster
        elif sym == 'Ci':
            cluster_type = SymCiCluster
        elif sym == 'P1':
            #cluster_type = NanoWire
            pass
        elif sym == 'B':
            cluster_type = Bilayer #临时
    cluster_type.config()
    return cluster_type
