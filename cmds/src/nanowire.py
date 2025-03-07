'''nanowire.py
该类实现了单质纳米线有关的基本操作
作者: 赛琳伟
创建时间: 2014-11-13
'''
from .cluster import Cluster
from .crystal import Crystal
from bicluster import BiCluster
from .point import Point
from .utils import write_log
from .atom import bond as BOND
from .abinitio import abinitio
from math import pi,sqrt
from random import random,randint,uniform,randrange
from copy import deepcopy
from configparser import ConfigParser
import os
import inspect
from functools import reduce


class NanoWire(Crystal):
    '''一元纳米线类.
    线的方向为z轴'''
    _element = 'C' #元素名
    _bond_len = BOND(_element, _element) #键长
    _max_atom_num = 30 #最大原子数
    _r = 1.


    def __init__(self):
        '''初始化'''
        Crystal.__init__(self)
        self.alpha = self.beta = self.gamma = 90.
        self.group = 'P1'


    @staticmethod
    def set_element(element):
        NanoWire._element = element
        NanoWire._bond_len = BOND(NanoWire._element, NanoWire._element) #键长


    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = ConfigParser()
        if os.path.isfile('config.ini'):
            cfg.read('config.ini')
        else:
            cfg.read(os.path.dirname(inspect.stack()[0][1])+'/config.ini')
        cls.set_element(cfg.get('cluster', 'element'))


    def read(self, file_name):
        Cluster.read(self, file_name)
        #self.simplify()


    def read_file_obj(self, f):
        line = f.readline().split() #读取第一行晶格常数
        assert len(line) == 7
        self.a = float(line[0])
        self.b = float(line[1])
        self.c = float(line[2])
        Cluster.read_file_obj(self, f)
        #self.simplify()


    def write_car(self, file_name):
        '''写car文件
        为了便于观察,将结构移动到晶胞中心'''
        c = self.deepcopy() #再定义一个变量,不改变自身坐标
        c.center()
        c.move(Point(self.a, self.b, self.c)/2)
        Crystal.write_car(c, file_name)


    def simplify(self):
        max_z = max([p.z for p in self.coordinate])
        self.calc_inertia()
        while True:
            for n in [2,3,5,7]:
                d = self.c / n
                c = deepcopy(self)
                for p in c.coordinate:
                    p.z += d
                    if p.z > max_z:
                        p.z -= self.c
                c.calc_inertia()
                if self.isomorphic(c):
                    break
            else:
                break
            write_log('simplify')
            self.coordinate = []
            for p in c.coordinate:
                if 0 <= p.z < d:
                    self.coordinate.append(p)


    def calc_cell(self):
        '''计算纳米线晶胞尺寸'''
        self.a = max([p.x for p in self.coordinate]) - \
            min([p.x for p in self.coordinate]) + 10
        self.b = max([p.y for p in self.coordinate]) - \
            min([p.y for p in self.coordinate]) + 10
        self.c = self._merge(self.coordinate, self.coordinate, NanoWire._bond_len, False, False).z


    def random(self):
        '''产生随机结构, 在圆柱形区域内'''
        atom_num = randint(2, NanoWire._max_atom_num) #原子个数
        self.c = sqrt(3)/8 * atom_num / self._r**2
        self.element = [NanoWire._element] * atom_num

        for i in range(2, atom_num):
            for it in range(100): #尝试100次
                p = Point.random_plane(NanoWire._r)
                p.z = uniform(0, self.c)
                for q in self.coordinate:
                    if p.dist(q) < Cluster.tolerance_min:
                        break
                    if p.dist(q+Point(0,0,self.c)) < Cluster.tolerance_min:
                        break
                    if p.dist(q+Point(0,0,-self.c)) < Cluster.tolerance_min:
                        break
                else:
                    self.coordinate.append(p)
                    break
            else:
                self.c += 0.2

        self.zoom(NanoWire._bond_len)
        self.a = self.b = NanoWire._r*2 + 10
        self.calc_cell()


    def perturbation(self, step=0.02):
        '''给结构微小扰动
        平均每个原子扰动100次,每次扰动不超过0.02'''
        branch = self.branchs() #连通分支个数
        for it in range(self.get_size()*100):
            index = randint(0, self.get_size()-1)
            old_pos = deepcopy(self.coordinate[index])
            self.coordinate[index] += Point.random(step*NanoWire._bond_len)
            #移动后的连通分支个数不得变少
            if self.branchs() > branch:
                self.coordinate[index] = old_pos
                continue
            #两原子不可太近
            for i in range(self.get_size()):
                if i == index:
                    continue
                if self.coordinate[index].dist(self.coordinate[i]) < Cluster.tolerance_min * NanoWire._bond_len or \
                    self.coordinate[index].dist(Point(self.coordinate[i].x,self.coordinate[i].y,self.coordinate[i].z+self.c)) < Cluster.tolerance_min * NanoWire._bond_len:
                    self.coordinate[index] = old_pos
                    break


    def mating(self, father, mother):
        '''杂交(串起来). 沿z轴,子代原子数不定.'''
        f = father.deepcopy()
        f.center()
        f.coordinate.sort(key=lambda p: p.z)
        f.coordinate = f.coordinate[int(uniform(f.c/4, f.c/2)):]

        m = deepcopy(mother)
        m.center()
        m.coordinate.sort(key=lambda p: p.z)
        m.rotate(0, 0, 2*pi*random())
        m.coordinate = mother.coordinate[:int(uniform(m.c/4, m.c/2))+1]
        self._merge(f.coordinate, m.coordinate, NanoWire._bond_len)
        self.element = [NanoWire._element] * self.get_size()
        self.calc_cell()


    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for it in range(100):
            src = randrange(0, self.get_size())
            dest = randrange(0, self.get_size())
            dest_pos = self.coordinate[dest] + Point.random(NanoWire._bond_len)
            if dest_pos.dist(self.coordinate[src]) < 2 * NanoWire._bond_len:
                continue
            for i in range(self._atom_num):
                if i != src and self.coordinate[i].dist(dest_pos) < \
                    Cluster.tolerance_min * NanoWire._bond_len:
                    break
            else:
                break

        self.coordinate[src] = dest_pos
        return src


    def mutation_more(self):
        best = self.deepcopy()
        best.set_energy(0.)
        for i in range(5):
            c = self.deepcopy()
            for j in range(randint(1,3)):
                active = []
                move = c.mutation_one()
                if not move in active:
                    active.append(move)
            c.write_fix_input(active) #写input文件
            abinitio(c, input_name='optimize_fix')
            if c.get_energy() < best.get_energy():
                for i in range(self.get_size()):
                    if c.coordinate[i].dist(self.coordinate[i]) > 2:
                        best = c.deepcopy()
                        break
        self = best.deepcopy()




class NanoWire2(NanoWire):
    '''二元纳米线类.
    线的方向为z轴'''
    _element2 = None #元素2名
    _atom_ratio = 1. #原子比例
    _bond_len11 = 1. #键长
    _bond_len12 = 1.
    _bond_len22 = 1.
    _simple1 = 0
    _simple2 = 0
    _max_atom_num = 13 #最大原子数
    _max_diameter = 5.5 #最大直径


    def __init__(self):
        '''初始化'''
        NanoWire.__init__(self)
        self._atom_num1 = self._atom_num2 = 0


    @staticmethod
    def set_element(element1, element2):
        NanoWire2._element, NanoWire2._element2 = element1, element2
        NanoWire2._bond_len11 = BOND(NanoWire2._element, NanoWire2._element)
        NanoWire2._bond_len12 = BOND(NanoWire2._element, NanoWire2._element2)
        NanoWire2._bond_len22 = BOND(NanoWire2._element, NanoWire2._element2)


    @staticmethod
    def set_ratio(natom1, natom2):
        NanoWire2._simple1 = natom1
        NanoWire2._simple2 = natom2


    @classmethod
    def config(cls):
        '''读取配置文件'''
        cfg = ConfigParser()
        if os.path.isfile('config.ini'):
            cfg.read('config.ini')
        else:
            cfg.read(os.path.dirname(inspect.stack()[0][1])+'/config.ini')
        e = cfg.get('cluster', 'element').split()
        if len(e) != 2:
            print('Please input two elements in config.ini.')
            raise
        cls.set_element(e[0], e[1])
        cls.set_element(e[0],e[1])
        atom_num = cfg.get('cluster', 'atom_num').split()
        if len(atom_num) != 2:
            print('Please input two atom num in config.ini.')
            raise
        cls.set_ratio(int(atom_num[0]), int(atom_num[1]))


    def deepcopy(self, other=None):
        if other is None:
            result = Crystal.deepcopy(self)
            result._atom_num1 = self._atom_num1
            result._atom_num2 = self._atom_num2
            return result
        else:
            Crystal.deepcopy(self, other)
            self._atom_num1 = other._atom_num1
            self._atom_num2 = other._atom_num2


    def calc_inertia(self):
        '''计算转动惯量'''
        self._inertia = 0.
        O2 = reduce(Point.__add__, self.coordinate[self._atom_num1:]) / self._atom_num2
        self._inertia = reduce(lambda x,y:x+y, [(p-O2).dot(p-O2) for p in self.coordinate[:self._atom_num1]])

        self._inertia2 = 0.
        O1 = reduce(Point.__add__, self.coordinate[:self._atom_num1]) / self._atom_num1
        self._inertia2 = reduce(lambda x,y:x+y, [(p-O1).dot(p-O1) for p in self.coordinate[self._atom_num1:]])


    def read(self, file_name):
        NanoWire.read(self, file_name)
        self._atom_num1 = self.element.count(NanoWire2._element)
        self._atom_num2 = self.element.count(NanoWire2._element2)
        assert self._atom_num1*self._simple2 == self._atom_num2*self._simple1


    def read_file_obj(self, f):
        NanoWire.read_file_obj(self, f)
        self._atom_num1 = self.element.count(NanoWire2._element)
        self._atom_num2 = self.element.count(NanoWire2._element2)


    def diameter(self):
        max_bond_len = max(NanoWire2._bond_len11, NanoWire2._bond_len12, NanoWire2._bond_len22)
        d = 0.
        for p in self.coordinate:
            for q in self.coordinate:
                if abs(p.z-q.z) < max_bond_len or abs(p.z-q.z) > self.c-max_bond_len:
                    d = max(d, sqrt((p.x-q.x)**2 + (p.y-q.y)**2))
        return d


    def calc_cell(self):
        '''计算纳米线晶胞尺寸'''
        self.a = max([p.x for p in self.coordinate]) - \
            min([p.x for p in self.coordinate]) + 10
        self.b = max([p.y for p in self.coordinate]) - \
            min([p.y for p in self.coordinate]) + 10
        max_bond_len = max(NanoWire2._bond_len11, NanoWire2._bond_len12, NanoWire2._bond_len22)
        self.c = self._merge(self.coordinate, self.coordinate, max_bond_len*Cluster.tolerance_min, False, False).z


    def best_cell(self):
        self.calc_cell()
        a = self.c - 0.7
        b = self.c + 0.1
        x1, x2 = None, None

        while b-a < 0.01:
            if x1 is None:
                x1 = b - (b-a)*0.618
                tmp = self.deepcopy()
                tmp.c = x1
                abinitio(tmp)
                e1 = tmp.get_energy()
                print(x1,e1)
                if tmp.diameter() > self._max_diameter:
                    b, x2 = x2, None
                    continue
            if x2 is None:
                x2 = a + (b-a)*0.618
                tmp = self.deepcopy()
                tmp.c = x2
                abinitio(tmp)
                e2 = tmp.get_energy()
                print(x2,e2)
                if tmp.diameter() > self._max_diameter:
                    b, x2 = x2, None
                    continue
            if e1 < e2 - 0.001:
                a, x1 = x1, None
            elif e2 < e1 - 0.001:
                b, x2 = x2, None
            else:
                a, x1 = x1, None
                b, x2 = x2, None
            self.deepcopy(tmp)


    def random(self):
        '''产生随机结构, 在圆柱形区域内'''
        n = randint(2, NanoWire2._max_atom_num/(NanoWire2._simple1+NanoWire2._simple2))
        BiCluster.set_element(NanoWire2._element, NanoWire2._element2)
        BiCluster.set_atom_num(n*NanoWire2._simple1, n*NanoWire2._simple2)
        c = BiCluster()
        c.random()
        self.coordinate = c.coordinate
        self.element = c.element
        self._atom_num1, self._atom_num2 = n*NanoWire2._simple1, n*NanoWire2._simple2
        self.best_cell()


    def perturbation(self, step=0.05):
        '''给结构微小扰动
        平均每个原子扰动100次,每次扰动不超过0.05'''
        n = self.get_size()
        self.coordinate = self.coordinate + [p+Point(0,0,self.c) for p in self.coordinate]
        self.element = self.element * 2
        count=0
        for it in range(n*100):
            idx = randrange(n)
            old_pos = deepcopy(self.coordinate[idx])
            self.coordinate[idx] += Point.random(step)
            for i in range(self.get_size()):
                if i != idx and not self._is_legal(i, idx):
                    self.coordinate[idx] = old_pos
                    count+=1
                    break
        self.coordinate = self.coordinate[:n]
        self.element = self.element[:n]


    def mating(self, father, mother):
        '''杂交(串起来). 沿z轴,子代原子数不定.'''
        f = Cluster()
        f.coordinate = deepcopy(father.coordinate)
        f.center()
        #AABB
        f.coordinate = f.coordinate[:father._atom_num1] + \
            [p+Point(0,0,father.c) for p in f.coordinate[:father._atom_num1]] + \
            f.coordinate[father._atom_num1:] + \
            [p+Point(0,0,father.c) for p in f.coordinate[father._atom_num1:]]

        m = Cluster()
        m.coordinate = deepcopy(mother.coordinate)
        m.center()
        m.coordinate = m.coordinate[:mother._atom_num1] + \
            [p+Point(0,0,mother.c) for p in m.coordinate[:mother._atom_num1]] + \
            m.coordinate[mother._atom_num1:] + \
            [p+Point(0,0,mother.c) for p in m.coordinate[mother._atom_num1:]]

        max_bond_len = max(NanoWire2._bond_len11, NanoWire2._bond_len12, NanoWire2._bond_len22)
        c = (father.c+mother.c)/2

        for it in range(200):
            face_left = uniform(-father.c/2, father.c/2)  #杂交面
            face_right = face_left + uniform(0.5, father.c)
            up1 = [p for p in f.coordinate[:2*father._atom_num1] if face_left < p.z < face_right]
            up2 = [p for p in f.coordinate[2*father._atom_num1:] if face_left < p.z < face_right]
            m.rotate(0, 0, 2*pi*random())
            face_left = uniform(-mother.c/2, mother.c/2)  #杂交面
            face_right = face_left + uniform(0.5, mother.c)
            down1 = [p for p in m.coordinate[:2*mother._atom_num1] if face_left < p.z < face_right]
            down2 = [p for p in m.coordinate[2*mother._atom_num1:] if face_left < p.z < face_right]
            if len(up1)+len(up2) == 0 or len(down1)+len(down2) == 0 or len(up2)+len(down2) == 0 or \
                len(up1)+len(up2)+len(down1)+len(down2) > NanoWire2._max_atom_num:
                continue
            if (len(up1)+len(down1)) * NanoWire2._simple2 == (len(up2)+len(down2)) * NanoWire2._simple1:
                gap = self._merge(down1+down2, up1+up2, max_bond_len, False, False)
                self.coordinate = down1 + [p+gap for p in up1] + down2 + [p+gap for p in up2]
                self.element = [NanoWire2._element]*(len(up1)+len(down1)) + \
                               [NanoWire2._element2]*(len(up2)+len(down2))
                self.best_cell()
                self._atom_num1 = len(up1)+len(down1)
                self._atom_num2 = len(up2)+len(down2)
                break
        else:
            print('mating failure')
            write_log('mating failure')
            self = father.deepcopy()
            self.perturbation() #改用变异


    def mutation_one(self):
        '''将一个原子移动到另外一个原子旁边。返回移动的原子序号'''
        for it in range(100):
            src = randrange(0, self.get_size())
            dest = randrange(0, self.get_size())
            if src < self._atom_num1:
                if dest < self._atom_num1:
                    dest_pos = self.coordinate[dest] + Point.random(NanoWire2._bond_len11)
                else:
                    dest_pos = self.coordinate[dest] + Point.random(NanoWire2._bond_len12)
                if dest_pos.dist(self.coordinate[src]) < 2 * NanoWire2._bond_len11:
                    continue
                try:
                    for i in range(self._atom_num1):
                        if i != src and self.coordinate[i].dist(dest_pos) < \
                            Cluster.tolerance_min * NanoWire2._bond_len11:
                            raise
                    for i in range(self._atom_num2):
                        if self.coordinate[self._atom_num1+i].dist(dest_pos) < \
                            Cluster.tolerance_min * NanoWire2._bond_len11:
                            raise
                except:
                    continue
            else:
                if dest < self._atom_num1:
                    dest_pos = self.coordinate[dest] + Point.random(NanoWire2._bond_len12)
                else:
                    dest_pos = self.coordinate[dest] + Point.random(NanoWire2._bond_len22)
                if dest_pos.dist(self.coordinate[src]) < 2 * NanoWire2._bond_len22:
                    continue
                try:
                    for i in range(self._atom_num1):
                        if self.coordinate[i].dist(dest_pos) < \
                            Cluster.tolerance_min * NanoWire2._bond_len12:
                            raise
                    for i in range(self._atom_num2):
                        if i != src and self.coordinate[self._atom_num1+i].dist(dest_pos) < \
                            Cluster.tolerance_min * NanoWire2._bond_len22:
                            raise
                except:
                    continue

        self.coordinate[src] = dest_pos
        return src


    def mutation_more(self):
        '''移动多个原子'''
        n = randrange(3)+1
        for i in range(1, n):
            self.mutation_one()




def cluster2wire(c):
    wire = NanoWire()
    wire.coordinate = deepcopy(c.coordinate)
    wire.element = deepcopy(c.element)
    wire.calc_cell()
    return wire
