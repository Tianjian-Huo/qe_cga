'''element.py
元素周期表的类
作者: 赛琳伟
创建时间: 2014-08-17
该文件仅需被import一次'''
import numpy as np
import os
import inspect
from .point import Point

#元素表
_table = ["Nul","H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fi","Mc","Lv","Ts","Og"]


#元素颜色表,可以不按顺序写.查看某原子的颜色的方法是在MS中打开结构,选中某原子,在左边的Properties Explorer里有显示
_color_table = {'Nul':'#ffffff',
                'H':'#ffffff',
                'He':'#d9ffffff',

                'Li':'#cd7eff',
                'Be':'#c5ff00',
                'B':'#ffb7b7',
                'C':'#828282',
                'N':'#0000ff',
                'O':'#ff0000',
                'F':'#b3ffff',
                'Ne':'#afe3f4',

                'Na':'#aa5ef2',
                'Mg':'#89ff00',
                'Al':'#ff64ff',
                'Si':'#ffc828',
                'P':'#e99cff',
                'S':'#ffeb00',
                'Cl':'#b3ffad',
                'Ar':'#81d1e4',

                'Ti':'#c0c3c6',
                'Ge':'#66ad63',

                'W':'#2794d6',
                'Au':'#ffd124',
                'Ag':'#99c6ff'}

_path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
_path = os.path.join(_path, 'data/bond_table.npy')
if not os.path.isfile(_path):
    _path = 'bond_table.npy'
    if not os.path.isfile(_path):
        print('error, no bond_table.npy')
        raise
_bond_table = np.load(_path)



def get_element(idx):
    if idx <= 0 or idx > 112:
        print('No such element.')
        return None
    return _table[idx]


def get_element_id(elem):
    try:
        return _table.index(elem)
    except:
        return -1


def get_color(element):
    if element in _color_table:
        return _color_table[element]
    else:
        print('No such element.')
        return '#ffffff'


def add_bond(element1, element2, bond):
    '''添加键长
    参数1:element1,元素1,字符串类型
    参数2:element2,元素2,字符串类型
    参数3:bond,键长,实数类型
    返回:无'''
    element1 = _table.index(element1)
    element2 = _table.index(element2)
    _bond_table[element1, element2] = bond
    _bond_table[element2, element1] = bond
    _path = os.path.split(inspect.getfile(inspect.currentframe()))[0]
    #_bond_table = np.load(os.path.join(_path, 'bond_table.npy'))
    np.save(os.path.join(_path, 'bond_table.npy'), _bond_table)


def bond_list():
    '''显示bond_table.npy中存储的全部键长，打印在屏幕上
    返回:无'''
    for i in range(119):
        for j in range(i, 119):
            if _bond_table[i,j] != 0.:
                print(_table[i], _table[j], _bond_table[i,j])


def bond(element1, element2):
    '''查询键长
    参数1:element1,元素1,字符串类型
    参数2:element2,元素2,字符串类型
    返回:键长'''
    element1 = _table.index(element1)
    element2 = _table.index(element2)
    if _bond_table[element1, element2] == 0.:
        print('can not found this bond, please add first.')
    return _bond_table[element1, element2]


def set_bond(elem1, elem2, length):
    _bond_table[_table.index(elem1), _table.index(elem2)] = length



tolerance_min = 0.75
tolerance_max = 1.25


class Atom(Point):
    zero = 0.0001

    def __init__(self, xx=0, yy=0, zz=0, e=''):
        if isinstance(xx, Point) and get_element_id(yy) != -1:
            self.x = xx.x
            self.y = xx.y
            self.z = xx.z
            self.elem = yy
        else:
            self.x = float(xx)
            self.y = float(yy)
            self.z = float(zz)
            self.elem = e


    def __repr__(self):
        '''转字符串'''
        return self.elem + '\t' + Point.__repr__(self)


    def coordinate(self):
        return Point(self.x, self.y, self.z)


    def __eq__(self, other):
        '''==运算符重载
        '''
        return self.elem == other.elem and Point.__eq__(self, other)


    def __add__(self, a):
        '''+运算符重载
        '''
        return Atom(self.x+a.x, self.y+a.y, self.z+a.z, self.elem)


    def __sub__(self, a):
        '''-运算符重载
        '''
        return Atom(self.x-a.x, self.y-a.y, self.z-a.z, self.elem)


    def __mul__(self, c):
        '''*运算符重载
        '''
        return Atom(c*self.x, c*self.y, c*self.z, self.elem)


    def __rmul__(self, c):
        '''*运算符重载,数在前,点在后.反射乘法
        '''
        return Atom(c*self.x, c*self.y, c*self.z, self.elem)


    def __div__(self, c):
        '''/运算符重载
        '''
        return Atom(self.x/c, self.y/c, self.z/c, self.elem)


    def __neg__(self):
        '''-(负号)运算符重载
        '''
        return Atom(-self.x, -self.y, -self.z, self.elem)


    def cross(self, a):
        '''外积'''
        return Atom(self.y*a.z - self.z*a.y,
                     self.z*a.x - self.x*a.z,
                     self.x*a.y - self.y*a.x, self.elem)


    def plane(self):
        '''返回平面点'''
        return Atom(self.x, self.y, 0, self.elem)


    def bond(self, other):
        return self.dist(other) < tolerance_max * bond(self.elem, other.elem)


    def is_legal(self, other):
        return self.dist(other) > tolerance_min * bond(self.elem, other.elem)


    def exchange_coordinate(self, other):
        self.x, other.x = other.x, self.x
        self.y, other.y = other.y, self.y
        self.z, other.z = other.z, self.z
