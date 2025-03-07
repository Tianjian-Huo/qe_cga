'''point.py
三维点类
作者: 赛琳伟
创建时间: 2014-05-04
'''
import math
from random import random, uniform


class Point:
    '''三维点的类
    '''
    zero = 0.0001
    
    def __init__(self, xx=0, yy=0, zz=0):
        self.x = float(xx)
        self.y = float(yy)
        self.z = float(zz)
        
    def set_coordinate(self, p):
        '''赋值会改变原Point的指针，若不想改变用该函数'''
        self.x = p.x
        self.y = p.y
        self.z = p.z

    def __repr__(self):
        '''转字符串
        用法: str(p), print p'''
        return '%f\t%f\t%f' %(self.x, self.y, self.z)

    def __eq__(self, p):
        '''==运算符重载
        '''
        return math.sqrt((self.x-p.x)*(self.x-p.x) + (self.y-p.y)*(self.y-p.y) + (self.z-p.z)*(self.z-p.z)) < Point.zero

    def __add__(self, p):
        '''+运算符重载
        '''
        return Point(self.x+p.x, self.y+p.y, self.z+p.z)

    def __iadd__(self, p):
        '''+=运算符重载
        '''
        self.x += p.x
        self.y += p.y
        self.z += p.z
        return self

    def __sub__(self, p):
        '''-运算符重载
        '''
        return Point(self.x-p.x, self.y-p.y, self.z-p.z)

    def __isub__(self, p):
        '''-=运算符重载
        '''
        self.x -= p.x
        self.y -= p.y
        self.z -= p.z
        return self

    def __mul__(self, c):
        '''*运算符重载
        '''
        return Point(c*self.x, c*self.y, c*self.z)

    def __rmul__(self, c):
        '''*运算符重载,数在前,点在后.反射乘法
        '''
        return Point(c*self.x, c*self.y, c*self.z)

    def __imul__(self, c):
        '''*=运算符重载
        '''
        self.x *= c
        self.y *= c
        self.z *= c
        return self

    def __truediv__(self, c):
        '''/运算符重载
        '''
        return Point(self.x/c, self.y/c, self.z/c)

    def __idiv__(self, c):
        '''/=运算符重载
        '''
        self.x /= c
        self.y /= c
        self.z /= c
        return self

    def __neg__(self):
        '''-(负号)运算符重载
        '''
        return Point(-self.x, -self.y, -self.z)

    def dot(self, p):
        '''内积'''
        return p.x*self.x + p.y*self.y + p.z*self.z

    def cross(self, p):
        '''外积'''
        return Point(self.y*p.z - self.z*p.y,
                     self.z*p.x - self.x*p.z,
                     self.x*p.y - self.y*p.x)

    def mix_prod(self, p, q):
        '''混合积: p×q·self'''
        return self.x * (p.y*q.z - p.z*q.y) +\
               self.y * (p.z*q.x - p.x*q.z) +\
               self.z * (p.x*q.y - p.y*q.x)
               
    def plane(self):
        '''返回平面点'''
        return Point(self.x, self.y)
    
    def norm(self):
        '''模'''
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def angle(self, p):
        '''夹角'''
        return math.acos(self.dot(p) / (self.norm()*p.norm()+1e-13))

    def angle_x(self):
        '''与x轴的夹角'''
        return math.acos(self.x / self.norm())

    def angle_y(self):
        '''与y轴的夹角'''
        return math.acos(self.y / self.norm())

    def angle_z(self):
        '''与z轴的夹角'''
        return math.acos(self.z / self.norm())
    
    def dist(self, p):
        return math.sqrt((self.x-p.x)**2 + (self.y-p.y)**2 + (self.z-p.z)**2)

    def area(self, p, q):
        '''self, p, q三点组成的面积'''
        return (p-self).cross(q-self).norm() / 2

    def volume(self, p, q):
        '''三棱锥体积'''
        return (p-self).mid_prod(q-self, q-self) / 6

    def unit(self):
        '''单位向量'''
        return self / self.norm()

    def rotate(self, theta, n, O=None):
        '''以O点为中心,n为法向,旋转theta'''
        if O is None:
            O = Point(0)
        n_unit = n.unit()
        prj = n_unit * (self-O).dot(n_unit)
        r = self - O - prj
        return r * math.cos(theta) - r.cross(n_unit) * math.sin(theta)  + prj + O

    @classmethod
    def random(cls, r=1.):
        '''类方法,返回模为1的随机点
        使用方式: a=Point.random()'''
        theta = random() * 2*math.pi
        phi =  uniform(-math.pi/2, math.pi/2)
        return Point(math.cos(phi) * math.sin(theta),
                         math.cos(phi) * math.cos(theta),
                         math.sin(phi)) * r

    @classmethod
    def random_plane(cls, r=1.):
        '''类方法,返回模为1的随机点
        使用方式: a=Point.random()'''
        theta = random() * 2*math.pi
        return Point(math.cos(theta), math.sin(theta)) * r
