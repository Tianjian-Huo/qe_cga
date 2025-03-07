'''产生特殊团簇结构
作者: 赛琳伟
版本: 1.0.3
创建时间: 2014-06-06
修改记录: 2014-06-09
'''
import sys
sys.path.append('../..')
from math import sqrt, sin, cos, pi
from cmds.src.point import Point
from cmds.src.atom import *
from cmds.src.cluster import Cluster
from cmds.src.optimize import Clean
#from crystal import Crystal
from copy import deepcopy
import numpy as np



def unique(c):
    '''去除重复元素'''
    result = Cluster()
    hash_ = set()

    for a in c.atoms:
        h = (int(a.x*10+0.5), int(a.y*10+0.5), int(a.z*10+0.5))
        if not h in hash_:
            result.add_atom(a)
            hash_.add(h)
    return result


def Ci(c):
    '''变为中心对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(-c.atoms[i], c.atoms[i].elem)
    return unique(result)


def sym_x(c):
    '''变为x轴对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(c.atoms[i].x, -c.atoms[i].y, -c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def sym_y(c):
    '''变为y轴对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(-c.atoms[i].x, c.atoms[i].y, -c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def sym_z(c):
    '''变为z轴对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(-c.atoms[i].x, -c.atoms[i].y, c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def sym_xy(c):
    '''变为xy平面对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(c.atoms[i].x, c.atoms[i].y, -c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def sym_xz(c):
    '''变为xz平面对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(c.atoms[i].x, -c.atoms[i].y, c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def sym_yz(c):
    '''变为yz平面对称'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        result.add_atom(Point(-c.atoms[i].x, c.atoms[i].y, c.atoms[i].z), c.atoms[i].elem)
    return unique(result)


def Cs(c):
    return sym_xy(c)


def C3(c):
    '''变为C3对称, 原子数最多变为3倍. x=y=z则只有1倍.'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        x, y, z = c.atoms[i].x, c.atoms[i].y, c.atoms[i].z
        result.add_atom(Point(y, z, x), c.atoms[i].elem)
        result.add_atom(Point(z, x, y), c.atoms[i].elem)
    return unique(result)


def Cn(c, n):
    result = deepcopy(c)
    for i in range(c.get_size()):
        for j in range(1,n):
            result.add_atom(c.atoms[i].rotate(2*pi*j/n, Point(0,0,1)), c.atoms[i].elem)
    return unique(result)


def Sn(c, n):
    result = deepcopy(c)
    for i in range(c.get_size()):
        for j in range(1,n):
            p = c.atoms[i].rotate(2*pi/n, Point(0,0,1))
            result.add_atom(Point(p.x,p.y,-p.z), c.atoms[i].elem)
    return unique(result)


def D2(c):
    '''变为D2对称, 原子数最多变为4倍. 在轴上则只有两倍'''
    result = deepcopy(c)
    for i in range(c.get_size()):
        x, y, z = c.atoms[i].x, c.atoms[i].y, c.atoms[i].z
        result.add_atom(Point(-x, -y, z), c.atoms[i].elem)
        result.add_atom(Point(-x, y, -z), c.atoms[i].elem)
        result.add_atom(Point(x, -y, -z), c.atoms[i].elem)
    return unique(result)


def D2h(c):
    '''变为D2h对称, 原子数最多变为8倍'''
    return Ci(D2(c))


def T(c):
    '''变为T对称, 原子数最多变为12倍'''
    return D2(C3(c))


def Td(c):
    '''变为Td对称, 24倍'''
    c2 = T(c)
    n = c2.get_size()
    for i in range(n):
        c2.add_atom(Point(c2.atoms[i].y, c2.atoms[i].x, c2.atoms[i].z), c2.atoms[i].elem)
    return unique(c2)


def Th(c):
    '''变为Th对称, 24倍'''
    return Ci(T(c))


def O(c):
    '''变为O对称, 24倍'''
    c2 = T(c)
    n = c2.get_size()
    for i in range(n):
        c2.add_atom(Point(c2.atoms[i].y, c2.atoms[i].x, -c2.atoms[i].z), c2.atoms[i].elem)
    return unique(c2)


def Oh(c):
    '''变为Oh对称, 48倍'''
    return Ci(O(c))


def I(c):
    '''变为I对称, 60倍'''
    axis5 = Atom(Point(0,(sqrt(5)-1)/2,1), 'Nul')
    result = Cluster()
    for a in c.atoms:
        for i in range(5):
            result.add_atom(a.rotate(2*pi*i/5, axis5))
    return T(unique(result))

def Ih(c):
    '''变为Ih对称, 最多120倍'''
    return Ci(I(c))


###############################################################################
#以下为产生特殊团簇的函数
###############################################################################
def ring(n, elem='H'):
    '''大小为n的环'''
    cluster = Cluster()
    r = 0.5 / sin(pi/n)

    for i in range(n):
        cluster.atoms.append(Point(r*cos(2*pi*i/n), r*sin(2*pi*i/n), 0))

    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def triangle(n, elem='H'):
    c = Cluster()
    c.atoms.append(Point(0))
    for i in range(1,n):
        for j in range(i+1):
            c.atoms.append(Point(-i/2.0+j, sqrt(3)/2*i))
    c.set_element(elem)
    c.zoom(bond(elem,elem))
    return c


def triangle_sheet6(a,b,c,na,elem='H'):
    '''边界为6条边的三角形网格
    横向来看，边长逐渐增大1，然后平稳，再减小1
    a：横向最长
    b：横向最上长度
    c:横向最下长度
    na：横向最大边平稳次数'''
    assert a>b>0 and a>c and na>0
    cl = Cluster()
    left_start = 0.
    for i in range(b,a):
        for j in range(i):
            cl.atoms.append(Point(left_start+j, sqrt(3)/2*(i-b)))
        left_start -= 0.5
    for i in range(na):
        for j in range(a):
            cl.atoms.append(Point(left_start+j, sqrt(3)/2*(a-b+i)))
        left_start -= 0.5
    left_start += 0.5
    for i in range(a-c):
        left_start += 0.5
        for j in range(a-i-1):
            cl.atoms.append(Point(left_start+j, sqrt(3)/2*(a-b+na+i)))
    cl.set_element(elem)
    cl.zoom(bond(elem,elem))
    return cl


def bending_boron(c):
    c.center()
    for a in c.atoms:
        a.z = a.norm()/5


def boron_planar_search():
    for na in range(1,6):
        for a in range(5,16):
            for b in range(2,a):
                for c in range(2,b+1):
                    cl = triangle_sheet6(a,b,c,na,'B')
                    bending_boron(cl)
                    print(na,a,b,c,cl.get_size())
                    cl.write_car('B_planar/B%d_planar%d%d%d%d.car'%(cl.get_size(),a,b,c,na))



###############################################################################
# 以下为5种正多面体
###############################################################################
def tetrahedron(elem='H'):
    '''正四面体
    原子数:4
    对称性:'''
    cluster = Cluster()
    cluster.add_atom(Point(1, 1, 1) * (sqrt(2)/4), elem)
    cluster = D2(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def cube(elem='H'):
    '''正方体
    原子数:6
    对称性:'''
    cluster = Cluster()
    cluster.add_atom(Point(1./2, 1./2, 1./2), elem)
    cluster = D2h(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def octahedron(elem='H'):
    '''正八面体
    原子数:6
    对称性:'''
    cluster = Cluster()
    cluster.add_atom(Point(1/sqrt(2), 0, 0), elem)
    cluster = C3(cluster)
    cluster = Ci(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def dodecahedron(elem='H'):
    '''正十二面体
    原子数:20
    对称性:'''
    cluster = Cluster()
    a = (1 + sqrt(5)) / 4
    cluster.add_atom(Point(a, a, a), elem)
    cluster = I(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster

    """cluster = Cluster()
    a = (1 + sqrt(5)) / 2
    cluster.add_atom(Point(0, 1./2, a*a/2), elem)
    cluster = T(cluster)

    cluster2 = Cluster()
    cluster2.add_atom(Point(1, 1, 1) * (a/2), elem)
    cluster2 = D2h(cluster2)

    cluster += cluster2
    cluster.zoom(bond(elem,elem))
    return cluster"""


def icosahedron(elem='H'):
    '''正二十面体
    原子数:12
    对称性:'''
    cluster = Cluster()
    cluster.add_atom(Point(0, 1./2, (1 + sqrt(5)) / 4), elem)
    cluster = T(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster



###############################################################################
# 以下为13种阿基米德多面体 和正棱柱,正反棱柱. 合成半正多面体
###############################################################################
def archimedes366(elem='H'):
    '''截角四面体
    原子数:12
    对称性: Td
    体积: 23/12*√2'''
    cluster = Cluster()
    cluster.add_atom(Point(sqrt(2)/4, sqrt(2)/4, 3*sqrt(2)/4), elem)
    cluster = T(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes3434(elem='H'):
    '''截半立方体
    原子数:12
    对称性: Oh
    体积: 5/3*√2'''
    cluster = Cluster()
    cluster.add_atom(Point(0, 1/sqrt(2), 1/sqrt(2)), elem)
    cluster = T(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes466(elem='H'):
    '''截角八面体
    原子数:24
    对称性: Oh
    体积: 8√2'''
    cluster = Cluster()
    cluster.add_atom(Point(0, 1/sqrt(2), sqrt(2)), elem)
    cluster = Td(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes388(elem='H'):
    '''截角立方体
    原子数:24
    对称性: Oh
    体积: (21 + 14√2) / 3'''
    cluster = Cluster()
    cluster.add_atom(Point((sqrt(2) + 1) / 2, (sqrt(2) + 1) / 2, 1./2), elem)
    cluster = O(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes3535(elem='H'):
    '''截半二十面体
    原子数:30
    对称性: Ih
    体积: (45 + 17√5) / 6'''
    cluster = Cluster()
    cluster.add_atom(Point((1 + sqrt(5)) / 4, 1./2, (1 + sqrt(5)) / 4 + 1./2), elem)
    cluster.add_atom(Point((1 + sqrt(5)) / 2, 0, 0), elem)
    cluster = Th(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes566(elem='H'):
    '''截角二十面体
    原子数:60
    对称性: Ih
    体积: (125 + 43√5) / 4
    外接球半径: (√(9(√5+1)/2 + 10)) / 2'''
    cluster = Cluster()
    a = (sqrt(5) + 1) / 2
    cluster.add_atom(Point(0, 1./2, 3*a/2), elem)
    cluster.add_atom(Point(1, 1./2+a, a/2), elem)
    cluster.add_atom(Point(1./2, 1+a/2, a), elem)
    cluster = Th(cluster)
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes3444(elem='H'):
    '''小斜方截半立方体
    原子数:24
    对称性: Oh
    体积: (12 + 10√2) / 3'''
    cluster = Cluster()
    cluster.add_atom(Point(1./2, 1./2, (sqrt(2) + 1) / 2), elem)
    cluster = Th(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes33334(elem='H'):
    '''扭棱立方体
    原子数:24
    对称性: O
    体积:
    有手性'''
    cluster = Cluster(elem='H')
    # xi^3 + xi^2 + xi = 1的根
    xi = ((17 + 3*sqrt(33)) ** (1./3) - (-17 + 3*sqrt(33)) ** (1./3) - 1) / 3
    cluster.add_atom(Point(1, xi, 1/xi) / sqrt(2*(xi-1/xi/xi)*(xi-1)), elem)
    cluster = O(cluster)
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes3aa(elem='H'):
    '''截角十二面体, 3_10_10
    原子数:60
    对称性: Ih
    体积: (99 + 47√5) * 5/12'''
    cluster = Cluster()
    a = (sqrt(5) + 1) / 2
    cluster.add_atom(Point(0, a-1, a+2), elem)
    cluster.add_atom(Point(a, 2*a, a-1), elem)
    cluster.add_atom(Point(a, 2, a+1), elem)
    cluster = Th(cluster)
    cluster.zoom(1/(2*a-2))
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes3454(elem='H'):
    '''小斜方截半二十面体rhombicosidodecahedron
    原子数:60
    对称性: Ih
    体积: (60 + 29√5) / 3'''
    cluster = Cluster()
    a = (sqrt(5) + 1) / 2
    cluster.add_atom(Point(1, 1, 2*a+1), elem)
    cluster.add_atom(Point(a, 2*a, a+1), elem)
    cluster.add_atom(Point(0, a+1, a+2), elem)
    cluster = Th(cluster)
    cluster.zoom(1./2)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes468(elem='H'):
    '''大斜方截半立方体
    原子数:48
    对称性: Oh
    体积: 22 + 14√2'''
    cluster = Cluster()
    cluster.add_atom(Point(1./2, (1 + sqrt(2)) / 2, (1 + 2*sqrt(2)) / 2), elem)
    cluster = Oh(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes33335(elem='H'):
    '''扭棱十二面体
    原子数:60
    对称性: I
    体积: (12xi^2(3tau+1) - xi(36tau+7) - (53tau+6)) / (6sqrt(3-xi^2)^2)
    有手性'''
    cluster = Cluster()
    tau = (1 + sqrt(5)) / 2
    # xi^3 - 2xi = tau的根
    xi = (tau/2 + sqrt(tau-5./27)/2) ** (1./3) + (tau/2 - sqrt(tau-5./27)/2) ** (1./3)
    a = xi - 1/xi
    b = xi*tau + tau*tau + tau/xi
    p = Point(a + b/tau + tau, -a*tau + b + 1/tau, a/tau + b*tau - 1)
    q = Point(a + b/tau - tau, a*tau - b + 1/tau, a/tau + b*tau + 1)
    r = Point(-a/tau + b*tau + 1, -a + b/tau - tau, a*tau + b - 1/tau)
    s = Point(-a/tau + b*tau - 1, a - b/tau - tau, a*tau + b + 1/tau)

    cluster.add_atom(Point(2*a, 2, 2*b), elem)
    cluster.add_atom(p, elem)
    cluster.add_atom(q, elem)
    cluster.add_atom(r, elem)
    cluster.add_atom(s, elem)

    cluster = T(cluster)
    cluster.zoom(1/p.dist(r))
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def archimedes46a(elem='H'):
    '''大斜方截半二十面体, 4_6_10
    原子数:120
    对称性: Ih
    体积: 95 + 50√5
    有手性'''
    cluster = Cluster()
    a = (1 + sqrt(5)) / 2
    cluster.add_atom(Point(1/a, 1/a, 3+a), elem)
    cluster.add_atom(Point(a, 2*a+1, 2/a), elem)
    cluster.add_atom(Point(a*a, 3*a-1, 1/a), elem)
    cluster.add_atom(Point(2, 2+a, 2*a-1), elem)
    cluster.add_atom(Point(a, 3, 2*a), elem)
    cluster = Th(cluster)
    cluster.zoom(1/(2*a-2))
    cluster.zoom(bond(elem,elem))
    return cluster


def prism(n, elem='H'):
    '''正棱柱
    '''
    assert n > 2
    n/=2

    polygon = ring(n, 'Nul')
    cluster = deepcopy(polygon)
    cluster.move(Point(0, 0, 0.5))
    polygon.move(Point(0, 0, -0.5))
    cluster += polygon
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def antiprism(n, elem='H'):
    '''正反棱柱
    '''
    assert n > 2
    n/=2

    polygon = ring(n,'')
    cluster = deepcopy(polygon)
    h = sqrt(1 - (1 - cos(pi/n)) / (2*sin(pi/n)*sin(pi/n))) / 2
    cluster.move(Point(0, 0, h))
    polygon.move(Point(0, 0, -h))
    polygon.rotate(0, 0, pi/n)
    cluster += polygon
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster



###############################################################################
# 以下为13种卡塔兰多面体, 它们是13种阿基米德多面体的对偶
###############################################################################
def rhombic_dodecahedron(elem='H'):
    '''菱形十二面体,截半立方体(3434)的对偶
    原子数:24
    对称性: O
    体积: 16/9*√3'''
    cluster = Cluster()
    cluster.add_atom(Point(1, 1, 1) / sqrt(3), elem)
    cluster.add_atom(Point(2, 0, 0) / sqrt(3), 'N')
    cluster = O(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def deltoidal_icositetrahedron(elem='H'):
    '''四角化二十四面体,小斜方截半立方体(3444)的对偶
    原子数:24
    对称性: Oh
    两种键长分别为:0.92,0.75'''
    cluster = Cluster()
    a = sqrt(2)/2 + 1./2
    b = sqrt(2)/4 + 1./2
    c = sqrt(2)/6 + 1./2
    cluster.add_atom(Point(0, 0, a), elem)
    cluster.add_atom(Point(0, b, b), 'N')
    cluster.add_atom(Point(c, c, c), elem)
    cluster = Th(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def deltoidal_hexecontahedron(elem='H'):
    '''五角化六十面体,小斜方截半二十面体(3454)的对偶
    原子数:62
    对称性: Ih
    两种键长分别为:sqrt(6*(a+2)) / 3, sqrt(2*(a+1))
    即1.553, 2.288'''
    cluster = Cluster()
    a = (1 + sqrt(5)) / 2
    cluster.add_atom(Point(0, 0, 2*a+1), 'N')
    cluster.add_atom(Point(0, a/3+1, (5*a+4)/3), elem)
    cluster.add_atom(Point(0, (9*a+3)/5, (3*a+6)/5), elem)
    cluster.add_atom(Point(1, 1, 1) * ((4*a+1)/3), elem)
    cluster.add_atom(Point((a+1)/2, a+1./2, 3*a/2+1), 'N')
    cluster = Th(cluster)
    cluster.zoom(bond(elem,elem))
    return cluster


def rhombic_triacontahedron(elem='H'):
    '''菱形三十面体,截半二十面体(3535)的对偶
    原子数:32
    对称性: Ih
    体积: 4*sqrt(5+2*sqrt(5))'''
    cluster = Cluster()
    a = (1 + sqrt(5)) / 2
    cluster.add_atom(Point(0, (a+2)/5, (3*a+1)/5), 'N')
    cluster.add_atom(Point(0, (2*a+1)/3, a/3), elem)
    cluster.add_atom(Point(1, 1, 1) * ((a+1)/3), elem)
    cluster = Th(cluster)
    cluster.zoom(1 / sqrt((a+1)/3))
    cluster.zoom(bond(elem,elem))
    return cluster


def C60():
    c = Cluster([Point( 3.05829,    1.636317,    0.757489)], ['C'])
    return I(c)


def H2O3454(elem='H'):
    '''小斜方截半二十面体rhombicosidodecahedron
    原子数:60
    对称性: Ih
    体积: (60 + 29√5) / 3'''
    c1 = Cluster()
    a = (sqrt(5) + 1) / 2
    p1=Point(1, 1, 2*a+1)
    p2=Point(a, 2*a, a+1)
    p3=Point(0, a+1, a+2)
    p4 = p1+p2-p3
    q1=(p3*0.97+p2*1.8)/2.77
    q2=(p1*0.97+p3*1.8)/2.77
    q3=(p3*0.97+p1*1.8)/2.77
    q4=(p2*0.97+p3*1.8)/2.77
    q5=(p4*0.97+p1*1.8)/2.77
    q6=(p1*0.97+p4*1.8)/2.77
    q7=(p4*0.97+p2*1.8)/2.77
    q8=(p2*0.97+p4*1.8)/2.77
    p5=Point(0,3.62598,5.01098)
    q9=(p5*0.97+p1*1.8)/2.77
    q10=(p1*0.97+p5*1.8)/2.77
    c1.add_atom(p1, 'O')
    c1.add_atom(p2, 'O')
    c1.add_atom(p3, 'O')
    c1 = Th(c1)

    c2 = Cluster()
    c2.add_atom(q1, elem)
    c2.add_atom(q2, elem)
    c2.add_atom(q3, elem)
    c2.add_atom(q4, elem)
    c2.add_atom(q5, elem)
    c2.add_atom(q6, elem)
    c2.add_atom(q7, elem)
    c2.add_atom(q8, elem)
    c2.add_atom(q9, elem)
    c2.add_atom(q10, elem)
    c2 = Th(c2)

    c = c1+c2
    c.zoom(1./2)
    c.zoom(2.77)
    c.write_xyz('H20_3454.xyz')
    return c


def tube(diameter, length, elem='H'):
    '''
    '''
    polygon = ring(diameter, 'Nul')
    cluster = deepcopy(polygon)
    h = sqrt(1 - (1 - cos(pi/diameter)) / (2*sin(pi/diameter)*sin(pi/diameter)))
    for i in range(1,length):
        polygon.move(Point(0, 0, h))
        polygon.rotate(0, 0, pi/diameter)
        cluster += deepcopy(polygon)
    cluster.set_element(elem)
    cluster.zoom(bond(elem,elem))
    return cluster


def multi_cube(n, elem='H'):
    '''多层立方体'''
    c = Cluster()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c.add_atom(Point(i,j,k), elem)
    c.zoom(bond(elem,elem))
    return c


def multi_tetrahedron(n, elem='H'):
    cluster = Cluster()
    p = Point(1,-1,-1)
    q = Point(-1,-1,1)
    r = Point(-1,1,-1)
    s = Point(1,1,1)
    cluster.add_atom(s, elem)
    for i in range(1,n+1):
        pp = (s*(n-i) + p*i) / n
        qq = (s*(n-i) + q*i) / n
        rr = (s*(n-i) + r*i) / n
        cluster.add_atom(pp, elem)
        for j in range(1,i+1):
            f = (pp*(i-j)+qq*j)/i
            cluster.add_atom(f, elem)
            for k in range(1,j+1):
                cluster.add_atom(f+(rr-qq)*k/i, elem)
    cluster.zoom(n/sqrt(2)/2*bond(elem,elem))
    return cluster

def multi_icosahedron(n_layer, elem):
    result = Cluster()
    result.add_atom(Point(0), elem)
    ico = icosahedron(elem)
    P,Q,R = ico.atoms[1],ico.atoms[2],ico.atoms[11]
    result += ico
    for l in range(2, n_layer+1):
        shell = Cluster()
        for i in range(l+1):
            for j in range(l+1-i):
                a = i*P + j*Q + (l-i-j)*R
                shell.add_atom(a, elem)
        result += Ih(shell)
    return result


#以下废弃
def read_contcar(file_name):
    def is_cluster():
        x_len = max([p.x for p in clst.coordinate]) - min([p.x for p in clst.coordinate])
        y_len = max([p.y for p in clst.coordinate]) - min([p.y for p in clst.coordinate])
        z_len = max([p.z for p in clst.coordinate]) - min([p.z for p in clst.coordinate])
        return a-x_len>5 and b-y_len>5 and c-z_len>5

    clst = Cluster()
    a,b,c = clst.read_contcar(total_path)
    if is_cluster():
        return clust
    cry = Crystal()
    cry.atoms = deepcopy(c.atoms)
    cry.a,cry.b,cry.c = a,b,c
    return cry


def cluster2crystal(c):
    cry = Crystal()
    cry.atoms = deepcopy(c.atoms)
    d = c.diameter()
    cry.a = cry.b = cry.c = d+10
    cry.center()
    return cry


def read_adj(file_obj):
    n = int(file_obj.readline())
    adj = np.zeros((n,n),'int')
    for i in range(n):
        line = file_obj.readline().split()
        for v in line:
            adj[i,v] = 1
    return adj


def read_adj_list(file_obj):
    n = int(file_obj.readline())
    adj_list = []
    for i in range(n):
        line = file_obj.readline().split()
        adj_list.append([int(v) for v in line])
    return adj_list


def adj_list2adj(adj_list):
    adj = np.zeros((len(adj_list),len(adj_list)),'int')
    for i in range(len(adj_list)):
        for v in adj_list[i]:
            adj[i,v] = 1
    return adj


def adj_list2cluster(adj_list):
    c = Cluster()
    current = [0]
    remain = list(range(1, len(adj_list)))
    z = 0
    c.atoms.append(Point(0))
    while len(remain) != 0:
        new_current = []
        for p in current:
            for v in adj_list[p]:
                if v in remain:
                    new_current.append(v)
                    remain.remove(v)
        z += 2
        current = new_current
        r = sin(pi/len(current))
        for i in range(len(current)):
            c.atoms.append(Point(r*cos(2*pi*i/len(current)), r*sin(2*pi*i/len(current)), z))
    return c


def adj2cluster(file_name, elem):
    f = open(file_name)
    adj_list = read_adj_list(f)
    c = adj_list2cluster(adj_list)
    c.set_element(elem)
    c.write_xsd('gen.xsd')
    clean = Clean(c, adj=adj_list2adj(adj_list))
    clean.run()
    return c


def hexagon_abc(n1,n2,n3,elem='H'):
    '''abcabc
    作废，并入triangle_sheet6'''
    n1,n2,n3 = sorted([n1,n2,n3],reverse=True)
    c = Cluster()
    x, y = 1./2, sqrt(3)/2
    for i in range(n3):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+i):
            c.add_atom(Point(x+j, y), elem)
    for i in range(n2-n3):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+n3-1):
            c.add_atom(Point(x+j, y), elem)
    for i in range(1,n3):
        x += 1./2
        y -= sqrt(3)/2
        for j in range(n1+n3-1-i):
            c.add_atom(Point(x+j, y), elem)
    c.zoom(bond(elem,elem))
    return c


def hexagon_ab(n1,n2,elem):
    '''ababab
    作废，并入triangle_sheet6'''
    c = Cluster()
    x, y = 1./2, sqrt(3)/2
    for i in range(n2):
        x -= 1./2
        y -= sqrt(3)/2
        for j in range(n1+i):
            c.add_atom(Point(x+j, y), elem)
    for i in range(1,n1):
        x += 1./2
        y -= sqrt(3)/2
        for j in range(n1+n2-1-i):
            c.add_atom(Point(x+j, y), elem)
    c.zoom(bond(elem,elem))
    return c




if __name__ == '__main__':
    c = Cluster()
    a = (1 + sqrt(5)) / 4
    c.add_atom(Point(a, a, a), 'B')
    c = I(c)
    c.write_xyz('I.xyz')
    raise

    c = Cluster()
    a = (1 + sqrt(5)) / 2
    c.add_atom(Point(0, 1./2, a*a/2), 'B')
    c=I(c)
    c.zoom(1.8)
    c.write_xyz('test.xyz')
    raise

    '''正二十面体
    原子数:12
    对称性:'''
    cluster = Cluster()
    cluster.add_atom(Point(0, 1./2, (1 + sqrt(5)) / 4), elem)
    cluster = T(cluster)
    cluster.zoom(bond(elem,elem))

    c = Cluster()
    c.read('3D Atomistic (2).car')
    c.move(Point(0.600404,-0.641730))
    cc = deepcopy(c)
    for p in c.coordinate:
        pp = deepcopy(p)
        pp.x = p.x*cos(pi*0.4) - p.y*sin(pi*0.4)
        pp.y = p.x*sin(pi*0.4) + p.y*cos(pi*0.4)
        cc.add_atom(pp,'B')
        pp = deepcopy(p)
        pp.x = p.x*cos(pi*0.8) - p.y*sin(pi*0.8)
        pp.y = p.x*sin(pi*0.8) + p.y*cos(pi*0.8)
        cc.add_atom(pp,'B')
        pp = deepcopy(p)
        pp.x = p.x*cos(pi*1.2) - p.y*sin(pi*1.2)
        pp.y = p.x*sin(pi*1.2) + p.y*cos(pi*1.2)
        cc.add_atom(pp,'B')
        pp = deepcopy(p)
        pp.x = p.x*cos(pi*1.6) - p.y*sin(pi*1.6)
        pp.y = p.x*sin(pi*1.6) + p.y*cos(pi*1.6)
        cc.add_atom(pp,'B')
    unique(cc)
    cc.write_xyz('B70_C5.xyz')
    raise
    c=archimedes466('B')
    #c.zoom(bond('B','P')/bond('H','H'))
    c.write_xyz('B24P24.xyz')
    raise
    c = triangle_sheet6(8,6,6,4,'B')
    bending_boron(c)
    print(c.get_size())
    c.write_xyz('B_planar%d_9771.xyz'%c.get_size())
    raise

    c = archimedes3434()
    c.zoom(2)
    #c.zoom(2.77)
    c.set_element('Cl')
    print(c.get_size())
    c.write_xyz('3454.xyz')
    c.display()
