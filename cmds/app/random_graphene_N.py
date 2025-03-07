# -*- coding: utf-8 -*-
"""
石墨烯随机替换N程序

原则：不能出现N-N键，N原子必须在某个六元环上
每次随机选择一个位置，然后随机执行以下2种操作中的1种，直到2种N原子的个数都达到给定个数：
1.将某位置的C原子替换为N原子
2.将某位置的原子删除，并将其相邻原子替换为N原子

操作1的规则：
该位置和相邻的3个原子必须都为C
操作2的规则：
该位置必须是C或吡啶N，该位置的邻位和非空邻位的邻位不能有N，N原子必须在某个六元环上（即该位置距离为3的位置不能是空）。为了让空洞尽量大，以一定概率（程序中设定为0.7）要求该位置是吡啶N

操作1：replace_N
操作2：pyridinic吡啶
操作3：pyridic吡咯，暂未实现

主程序：CN3, CN3s
"""
import sys
sys.path.append('../src')
from point import Point
from atom import bond,set_bond
from crystal import Crystal
from math import sqrt
from random import randrange,shuffle,random
import numpy as np
from copy import deepcopy


set_bond('C', 'C', 1.421)

def array2graphene(a):
    '''将数组转成结构
    参数a：numpy类型的二维数组
    返回：该数组对应的结构
    数组a的含义：3C, <0N, 0空位, -3替换N，-2吡啶N，-1吡咯N
    '''
    h, w = a.shape
    c = Crystal()
    for i in xrange(h):
        for j in xrange(w):
            if a[i,j] == 0:
                continue
            elem = 'C' if a[i,j]>0 else 'N'
            if a[i,j] == -1: #pyrrole
                pass
            else:
                c.add_atom(Point(j*sqrt(3)/2, i*1.5+(i+j)%2*0.5, 10), elem)
    c.a, c.b, c.c = sqrt(3)/2*w, 1.5*h, 20/bond('C','C') #晶格常数
    c.zoom(bond('C','C'))
    return c
    
def ortho(x,y,h,w):
    '''邻位'''
    f = 1 if (x+y)%2 else -1
    return [(x,(y-1)%w), (x,(y+1)%w), ((x+f)%h,y)]
    
def meta(x,y,h,w):
    '''间位'''
    f = 1 if (x+y)%2 else -1
    return [(x,(y-2)%w), (x,(y+2)%w), ((x+f)%h,(y-1)%w), ((x+f)%h,(y+1)%w), ((x-f)%h,(y-1)%w), ((x-f)%h,(y+1)%w)]
    
def para(x,y,h,w):
    '''对位'''
    f = 1 if (x+y)%2 else -1
    return [((x+f)%h,(y-2)%w), ((x+f)%h,(y+2)%w), ((x-f)%h,y)]
    
def ortho_para(x,y,h,w): 
    '''邻位的对位和对位的邻位'''
    f = 1 if (x+y)%2 else -1
    return [((x-2*f)%h,y), ((x+f)%h,(y-3)%w), ((x+f)%h,(y+3)%w),
            ((x+2*f)%h,y), ((x-f)%h,(y-3)%w), ((x-f)%h,(y+3)%w)]

def replace_N(a, x, y):
    h, w = a.shape
    for p in ortho(x,y,h,w): #相邻3个原子必须都是C
        if a[p] <= 0:
            return False
    a[x,y] = -3
    return True

def pyridine(a, x, y):
    '''
    返回：0-3个原子吡啶原子'''
    h, w = a.shape
    if a[x,y] == 0 or a[x,y] == -1: #该位置必须是C或吡啶N
        return 0
    #邻位全为C或者有2个吡啶N
    neighbour = ortho(x,y,h,w) #邻位
    vacancy = None #邻位中的空缺位
    origin = [(x,y)] #六元环的情况，原始点为3个，包含邻接原子中的2个吡啶N位
    for p in neighbour:
        if a[p] == -3:
            return 0
        elif a[p] == -2:
            origin.append(p)
        elif a[p] == 0:
            vacancy = p
    if len(origin) == 2: #只能为1或3
        return 0
    if vacancy:
        neighbour.remove(vacancy) #剔除空缺位
    if len(origin) == 1:
        add_N = 0 if a[x,y]==3 else -1 #增加的吡啶N个数
    else: #==3，六元环情况
        neighbour = reduce(list.__add__, [ortho(p[0],p[1],h,w) for p in origin]) #这3个点的邻位点
        neighbour = filter(lambda p:a[p]!=0 and p not in origin, neighbour) #去除空位
        add_N = -2 #先减少2个吡啶N，再增加3个
    #非空邻位的邻位不能是N原子
    for p in neighbour:
        for q in ortho(p[0],p[1],h,w):
            if q not in origin and a[q] < 0:
                return 0
    #N原子必须在某个六元环上
    for p in origin:
        for q in ortho_para(p[0],p[1],h,w):
            if a[q] == 0:
                return 0
    #删除原始原子
    for p in origin:
        a[p] = 0
    for p in neighbour: #将邻接原子标记为吡啶N
        a[p] = -2
        add_N += 1
    return add_N

def pyrrole(a, x, y):
    pass
    
def CN3(h, w, n1, n2, n3, file_name):
    '''产生1个结构
    w：石墨烯片的宽。必须为偶数！
    h：石墨烯片的高。必须为偶数！
    n1：每个随机结构中单N的个数
    n2：每个随机结构中pyridinic吡啶的个数
    n3：每个随机结构中pyrrolic的个数'''
    #初始石墨烯网格。六元环被压成矩形
    # ■■☐☐■■☐☐
    # ☐■■☐☐■■☐
    # ■■☐☐■■☐☐
    # 每个原子左右都有连接，但上下只是交替有
    a = np.ones((h,w),'int') * 3
    
    total_pos = range(w*h)
    shuffle(total_pos) #将网格中所有点打乱顺利
    for p in total_pos*3: #遍历网格中所有点
        if n1+n2+n3 == 0:
            break
        x,y = p/w, p%w
        idx = randrange(n1+n2+n3) #随机决定执行哪种操作
        if idx < n1:
            if replace_N(a, x, y):
                n1 -= 1
        elif idx < n1+n2:
            if random() < 0.1: #为了让空位尽可能相连，以一定的概率要求该位置必须是吡啶N或邻接位置有2个吡啶N
                total_pos2 = range(w*h)
                shuffle(total_pos2)
                for p in total_pos2:
                    x,y = p/w, p%w
                    if a[x,y] == -2 or [a[p] for p in ortho(x,y,h,w)].count(-2) == 2:
                        break
            n2 -= pyridine(a, x, y)
            if n2 < 0:
                n2 = 0
        else:
            if pyrrole(a, x, y):
                n3 -= 1
    c = array2graphene(a)
    c.write_car(file_name)
    print n1,n2,n3 #打印未实现的个数
    return n1,n2,n3
           
def CN3s(h, w, n, n1, n2, n3):
    '''产生n个结构'''
    for i in range(n):
        CN3(w, h, n1, n2, n3, 'cn3_%d.car'%i)
       
       
if __name__ == '__main__':
    CN3(20, 60, 0, 400, 0, 'test.car')
    #CN3s(20, 60, 1000, 10, 200, 0)
