'''basin_hopping.py
basin hopping类
作者: 赛琳伟
创建时间: 2014-08-19
'''
import sys
sys.path.append('../..')
from math import exp
from random import random
from configparser import ConfigParser
import os
from copy import deepcopy
from cmds.src.atom import get_element_id
from cmds.src.optimize import optimize
from cmds.src import fix_cluster
from cmds.src.io import read_cga_obj, write_car, append_arc



class BasinHopping(object):
    '''basin hopping类
    调用方式:
    bs = BasinHopping() #参数请参阅__init__函数
    bs()
    中间接受的结构保存在arc文件中'''
    def config(self):
        cfg = ConfigParser()
        cfg.read('config.ini')
        self.max_gen = cfg.getint('genetic algorithm', 'max_gen')
    
    def __init__(self, T=300, disp=True):
        '''构造函数
        参数structure: 团簇类型,要优化的结构
        参数arc_name: 要保存的arc文件名
        参数T: 温度. 300K时能量高0.00778时,以50%的概率接受
        参数disp: 是否打印信息'''
        self.interval = 50 #每多少步调整一次步长
        self.target_accept_rate = 0.5 #以该接受率为界调整步长.大于则增大,小于则减小
        self.factor = 0.9 #每次步长调整为原来的factor或1/factor倍 
        self._beta = 11604.8 / T
        self.step = 0.5 #每次微扰的步长
        self.iter = 0 #当前迭代次数
        self.accept_num = 0 #interval内接受的次数

        #读取配置文件
        self.config()
        self.cluster_type = fix_cluster.config()

        #初始结构
        self.x = self.cluster_type()
        if os.path.isfile('recover.xyz'): #接上次计算
            print('continue basin hopping...')
            f = open('recover.xyz')
            f.readline()
            stru, line = read_cga_obj(f)
            line = line.split()
            self.iter = int(line[0])
            self.accept_num = int(line[1])
            self.step = float(line[2])
            self.best = self.cluster_type()
            self.best.set_energy(float(line[3]))
            self.x.deepcopy(stru)
        else: #新计算
            self.x.random()
            optimize(self.x) #计算初始结构的能量
            self.best = deepcopy(self.x)
            if os.path.isfile('basin_hopping.arc'):
                os.remove('basin_hopping.arc')
            
            f = open('record.txt','w')
            f.write('iter=0\n')
            f.write(str(self.x))
            f.write('\n')
            
            open('recover.xyz','w').write('%d\n'%self.x.get_size() + str(self.x) + \
                '%d %d %d %f\n'%(self.iter, self.accept_num, self.step, \
                                 self.best.get_energy())) #写恢复文件
            
            append_arc(self.x, 'basin_hopping.arc') #将结构写入arc文件
            
            if os.path.isfile('log.txt'): #删除之前的log文件
                os.remove('log.txt')
        print('iter=%-5d, energy=%-15.7f' %(self.iter, self.x.get_energy()))
            
    def accept(self):
        print(', accept', end='')
        self.record.write('    accept\n')
        self.accept_num += 1 #记录接受次数
        self.x = deepcopy(self.next)
        append_arc(self.x, 'basin_hopping.arc')
        open('recover.xyz','w').write('%d\n'%self.x.get_size() + str(self.x) + \
            '%d %d %d %f\n'%(self.iter, self.accept_num, self.step, \
                             self.best.get_energy())) #写恢复文件
        if self.next.get_energy() < self.best.get_energy(): #更新最优结构
            self.best = deepcopy(self.next)
            print(', new best', end='')
            self.record.write('    new best\n')
            write_car(self.best, 'best.car')
            
    def adjust_step(self):
        accept_rate = float(self.accept_num) / self.interval
        if accept_rate > self.target_accept_rate:
            self.step /= self.factor #接受频率大,可能陷入了一个basin,增大step跳出
        else:
            self.step *= self.factor #接受频率小,减小step
        self.accept_num = 0
        print('step=%6.4f,  accept rate: %6.2f' %(self.step, accept_rate))
        self.record.write('step=%6.4f,  accept rate: %6.2f\n' %(self.step, accept_rate))
        
    def __call__(self):
        '''以对象()的形式调用
        返回最终的最优结构'''       
        for self.iter in range(self.iter+1, self.max_gen):
            #在原结构x上随机微扰一下,得到新结构next
            self.next = deepcopy(self.x)
            self.next.perturbation(self.step)
            
            #优化新结构
            optimize(self.next)
            e = self.next.get_energy()
            print('iter=%-5d, energy=%-15.7f' %(self.iter,e), end='')
            self.record = open('record.txt', 'a')
            self.record.write('iter=%d\n' %self.iter)
            self.record.write(str(self.next))
            
            #以一定概率接受新结构
            try:
                ratio = exp(-(e - self.x.get_energy()) * self._beta)
            except:
                ratio = 0. if e > self.x.get_energy() else 1.
            if ratio > random(): #接受结构
                self.accept()
            print()
                
            #调整step
            if self.iter % self.interval == 0:
                self.adjust_step()
            
            self.record.write('\n')
            self.record.close()
        return self.best


if __name__ == '__main__':
    '''主函数'''
    bs = BasinHopping()
    bs()
