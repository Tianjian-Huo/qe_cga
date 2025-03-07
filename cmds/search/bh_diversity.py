'''basin_hopping.py
basin hopping类
作者: 赛琳伟
创建时间: 2020-09-04
'''
import sys
sys.path.append('../..')
from math import exp
from random import random
from configparser import ConfigParser
import os
from copy import deepcopy
from atom import get_element_id
from abinitio import abinitio
from fix_cluster import CLUSTER



class BasinHopping_Diversity(object):
    '''basin hopping类，考虑结构之间的不同构，为了产生更多有效结构'''
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
        self.step = 0.2 #每次微扰的步长
        self.iter = 0 #当前迭代次数
        self.accept_num = 0 #interval内接受的次数
        self.history_stru = []

        #读取配置文件
        self.config()

        #初始结构
        self.x = CLUSTER()
        if os.path.isfile('recover.txt'): #接上次计算
            print('continue basin hopping...')
            print('last iter:')
            f = open('recover.txt')
            line = f.readline().split()
            self.iter = int(line[0])
            self.accept_num = int(line[1])
            self.best = CLUSTER()
            self.best.set_energy(float(line[2]))
            self.x.read_file_obj(f)
        else: #新计算
            self.x.random()
            abinitio(self.x) #计算初始结构的能量
            self.x.center()
            self.x.calc_fingerprint()
            self.best = deepcopy(self.x)
            
            f = open('record.txt','w')
            f.write('iter=0\n')
            f.write(str(self.x))
            f.write('\n')
            
            open('recover.txt','w').write('%d %d %f\n'%(self.iter,self.accept_num, \
                self.best.get_energy()) + str(self.x)) #写恢复文件
            
            arc = open('basin_hopping.arc', 'w')
            arc.write('!BIOSYM archive 3\nPBC=OFF\n') #写arc文件头
            self.x.append_arc('basin_hopping.arc') #将初始结构写入arc文件
            
            if os.path.isfile('log.txt'): #删除之前的log文件
                os.remove('log.txt')
        print('iter=%-5d, energy=%-15.7f' %(self.iter, self.x.get_energy()))
        self.history_stru.append(self.x)
            
    def accept(self):
        print('    accept')
        self.record.write('    accept\n')
        self.accept_num += 1 #记录接受次数
        self.x = deepcopy(self.next)
        self.x.append_arc('basin_hopping.arc')
        open('recover.txt','w').write('%d %d %f\n'%(self.iter,self.accept_num, \
            self.best.get_energy()) + str(self.x)) #写恢复文件
        if self.next.get_energy() < self.best.get_energy(): #更新最优结构
            self.best = deepcopy(self.next)
            print('new best')
            self.record.write('    new best\n')
            self.best.write_xyz('best.xyz')
            
    def adjust_step(self):
        accept_rate = float(self.accept_num) / self.interval
        if accept_rate > self.target_accept_rate:
            self.step /= self.factor #接受频率大,可能陷入了一个basin,增大step跳出
        else:
            self.step *= self.factor #接受频率小,减小step
        self.accept_num = 0
        print('step=%6.4f,  accept rate: %6.2f' %(self.step, accept_rate))
        self.record.write('step=%6.4f,  accept rate: %6.2f\n' %(self.step, accept_rate))
    
    def run(self):
        '''执行basin hopping优化'''        
        for self.iter in range(self.iter+1, self.max_gen):
            #在原结构x上随机微扰一下,得到新结构next
            self.next = deepcopy(self.x)
            self.next.perturbation(self.step)
            for c in self.history_stru: #不得跟以往同构
                if self.next.isomorphic(c):
                    break
            else:
                self.accept_num += 1 #同构太多或接受太多都需要增大setp
                continue
            
            #优化新结构
            abinitio(self.next)
            self.next.center()
            self.next.calc_fingerprint()
            self.history_stru.append(self.next)
            e = self.next.get_energy()
            print('iter=%-5d, energy=%15.7f' %(self.iter,e))
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
                
            #调整step
            if self.iter % self.interval == 0:
                self.adjust_step()
            
            self.record.write('\n')
            self.record.close()
        
    def __call__(self):
        '''以对象()的形式调用
        返回最终的最优结构'''
        self.run()
        return self.best


if __name__ == '__main__':
    '''主函数'''
    bs = BasinHopping_Diversity()
    bs()
