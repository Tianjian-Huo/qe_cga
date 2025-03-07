'''ga.py
遗传算法类
作者: 赛琳伟
创建时间: 2014-05-04
修改记录:
'''
from configparser import ConfigParser
import os
import sys
import time
from copy import deepcopy
from random import random
import numpy as np
from ..src import utils
from ..src import fix_cluster
from ..src.io import write_car,read_cga_obj
from ..src.optimize import optimize
            
            
class CGA:
    '''遗传算法类'''

    def __init__(self):
        '''初始化'''
        #读取配置文件
        self.config()

        #初始化成员变量
        self.pop = [] #种群
        self.gen = 0 #代数
        self.father = None
        self.mother = None
        self.breed_type = None
        self.son = self.cluster_type()
        self.father = self.cluster_type()
        self.mother = self.cluster_type()
        self.iso_pair = None #同构的结构对

    def config(self):
        '''读取配置文件'''
        config = ConfigParser()
        if os.path.isfile('config.ini'):
            config.read('config.ini')
        else:
            #config.read(os.path.join(os.path.split(__file__)[0],'config.ini'))
            print('请提供config.ini文件')
            sys.exit(0)

        #团簇类型配置
        self.cluster_type = fix_cluster.config()
        #读取遗传算法参数: 种群大小,迭代次数
        self.pop_size = config['genetic algorithm'].getint('pop_size')
        self.max_gen = config['genetic algorithm'].getint('max_gen')
        #读取产生后代方式及概率
        self.method_str = config['genetic algorithm'].get('breed_method').split(',')
        self.method_prob = config['genetic algorithm'].get('breed_probability').split(',')
        assert len(self.method_str) == len(self.method_prob)
        self.method = []
        for i in range(len(self.method_str)):
            self.method_prob[i] = float(self.method_prob[i])
            self.method_str[i] = self.method_str[i].strip()
            if not hasattr(self.cluster_type, self.method_str[i]):
                print('no function %s in' %self.method_str[i], self.cluster_type)
            self.method.append(getattr(self.cluster_type, self.method_str[i]))
        assert 0.9 < sum(self.method_prob) < 1.1
        #转换概率
        sum_ = sum(self.method_prob)
        self.method_prob = [p/sum_ for p in self.method_prob]
        if 'breed_prob_end' in config['genetic algorithm']:
            prob_end = config['genetic algorithm'].get('breed_prob_end').split(',')
            assert len(self.method_str) == len(prob_end)
            prob_end = [float(p) for p in prob_end]
            prob_end = [p/sum(prob_end) for p in prob_end]
            self.method_prob_increment = [(end-start)/self.max_gen for start,end in zip(self.method_prob,prob_end)]
        else:
            self.method_prob_increment = None
            
    def read_recover(self):
        '''从recover中读取种群'''
        f = open('recover.txt')
        for i in range(self.pop_size): #读结构
            f.readline()
            c = self.cluster_type()
            stru, _ = read_cga_obj(f)
            c.deepcopy(stru)
            if not c.check():
                print('Please check pop %d in recover.txt.' %i)
                os._exit(1)
            self.pop.append(c)

    @staticmethod
    def get_cluster_info(c):
        info ='energy=%-.7f, fp=' %c.get_energy()
        for fp in c.get_fingerprint():
            info += '%.2f '%fp
        return info

    def init_pop(self):
        '''产生初始种群'''
        #打印CPU使用信息
        try:
            from multiprocessing import cpu_count
            print('total CPU number:',cpu_count())
            import psutil
            print('CPU use: %.1f%%' %psutil.cpu_percent(1))
        except:
            pass

        if os.path.exists('recover.txt'): #接上次计算
            print('continue genetic algorithm study.')
            utils.write_log('task continue at: '+time.ctime())
            #从record文件里读代数
            if os.path.isfile('record.txt'):
                f = open('record.txt')
                f.seek(max(os.path.getsize('record.txt')-25000,0), 0)
                while True:
                    line = f.readline()
                    if line == '':
                        break
                    if line.startswith('gen'):
                        self.gen = int(line.split()[1])
            else:
                self.gen = 0
            self.read_recover()
        else: #新计算
            print('start new genetic algorithm study.')
            print('initial population:')
            utils.write_log('task start at: '+time.ctime())
            if os.path.isfile('energy.txt'): #删除之前的energy文件
                os.remove('energy.txt')
            if os.path.isfile('record.txt'): #删除之前的record文件
                os.remove('record.txt')
            if os.path.isfile('log.txt'): #删除之前的log文件
                os.remove('log.txt')

            while len(self.pop) < self.pop_size:
                utils.write_log('init %d'%len(self.pop))
                c = self.cluster_type()
                c.random()
                optimize(c)
                print('%14d, '%len(self.pop) + self.get_cluster_info(c))
                write_car(c, '%02d_0000.car' %len(self.pop)) #保存
                open('record.txt', 'a').write('\ninit%3d\n' %len(self.pop) + str(c)) #写record文件
                self.pop.append(c)
            self._save() #写recover文件

        self._write_energy()

    def _save(self):
        '写恢复文件'
        f = open('recover.txt', 'w')
        for i in range(self.pop_size):
            f.write('pop%4d\n'%i)
            f.write(str(self.pop[i]))
            f.write('\n')
        f.close()

    @staticmethod
    def roulette(p):
        '''轮盘赌选择。要求sum(p)==1'''
        needle = random()
        cp = 0
        for i in range(len(p)):
            cp += p[i]
            if needle <= cp:
                return i

    def _select(self):
        '''轮盘赌方式选择父代'''
        e = [c.get_energy() for c in self.pop]
        worst = max(e)
        diff = np.array([worst-x for x in e])
        p = diff + np.sum(diff)/self.pop_size/5 #最差的被选择的概率为平均概率的1/5
        p /= p.sum()
        #roulette = stats.rv_discrete(values=(list(range(self.pop_size)),p))
        #idx = roulette.rvs()
        idx = self.roulette(p)
        parent = deepcopy(self.pop[idx])
        parent_info = '    parent:%3d, '%idx + self.get_cluster_info(parent)
        print(parent_info)
        record = open('record.txt', 'a')
        record.write(parent_info+'\n')
        return parent

    def _breed(self):
        '''产生子代'''
        print('%d:' %self.gen, end=' ')
        record = open('record.txt', 'a')
        record.write('\ngen%7d\n' %self.gen)
        i = self.roulette(self.method_prob) #选择产生子代方式
        self.breed_type = self.method_str[i] #记录产生子代方式
        print('    method=',self.breed_type)
        record.write('method=%s\n' %self.method_str[i])
        record.flush()
        self.father = self._select()
        if self.method[i].__code__.co_argcount == 3: #杂交
            self.mother = self._select()
            self.method[i](self.son, self.father, self.mother) #产生子代
        else: #变异
            self.son = deepcopy(self.father)
            self.method[i](self.son) #产生子代
        #记录子代
        record.write('son\n')
        record.write(str(self.son))
        record.close()

    def _relax(self):
        '''优化结构'''
        try:
            optimize(self.son)
        except utils.SymmetryErr: #dmol可能将对称结构优化为没有对称性!!!
            self.son = deepcopy(self.pop[0])
            self.son.set_energy(0.)
        finally:
            pass
        print('       son:   , ' + self.get_cluster_info(self.son))
        open('record.txt', 'a').write('optimized\n' + str(self.son)) #写记录文件

    def _get_replaced(self):
        replaced, reason = None, None #被替换的个体编号,原因
        for i in range(self.pop_size):
            if self.son.isomorphic(self.pop[i]): #同构
                if self.son.get_energy() < self.pop[i].get_energy(): #同构并且能量更小,替换
                    replaced = i #替换
                    reason = 'isomorphic'
                else:
                    reason = 'isomorphic with %d' %i
                break
        else: #不同构
            #检查种群中是否有重复的
            try:
                for i in range(self.pop_size): #则替换种群中同构的.效果未知,暂用
                    for j in range(i):
                        if self.pop[j].isomorphic(self.pop[i]):
                            raise
            except: #有重复的
                if self.pop[i].get_energy() < self.pop[j].get_energy(): #i,j中能量高的被替换
                    i,j = j,i
                replaced = i #替换
                reason = 'repeated'
            else: #种群中无同构的,则替换最差的
                worst = np.argmax([c.get_energy() for c in self.pop])
                if self.son.get_energy() < self.pop[worst].get_energy():
                    replaced = worst #替换
                reason = 'worst'
        return replaced, reason

    def _replace(self, replaced, reason):
        '''替换种群'''
        if replaced is None:
            return
        record = open('record.txt', 'a')
        #记录信息
        if replaced == None:
            print('    not replace, because',reason)
            record.write('not replace, because %s\n' %reason)
            record.close()
            return None
        utils.write_log('    replace ' + reason + ' %d'%replaced)
        print('    replace ' + reason + ' %d'%replaced)
        record.write('replace ' + reason + ' %d\n'%replaced)
        record.close()
        #替换种群
        self.pop[replaced] = deepcopy(self.son)
        #更新结构文件
        write_car(self.son, '%02d_%04d.car' %(replaced,self.gen))
        #写energy文件
        self._write_energy()
        #保存状态
        self._save()

    def _write_energy(self):
        '''更新energy.txt文件'''
        worst = np.argmax([c.get_energy() for c in self.pop])
        best = np.argmin([c.get_energy() for c in self.pop])
        #更新energy文件
        f = open('energy.txt', 'a')
        #(代数,最低能量,最高能量,最低能量个体)
        f.write('%5d:%3d%16.7f%16.7f %s\n' %(self.gen, best,
                                             self.pop[best].get_energy(),
                                             self.pop[worst].get_energy(),
                                             self.breed_type))
        f.close()

    def iterator_one(self):
        '''迭代一次'''
        utils.write_log('step %6d'%self.gen)
        self._breed() #产生子代
        self._relax() #优化子代
        replaced, reason = self._get_replaced()
        self._replace(replaced, reason) #替换种群
        if self.method_prob_increment is not None:
            self.method_prob = [p+inc for p,inc in zip(self.method_prob,self.method_prob_increment)]

    def __call__(self):
        '''运行遗传算法'''
        self.init_pop()
        for self.gen in range(self.gen+1, self.max_gen):
            self.iterator_one()


class CGA_Life(CGA):

    def init_pop(self):
        super(CGA_Life, self).init_pop()
        e = np.array([x.get_energy() for x in self.pop])
        idx = e.argsort()
        for i in range(self.pop_size):
            self.pop[i].life = self.pop_size - idx[i]

    def _breed(self):                          
        super()._breed()
        self.son.life = 6

    def iterator_one(self):
        super().iterator_one()
        #替换生命结束的
        for i,p in enumerate(self.pop):
            p.life -= 1
            if p.life == 0:
                self._breed()
                self._relax()               
                self._replace(i, 'died')

if __name__ == '__main__':
    '''主函数'''
    #if not utils.check():
    #    os._exit(0)
    ga = CGA()
    ga()
    #os.system('pause')

