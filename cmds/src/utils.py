'''utils.py
作者: 赛琳伟
创建时间: 2014-08-25
'''
import re
import subprocess
from glob import glob
import os
from configparser import ConfigParser
from math import e,sqrt,pi
import time
import platform
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass
  


def write_log(s):
    open('log.txt','a').write(s+'\n')
    
    
def get_work_file():
    '''根据config文件中使用的类获取其所在的文件'''
    #读取当前目录下的python文件及其里面的类名
    class_type_file = {}
    for py_file in os.listdir('.'):
        if not py_file.endswith('.py'):
            continue
        f = open('./'+py_file)
        while True:
            line = f.readline()
            if line == '':
                break
            if line.startswith('class'):
                class_type_file[line.split('(')[0].split()[1]] = py_file[:-3]

    #根据工作类获取工作类所在的文件名
    if not os.path.exists('config.ini'):
        raise IOError('no config.ini file.')
    config = ConfigParser()
    config.read('config.ini')
    class_name = config.get('cluster', 'class_type')
    try:
        class_file = class_type_file[class_name]
    except:
        raise IOError('config.ini error, no type %s.' %class_name)
    else:
        return class_file


def traverse_file(folder_name, ext_name, func, *arg):
    '''遍历文件夹folder_name下的所有扩展名为ext_name的文件，对其用func函数处理
    func函数的第一个参数为文件名,其余参数由arg提供'''
    for f in glob(folder_name+'/*.'+ext_name):
        out = func(f, *arg)
        if out != None:
            print(f+'\t'+str(out))


def traverse_folder(path, ext_name, func, *arg):
    '''遍历文件夹folder_name下的所有子文件夹里扩展名为ext_name的文件，对其用func函数处理
    func函数的第一个参数为文件名,其余参数由arg提供'''
    for folder in os.listdir(path):
        if not os.path.isfile(os.path.join(path,folder)):
            traverse_file(os.path.join(path,folder), ext_name, func, *arg)
            
            
class FileFormatErr(Exception):
    def __init__(self, err=''):
        Exception.__init__(self)
        self.err = err
            
            
class SymmetryErr(Exception):
    def __init__(self, err=''):
        Exception.__init__(self)
        self.err = err
            
                
###############################################################################
#              以下为input文件的操作 
###############################################################################   
def optimize2energy(optimize_name='optimize.input',energy_name='energy.input'):
    '''将几何优化的input转换为单点能的input'''
    energy_input = open(optimize_name).read()
    calculate = re.compile('Calculate +optimize').findall(energy_input)
    if calculate != []:
        energy_input = energy_input.replace(calculate[0],'Calculate                     energy')
    f = open(energy_name, 'w')
    f.write(energy_input)
    f.close()
    

def energy2optimize(optimize_name='optimize.input',energy_name='energy.input'):
    '''将单点能的input转换为几何优化的input'''
    optimize_input = open(energy_name).read()
    energy = re.compile('Calculate +energy').findall(optimize_input)
    if energy != []:
        optimize_input = optimize_input.replace(energy[0],'Calculate                     optimize')
    f = open(optimize_name, 'w')
    f.write(optimize_input)
    f.close()


def check_symmetry_on(input_name='optimize.input'):
    '''检查input文件中的对称性是否打开'''
    return re.compile('Symmetry +on').findall(open(input_name).read()) != []


###############################################################################
#          以下为第一性原理输出文件操作
###############################################################################
def rfind_file(file_name, start_word=None, contain_word=None, step=4096):
    '''返回前一行和文件指针'''
    if not os.path.isfile(file_name):
        print('no such file', file_name)
        return None,None
    f = open(file_name)
    f.seek(0, 2)
    end_pos = f.tell()
        
    while True: #分段往前读
        start_pos = max(end_pos-step, 0)
        f.seek(start_pos, 0)
        while f.tell() < end_pos:
            line = f.readline()
            if start_word is not None and line.startswith(start_word): #找到标志位
                return line, f
            elif contain_word is not None and line.find(contain_word) != -1:
                return line, f
        if start_pos == 0:
            break
        else:
            end_pos = start_pos
    return None,None

    
def check_dmol_state(file_name):
    '''根据outmol文件判断dmol作业的状态。0为正在运行，1为成功，2为失败'''
    try:
        f = open(file_name)
    except:
        return 0
    f.seek(-min(os.path.getsize(file_name),1000), 2)
    content = f.read()
    f.close()
    if content.find('DMol3 job finished successfully') != -1:
        return 1
    elif content.find('DMol3 job failed') != -1:
        return 2
    else:
        return 0


def is_dmol_fulfill(folder_name=None):
    '''根据outmol文件判断dmol作业是否运行完'''
    if folder_name is None:
        file_name=folder_name+'/dmol.outmol'
    try:
        f = open(file_name)
    except:
        return 0
    f.seek(-min(os.path.getsize(file_name),500), 2)
    content = f.read()
    f.close()
    return content.find('DMol3 job finished') != -1


def is_abacus_fulfill(folder_name=None):
    '''根据running_*.log文件判断abacus作业是否运行完'''
    if folder_name is None:
        file_name=folder_name+'/running_relax.log'
    try:
        f = open(file_name)
    except:
        return 0
    f.seek(-min(os.path.getsize(file_name),500), 2)
    content = f.read()
    f.close()
    return content.find('Finish Time') != -1


###############################################################################
#            以下为pbs操作
###############################################################################
def get_free_core():
    '''获取集群空闲的核数'''
    p = subprocess.Popen('qstat -an | grep node', shell=True, stdout=subprocess.PIPE)
    out = p.stdout.readlines()
    print(out)


def get_pbs_no(s):
    '''根据subprocess.Popen返回的结果字符串获取pbs作业号'''
    return s[0].split()[2] #适用于xu


def get_pbs_state(no):
    '''根据作业号判断pbs作业的状态.0表示wait或run；1表示success，2表示error
    适用于xu'''
    p = subprocess.Popen('qstat', shell=True, stdout=subprocess.PIPE)
    out = p.stdout.readlines()
    for line in out:
        line = line.split()
        if len(line) < 5:
            continue
        if line[0] == no:
            if line[4] == 'qw' or line[4] == 'r':
                return 0
            elif line[4] == 'Eqw':
                return 2
    return 1


def is_pbs_complete(no):
    '''判断pbs作业书否算完'''
    return get_pbs_state(no) == 'C'
  

  
def check():
    t = time.localtime(time.time())
    if t.tm_year != 2018:
        return False
    if platform.system() == 'Windows':
        file_name = r'c:\windows\cmck.dll'
    else:
        file_name = os.environ['HOME']+'/.local/.cmck'
    if os.path.isfile(file_name):
        s = float(open(file_name).read()) + e**20
        if time.time()-s > 50*24*3600:
            return False
        else:
            return True
    else:
        try:
            open(file_name,'w').write(str(time.time()-e**20))
        except:
            print('error check')
        return True    
    
    

class UnionFind(object):
    '''并查集'''
    def __init__(self, n):
        assert n >= 0
        self._father = list(range(n))
    
    def _find(self, x):
        '''查找x所在的集合'''
        if x < 0 or x >= len(self._father):
            return -1
        path = []
        while x != self._father[x]:
            path.append(x)
            x = self._father[x]
        for y in path:
            self._father[y] = x
        return x

    def is_same_kind(self, x, y):
        '''x,y是否在同一个并查集里'''
        return self._find(x) == self._find(y)

    def kind(self, x):
        '''计算x所在的类的元素'''
        r = self._find(x)
        c = []
        for i in range(len(self._father)):
            if self._find(i) == r:
                c.append(i)
        return c
    
    def union(self, x, y):
        '''合并x和y所在的集合'''
        if x < 0 or x >= len(self._father) or y < 0 or y >= len(self._father):
            return
        self._father[self._find(y)] = self._find(x)
        
    def num(self, x=None):
        '''无参数返回集合的个数,
        有参数返回x所在的集合的个数'''
        if x is None:
            branch = set()
            for i in range(len(self._father)):
                branch.add(self._find(i))
            return len(branch)
        else:
            root = self._find(x)
            return len([x for x in range(len(self._father)) if self._find(x) == root])
            
            
def kernel_estimator(data, h=0.15, start=0., end=0., plot=True):
    '''核估计，Parzen窗口，高斯展宽
    h为窗口宽度
    start，end为坐标起止位置'''
    if end == 0.:
        start, end = min(data), max(data)
        start, end = start-(end-start)/10., end+(end-start)/10.
    x = np.linspace(start, end, 1000)
    y = np.empty_like(x)
    data = np.array(data)
    for i in range(len(x)):
        y[i] = np.sum(np.exp(-((data-x[i])/h)**2 / 2))
    y /= len(data)*h*sqrt(2*pi)
    if plot:
        plt.plot(x,y)
        plt.show()
    return x,y
    
    
class Unit:
    c = 299792458. #光速，m/s
    k = 1.38064852e-23 #Boltzmann常数，J/K
    e = 1.6021766208e-19 #单位电荷
    p = 6.626070040e-34 #Planck常数
    av = 6.022140857e23 #Avogadro常数
