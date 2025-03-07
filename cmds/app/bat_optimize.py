import sys
sys.path.append('../..')
from glob import glob
import os
from cmds.src.material_analysis import read
from cmds.src.optimize import *


def bat_optimize(directory, method):
    '''directory：结构所在文件夹，并将optimize.input放到该文件夹！
    method：第一性原理计算法方法,dmol3,vasp,gaussian等
    当前计算的结构会被删除，计算结束前放入结构也会被计算'''
    if not os.path.isdir(directory):
        print('no such directory:',directory)
        sys.exit()
    else:
        os.chdir(directory)
        if not os.path.isfile('optimize.input'):
            print('Please put optimize.input in',directory)
            sys.exit()
            
    set_optimize_method(method)
        
    loop = True
    while loop:
        loop = False
        for file_name in glob('*.car') + \
            glob('*.xyz') + \
            glob('*.xsd'):
            print(file_name,end='')
            loop = True
            c = read(f)
            os.remove(file_name)
            file_name, _ = os.path.splitext(file_name)
            optimize(c, folder_name=file_name, task_name=file_name)
            print('\t%f' %c.get_energy())
                        
                        
if __name__ == '__main__':
    bat_optimize('.', 'dmol3')
