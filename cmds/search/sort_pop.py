import sys
sys.path.append('../..')
from cmds.src.io import read_cga_obj
import os
import numpy as np


def sort_pop(folder_name, prefix=''):
    '''种群排序
    给定recover.txt，对里面的种群按能量从低到高排序。在屏幕上打印出每个结构的能量，并将结构写入文件。
    文件名形如prefix_序号.car
    '''
    pop = []
    f = open(folder_name+'/recover.txt')
    size = os.path.getsize(folder_name+'/recover.txt')

    while f.tell() < size:
        f.readline()
        _, c = read_cga_obj(f)
        pop.append(c)

    for i in range(len(pop)):
        idx = np.argmin([c.get_energy() for c in pop])
        pop[idx].write_car(prefix+'_%02d.car' %i)
        print(pop[idx].get_energy())
        del pop[idx]


if __name__ == '__main__':
    sort_pop(foler_name=r'E:\cluster\alloy\CuRh\Cu13Rh2-ok', prefix='Cu13Rh2')
