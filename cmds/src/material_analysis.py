'''material.py
作者: 赛琳伟
创建时间: 2017-07-20
'''
from . import utils
from .io import read
from .atom import bond as BOND
from math import pi
import os
try: #这些模块有些机器可能没装
    import matplotlib.pyplot as plt #display函数使用
    from mpl_toolkits.mplot3d import Axes3D #display函数使用
except:
    pass



###############################################################################
#        以下为第一性原理计算结果分析
###############################################################################
def photo_electron_spectrum(file_name, VDE, start, end):
    '''可以从outmol或OUTCAR或gauss的LOG文件中读取数据
    outmol文件为负离子团簇优化后的文件
    VDE为阴离子的能量-阴离子结构的中性能量
    start，end为坐标起止位置'''
    def read_outmol():
        step = 4096
        f = open(file_name)
        f.seek(0,2)
        size = f.tell()
        f.seek(max(size-step-1,0), 0)
        last_start = f.tell()

        while True: #分段往前读
            line = f.readline()
            if line == '':
                raise IOError('file: %s format error.' %file_name)
            if line.find('    state                         eigenvalue        occupation') != -1:
                break
            if f.tell() >= last_start + step:
                if f.tell() >= 2*step:
                    f.seek(-2*step, 1)
                    last_start = f.tell()
                else:
                    raise IOError('file: %s format error.' %file_name)

        f.readline(), f.readline()
        energyLevel = []
        while True:
            line = f.readline()
            if len(line) < 4:
                break
            line = line.split()
            if float(line[-1]) < 0.5:
                break
            energyLevel.append(float(line[-2]))
        return [VDE + energyLevel[-1] - e for e in energyLevel]

    def read_outcar():
        band_energy1 = []
        band_energy2 = []
        f = open('OUTCAR')
        f.seek(0,2)
        size = f.tell()
        f.seek(max(size-10000,0), 0)
        line = f.readline()
        while line != '':
            if line == ' spin component 1\n':
                f.readline(), f.readline(), f.readline()
                while True:
                    line = f.readline()
                    if line == '\n':
                        break
                    data = [float(s) for s in line.split()]
                    band_energy1.append(data[1])
                    if data[2] > 0.5:
                        be1 = data[1]
                        occupy1 = data[2]
            if line == ' spin component 2\n':
                f.readline(), f.readline(), f.readline()
                while True:
                    line = f.readline()
                    if line == '\n':
                        break
                    data = [float(s) for s in line.split()]
                    band_energy2.append(data[1])
                    if data[2] > 0.5:
                        be2 = data[1]
                        occupy2 = data[2]
            line = f.readline()

        be = be1 if occupy1<occupy2 else be2
        shift = VDE + be
        band_energy1 = [shift-e for e in band_energy1]
        band_energy2 = [shift-e for e in band_energy2]
        return band_energy1 + band_energy2

    def read_log():
        '''
         The electronic state is语句下面是
         Alpha  occ.
         Alpha virt.
          Beta  occ.
          Beta virt.
        '''
        _, f = utils.rfind_file(file_name,' The electronic state is')
        band_energy = []
        while True:
            line = f.readline().split()
            if line[1] == 'virt.':
                break
            for v in line[4:]:
                while True:
                    #数据之间没空格了
                    second_neg = v.rfind('-')
                    if second_neg == 0:
                        band_energy.append(float(v))
                        break
                    elif second_neg != -1:
                        band_energy.append(float(v[second_neg:]))
                        v = v[:second_neg]
                    else:
                        band_energy.append(float(v))
                        break
        occupy = max(band_energy)*27.2114
        shift = VDE + occupy
        band_energy = [shift-e*27.2114 for e in band_energy]
        return band_energy

    if file_name.endswith('outmol'):
        energy_level = read_outmol()
    elif file_name == 'OUTCAR':
        energy_level = read_outcar()
    elif file_name.endswith('LOG') or file_name.endswith('out'):
        energy_level = read_log()
    else:
        print('do not support this file type.')
        return
    x,y = utils.kernel_estimator(energy_level, 0.1, start, end)
    if len(file_name.split(os.path.sep)) == 1:
        write_name = file_name.split('.')[0]
    else:
        write_name = file_name.split(os.path.sep)[-2]
    f = open(write_name+'.txt','w')
    for xx,yy in zip(x,y):
        f.write('%f\t%f\n' %(xx,yy))

def gauss_dos(file_name, VDE, start, end):
    energy_level = read_log()
    x,y = utils.kernel_estimator(energy_level, 0.1, start, end)
    if file_name.split(os.path.sep) == 1:
        write_name = file_name.split('.')
    else:
        write_name = file_name.split(os.path.sep)[-2]
    f = open(write_name+'.txt','w')
    for xx,yy in zip(x,y):
        f.write('%f\t%f\n' %(xx,yy))

###############################################################################
#                           以下为团簇分析
###############################################################################
def bond_distribution(c):
    d = []
    for i in range(c.get_size()):
        for j in range(i):
            d.append(c.atoms[i].dist(c.atoms[j]))
    x,y = utils.kernel_estimator(d)
    f = open('bond_distribution.txt','w')
    for xx,yy in zip(x,y):
        f.write('%f\t%f\n' %(xx,yy))

def angle_distribution(c):
    a = []
    for i in range(c.get_size()):
        for j in range(c.get_size()):
            if c.atoms[i].bond(c.atoms[j]) and j!=i:
                for k in range(j):
                    if c.atoms[i].bond(c.atoms[k]) and k != i:
                        angle = (c.atoms[i]-c.atoms[j]).angle(c.atoms[i]-c.atoms[k]) * 180/pi
                        if angle < 0:
                            angle = 180-angle
                        a.append(angle)
    x,y = utils.kernel_estimator(a, 5, start=0, end=180)
    f = open('angle_distribution.txt','w')
    for xx,yy in zip(x,y):
        f.write('%f\t%f\n' %(xx,yy))
    a.sort()
    open('angle.txt','w').write('\n'.join([str(_) for _ in a]))

def aver_bond_len(c, cutoff):
    '''平均键长
    c为要计算的团簇
    单质时cutoff为键长的截断；化合物为容忍度
    可通过bond_distribution函数分析出cutoff'''
    if len(c.get_elements_count()) == 1:
        cutoff /= BOND(c.atoms[0].elem, c.atoms[0].elem)

    bonds = {}
    for i in range(c.get_size()):
        for j in range(i):
            d = c.atoms[i].dist(c.atoms[j])
            if d < BOND(c.atoms[i].elem, c.atoms[j].elem) * cutoff:
                if (c.atoms[i].elem, c.atoms[j].elem) in bonds:
                    bonds[(c.atoms[i].elem, c.atoms[j].elem)].append(d)
                elif (c.atoms[j].elem, c.atoms[i].elem) in bonds:
                    bonds[(c.atoms[j].elem, c.atoms[i].elem)].append(d)
                else:
                    bonds[(c.atoms[i].elem, c.atoms[j].elem)] = [d]
    for key,value in list(bonds.items()):
        print('%s-%s\t%.3f'%(key[0], key[1], sum(value)/len(value)))

def aver_degree(c, cutoff):
    '''平均配位数
    c为要计算的团簇
    单质时cutoff为截断；化合物为容忍度'''
    elem_count = c.get_elements_count()
    if len(elem_count) == 1:
        cutoff /= BOND(c.atoms[0].elem, c.atoms[0].elem)

    coords = {}
    for i in range(c.get_size()):
        if c.atoms[i].elem not in coords:
            coords[c.atoms[i].elem] = 0.
        for j in range(c.get_size()):
            if j == i:
                continue
            if c.atoms[i].dist(c.atoms[j]) < BOND(c.atoms[i].elem, c.atoms[j].elem) * cutoff:
                coords[c.atom[i].elem] += 1.
    for key,value in list(coords.items()):
        print('%s\t%.2f'%(key, value/elem_count[key]))

def bond_order(file_name, cutoff):
    '''计算键级
    file_name: outmol文件名（dmol）或LOG文件名（Gauss）
    单质时cutoff为键长的截断；化合物为容忍度'''
    c = read(file_name)
    if len(c.get_elements_count()) == 1:
        cutoff /= BOND(c.atoms[0].elem, c.atoms[0].elem)
    if os.path.splitext(file_name)[1] == '.outmol':
        bond_orders = Dmol3.bond_order(file_name)
    elif os.path.splitext(file_name)[1] == '.LOG':
        bond_orders = Gaussian.bond_order(file_name)
    else:
        print('file format error.')
        raise

    bond_order_dict = {}
    for i in range(c.get_size()):
        for j in range(i):
            d = c.atoms[i].dist(c.atoms[j])
            if d < BOND(c.atoms[i].elem, c.atoms[j].elem) * cutoff:
                if (c.atoms[i].elem, c.atoms[j].elem) in bond_order_dict:
                    bond_order_dict[(c.atoms[i].elem, c.atoms[j].elem)].append(bond_orders[j,i])
                elif (c.atoms[j].elem, c.atoms[i].elem) in bond_order_dict:
                    bond_order_dict[(c.atoms[j].elem, c.atoms[i].elem)].append(bond_orders[j,i])
                else:
                    bond_order_dict[(c.atoms[i].elem, c.atoms[j].elem)] = [bond_orders[j,i]]
    for key,value in list(bond_order_dict.items()):
        print('%s-%s\t%.3f'%(key[0], key[1], sum(value)/len(value)))

###############################################################################
#                           以下为结构显示
###############################################################################
def display(c, adj=None):
    '''显示团簇.
    参数adj: 矩阵类型,团簇的邻接矩阵.若不指定,则在该函数中计算
    返回: 无'''
    if adj is None:
        adj = c.adjacent()
    fig = plt.figure(figsize=(8,8))
    ax = Axes3D(fig) #ax = fig.gca(projection='3d')
    ax.patch.set_facecolor('black') #背景色
    ax.set_axis_off() #不显示坐标轴

    for i in range(c.get_size()):
        ax.scatter(c.atoms[i].coo.x, c.atoms[i].coo.y, c.atoms[i].coo.z,
                   s=400, c=get_color(c.atoms[i].elem)) #画原子
        for j in range(i): #画键
            if adj[i,j] == 0:
                continue
            x = [c.atoms[i].coo.x, c.atoms[j].coo.x]
            y = [c.atoms[i].coo.y, c.atoms[j].coo.y]
            z = [c.atoms[i].coo.z, c.atoms[j].coo.z]
            ax.plot(x, y, z, lw=4, c='#7f7f7f')

    #设置坐标轴显示范围
    left_bound = 0.8*min(min([p.coo.x for p in c.atoms]),
                     min([p.coo.y for p in c.atoms]),
                     min([p.coo.z for p in c.atoms]))
    right_bound = 0.8*max(max([p.coo.x for p in c.atoms]),
                     max([p.coo.y for p in c.atoms]),
                     max([p.coo.z for p in c.atoms]))
    ax.set_xlim(left_bound, right_bound)
    ax.set_ylim(left_bound, right_bound)
    ax.set_zlim(left_bound, right_bound)

    plt.show()

def write_png(c,file_name):
    '''保存为图片'''
    adj = c.adjacent()
    fig = plt.figure(figsize=(8,8))
    ax = Axes3D(fig) #ax = fig.gca(projection='3d')
    ax.patch.set_facecolor('black') #背景色
    ax.set_axis_off() #不显示坐标轴

    for i in range(c.get_size()):
        ax.scatter(c.atoms[i].x, c.atoms[i].y, c.atoms[i].z,
                   s=400, c=get_color(c.atoms[i].elem)) #画原子
        for j in range(i): #画键
            if adj[i,j] == 0:
                continue
            x = [c.atoms[i].coo.x, c.atoms[j].coo.x]
            y = [c.atoms[i].coo.y, c.atoms[j].coo.y]
            z = [c.atoms[i].coo.z, c.atoms[j].coo.z]
            ax.plot(x, y, z, lw=4, c='#7f7f7f')

    #设置坐标轴显示范围
    left_bound = 0.8*min(min([p.coo.x for p in c.atoms]),
                     min([p.coo.y for p in c.atoms]),
                     min([p.coo.z for p in c.atoms]))
    right_bound = 0.8*max(max([p.coo.x for p in c.atoms]),
                     max([p.coo.y for p in c.atoms]),
                     max([p.coo.z for p in c.atoms]))
    ax.set_xlim(left_bound, right_bound)
    ax.set_ylim(left_bound, right_bound)
    ax.set_zlim(left_bound, right_bound)

    plt.savefig(file_name)

if __name__=="__main__":
    c = read('B96.car')
    angle_distribution(c)
    raise
    photo_electron_spectrum(r'E:\cluster\B\gauss_opt\B20-\anion_opt\B20-_03/temp.out',3.18,2,6)
    raise
    ga = GA_Detail()
    ga.read_record(r'E:\cluster\B\search\43\B43_C3')
    pop = ga.top_best(16)
    f = open('recover_test.txt','w')
    for i,c in enumerate(pop):
        f.write('pop %d\n'%i+str(c)+'\n')
