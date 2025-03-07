import sys
sys.path.append('../..')
import os
from shutil import copy
from cmds import *


def extract_result(directory='.'):
    '''这里为文件夹所在的文件夹路径
    如果路径中包含中文，应该将路径前面的r改为u，并把路径中的\改为\\.
    不包含中文路径，直接拷路径，并在前面加r即可（linux路径不必加r）'''
    if not os.path.isdir(directory):
        print('Please check directory is right.')
        sys.exit(1)

    for folder in os.listdir(directory):
        path = os.path.join(directory,folder)
        if os.path.isfile(path): #遍历outmol和LOG文件
            if path.endswith('outmol'):
                print('%s\t%f' %(os.path.splitext(os.path.split(path)[1])[0], Dmol3).read_energy(path))
            elif path.endswith('LOG'):
                print('%s\t%f' %(os.path.splitext(os.path.split(path)[1])[0], Gaussian.read_energy(path)))
        else: #遍历子文件夹
            for file_name in os.listdir(path):
                if file_name.endswith('car'):
                    copy(os.path.join(path,file_name), directory)
                elif file_name.endswith('outmol'):
                    print('%s\t%f' %(file_name[:-7], Dmol3.read_energy(os.path.join(path,file_name))))
                elif file_name == 'CONTCAR':
                    c = read(os.path.join(path,file_name))
                    write_car(c, path+'.car')
                    print('%s\t%f' %(folder, Vasp.read_energy(path)))
                elif file_name.endswith('out') or file_name.endswith('LOG'):
                    c = read_gaussian_stru(os.path.join(path,file_name))
                    write_xyz(c, path+'.xyz')
                    print(os.path.basename(path) + '\t',end='')
                    e = Gaussian.read_energy(os.path.join(path,file_name))
                    result = Gaussian.zpe(os.path.join(path,file_name))
                    if result != None:
                        print('%f\t%f\t%f' %(e,result[0],result[1]))
                    else:
                        print(e)
                        
                        
if __name__ == '__main__':
    extract_result(r'E:\360disk\Cluster\B\B51++_Files\Documents\58')