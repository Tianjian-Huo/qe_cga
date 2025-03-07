'''io.py
该文件实现了常用文件的读写
作者: 赛琳伟
创建时间: 2020-11-14
'''
from .point import Point
from .atom import Atom,get_element,get_element_id
from . import utils
from .cluster import Cluster
from .crystal import Crystal,param6_cell,cell_param6
import os
import numpy as np
import time
from xml.sax.handler import ContentHandler
from xml.sax import parse

__all__ = ['read_cga_obj', 'write_xyz', 'read_xyz', 
           'write_car', 'read_car', 'append_arc', 'read_arc', 
           'write_xsd', 'read_xsd', 'read_outmol', 'outmol_structures',
           'write_poscar', 'write_contcar', 'read_vasp_stru', 'read_contcar',
           'write_gaussian_stru', 'read_gaussian_stru',  'read_gaussian_out_stru',
           'write_lammps_stru', 'read_lammps_stru',
           'write_cell', 'read_castep_stru',
           'write_abacus_stru', 'read_abacus_stru', 
           'read_cif', 
           'read']
                    
    
###############################################################################
#                 内部格式
###############################################################################
def read_cga_obj(f):
    '''从文件对象中读取CGA格式结构
    返回：团簇，下一行'''
    try:
        #读晶格常数
        line = f.readline()
        if line.startswith('pbc:'):
            line = line.split()
            c = Crystal()
            c.cell = param6_cell(*(float(x) for x in line[1:-1]))
            c.group = line[-1]
            line = f.readline()
        else:
            c = Cluster()
        #读能量,指纹
        content = line.split()
        if len(content) != 9:
            raise
        c.set_energy(float(content[0]))
        c.fingerprint = np.array([float(x) for x in content[1:]])
        #读坐标
        while True:
            line = f.readline()
            content = line.split()
            if len(content) != 4 or len(content[0]) > 2:
                break
            c.add_atom(Atom(float(content[1]), float(content[2]), float(content[3]), content[0]))
    except:
        print('file format error! error at line: ', line)
        raise utils.FileFormatErr(line)
    return c, line
    
###############################################################################
#                 Dmol
###############################################################################
def write_xyz(c, file_name):
    '''写xyz文件'''
    f = open(file_name, 'w')
    f.write(str(c.get_size()))
    f.write('\n\n')
    for a in c.atoms:
        f.write('%s\t%f\t%f\t%f\n' %(a.elem, a.x, a.y, a.z))
        
def read_xyz_obj(f):
    '''从xyz格式的文件对象中读取结构'''
    atom_num = int(f.readline())
    f.readline() #读空行
    atoms = []

    for _ in range(atom_num):
        line = f.readline().split()
        if len(line) != 4:
            raise
        atoms.append(Atom(float(line[1]), float(line[2]), float(line[3]), line[0]))
        
    return Cluster(atoms)       

def read_xyz(file_name):
    '''读取xyz文件
    返回：团簇'''
    if not os.path.isfile(file_name):
        print('no such file',file_name)
        return None

    try:
        return read_xyz_obj(open(file_name))
    except:
        print(file_name+' format error.')
        return None

def _write_car_obj(c, f):
    '''写car文件前3行之后的内容'''
    f.write('!DATE     ' + time.asctime()[4:] + '\n')
    if hasattr(c, 'cell'):
        va,vb,vc,alpha,beta,gamma = cell_param6(c.cell)
        f.write('PBC %9.4f%10.4f%10.4f%10.4f%10.4f%10.4f (%s)\n'
                %(va,vb,vc,alpha,beta,gamma,c.group))
    atom_format = '%s%-widthd%15.9f%15.9f%15.9f XXXX 1      xx      %-3s 0.000\n'
    count = {}
    for a in c.atoms:
        if count.get(a.elem) != None:
            count[a.elem] += 1
        else:
            count[a.elem] = 1
        f.write(atom_format.replace('width', str(5-len(a.elem)))
                %(a.elem, count[a.elem], a.x, a.y, a.z, a.elem))
    f.write('end\n')
    f.write('end\n')
    
def write_car(c, file_name):
    '''写car文件'''
    f = open(file_name, 'w')
    f.write('!BIOSYM archive 3\n')
    if hasattr(c, 'cell'):
        f.write('PBC=ON\n')
    else:
        f.write('PBC=OFF\n')
    f.write('CGA Generated CAR File\n')
    _write_car_obj(c, f)
    f.close()

def read_car_obj(f, pbc):
    '''读取car文件
    pbc：True or False
    返回：团簇或晶体'''
    line = f.readline()
    if line == '':
        return None
    line = f.readline()
    if line == '':
        return None
    if pbc: #读取晶格常数
        line = f.readline().split()
        cell = param6_cell(*(float(x) for x in line[1:7]))
        group = line[7][1:-1]
    else:
        cell = group = None
        
    atoms = []
    line = f.readline()
    while line:
        if line.startswith('end'):
            line = f.readline()
            if line.startswith('end'):
                break
        line = str(line).split()
        try:
            atoms.append(Atom(float(line[1]), float(line[2]), float(line[3]), line[-2]))
        except:
            print('file format error.')
            raise utils.FileFormatErr(line)
        line = f.readline()
        
    if pbc:
        return Crystal(atoms, cell, group)
    else:
        return Cluster(atoms)
    
def read_car(file_name):
    '''读取car文件
    返回：团簇或晶体'''
    if not os.path.isfile(file_name):
        print('no such file',file_name)
        return None

    f = open(file_name)
    f.readline()
    line = f.readline().split('=')
    pbc = True if line[1] == 'ON' else False
    return read_car_obj(f, pbc) 


def read_qe_obj(f, pbc):
    '''读取qe文件，返回团簇对象，其中包含原子坐标和名称'''
    
    atoms = []  # 用于存储原子信息
    line = f.readline()
    
    # 如果文件为空，返回 None
    if line == '':
        print("[警告] 文件为空，未读取到任何内容。")
        return None
    
    # 查找 "Atomic positions" 部分
    while line:
        line = f.readline()  # 跳过注释行
        while line.strip():  # 读取每个原子坐标
            data = line.split()
            if len(data) >= 4:
                element = data[0]  # 元素符号
                try:
                    # 检查 x, y, z 是否是有效的浮动数值
                    x = float(data[1])
                    y = float(data[2])
                    z = float(data[3])
                    
                    # 创建 Atom 对象
                    atom = Atom(x, y, z, element)
                    atoms.append(atom)  # 添加到原子列表
                except ValueError as e:
                    print(f"Error parsing line {line}: {e}")
                    raise ValueError("Error in parsing atomic positions")
            line = f.readline()
        break
    line = f.readline()

    return Cluster(atoms)


def read_qe(file_name):
    '''读取qe文件，返回包含原子坐标和名称的Cluster对象'''
    
    # 判断文件是否是 .out 文件，若是，改成 .qe 文件
    if file_name.endswith('.out'):
        file_name = file_name[:-4] + '.qe'  # 把 .out 改为 .qe
    
    if not os.path.isfile(file_name):
        print(f"no such file {file_name}")
        return None

    # 打开文件并传递给 `read_qe_obj`，不使用 `with open` 来避免自动关闭文件
    f = open(file_name, 'r')
    
    # 检查是否启用PBC
    line = f.readline().split('=')
    pbc = False  # 默认值为 False
    
    # 读取文件内容并传递给 `read_qe_obj`
    result = read_qe_obj(f, pbc)
    
    # 关闭文件
    f.close()
    
    return result





def append_arc(c, file_name):
    '''向arc文件追加一个结构'''
    if os.path.isfile(file_name):
        f = open(file_name, 'a')
    else: #文件不存在,先写文件头
        f = open(file_name, 'w')
        f.write('!BIOSYM archive 3\n')
        if hasattr(c, 'cell'):
            f.write('PBC=ON\n')
        else:
            f.write('PBC=OFF\n')
            
    f.write('%80f\n' %c.get_energy()) #写能量
    _write_car_obj(c, f)
    f.close()


def read_arc(file_name):
    '''返回：团簇或晶体列表'''
    if not os.path.isfile(file_name):
        print('no such file',file_name)
        return None

    f = open(file_name)
    f.readline()
    line = f.readline().split('=')
    pbc = True if line[1] == 'ON' else False
        
    stru_list = []    
    while True:
        c = read_car_obj(f, pbc)
        if c is None:
            break
        else:
            stru_list.append(c)
    
    return stru_list

class _XsdParser(ContentHandler):
    '''解析xsd文件的类'''
    def __init__(self, atom, bond):
        self.atom = atom
        self.bond = bond
    def startElement(self, name, attrs):
        if name == 'Atom3d':
            self.atom.append(attrs)
        elif name == 'Bond':
            self.bond.append(attrs)

def write_xsd(c, file_name, bonds=None):
    '''写xsd文件
    bonds: 成键矩阵。参见get_bonds函数'''
    properties = [('AngleAxisType','AngleBetweenPlanesBender','Enumerated'),
                ('AngleEnergy','ClassicalEnergyHolder','Double'),
                ('BeadDocumentID','MesoMoleculeSet','String'),
                ('BendBendEnergy','ClassicalEnergyHolder','Double'),
                ('BendTorsionBendEnergy','ClassicalEnergyHolder','Double'),
                ('BondEnergy','ClassicalEnergyHolder','Double'),
                ('EFGAsymmetry','Atom','Double'),
                ('EFGQuadrupolarCoupling','Atom','Double'),
                ('ElectrostaticEnergy','ClassicalEnergyHolder','Double'),
                ('FaceMillerIndex','GrowthFace','MillerIndex'),
                ('FacetTransparency','GrowthFace','Float'),
                ('FermiLevel','ScalarFieldBase','Double'),
                ('Force','Matter','CoDirection'),
                ('FrameFilter','Trajectory','String'),
                ('HarmonicForceConstant','HarmonicRestraint','Double'),
                ('HarmonicMinimum','HarmonicRestraint','Double'),
                ('HydrogenBondEnergy','ClassicalEnergyHolder','Double'),
                ('ImportOrder','Bondable','UnsignedInteger'),
                ('InversionEnergy','ClassicalEnergyHolder','Double'),
                ('IsBackboneAtom','Atom','Boolean'),
                ('IsChiralCenter','Atom','Boolean'),
                ('IsOutOfPlane','Atom','Boolean'),
                ('IsRepeatArrowVisible','ElectrodeWire','Boolean'),
                ('KineticEnergy','ClassicalEnergyHolder','Double'),
                ('LineExtentPadding','BestFitLineMonitor','Double'),
                ('LinkageGroupName','Linkage','String'),
                ('ListIdentifier','PropertyList','String'),
                ('NMRShielding','Atom','Double'),
                ('NonBondEnergy','ClassicalEnergyHolder','Double'),
                ('NormalMode','Bondable','Direction'),
                ('NormalModeFrequency','Bondable','Double'),
                ('NumScanSteps','LinearScan','UnsignedInteger'),
                ('OrbitalCutoffRadius','Bondable','Double'),
                ('PlaneExtentPadding','BestFitPlaneMonitor','Double'),
                ('PotentialEnergy','ClassicalEnergyHolder','Double'),
                ('QuantizationValue','ScalarFieldBase','Double'),
                ('RelativeVelocity','Matter','Direction'),
                ('RepeatArrowScale','ElectrodeWire','Float'),
                ('RestraintEnergy','ClassicalEnergyHolder','Double'),
                ('ScanEnd','LinearScan','Double'),
                ('ScanStart','LinearScan','Double'),
                ('SeparatedStretchStretchEnergy','ClassicalEnergyHolder','Double'),
                ('SimulationStep','Trajectory','Integer'),
                ('StretchBendStretchEnergy','ClassicalEnergyHolder','Double'),
                ('StretchStretchEnergy','ClassicalEnergyHolder','Double'),
                ('StretchTorsionStretchEnergy','ClassicalEnergyHolder','Double'),
                ('Temperature','ClassicalEnergyHolder','Double'),
                ('ThreeBodyNonBondEnergy','ClassicalEnergyHolder','Double'),
                ('TorsionBendBendEnergy','ClassicalEnergyHolder','Double'),
                ('TorsionEnergy','ClassicalEnergyHolder','Double'),
                ('TorsionStretchEnergy','ClassicalEnergyHolder','Double'),
                ('TotalEnergy','ClassicalEnergyHolder','Double'),
                ('Units','ScalarFieldBase','String'),
                ('ValenceCrossTermEnergy','ClassicalEnergyHolder','Double'),
                ('ValenceDiagonalEnergy','ClassicalEnergyHolder','Double'),
                ('VanDerWaalsEnergy','ClassicalEnergyHolder','Double'),
                ('_Stress','MatterSymmetrySystem','Matrix'),
                ('_TrajectoryStress','MatterSymmetrySystem','Matrix')]
    property_template = '\t\t<Property Name="%s" DefinedOn="%s" Type="%s"/>\n'
    atom_template = '\t\t\t<Atom3d ID="%d" Name="%s%d" UserID="%d"' +\
        ' DisplayStyle="Ball and Stick" XYZ="%f,%f,%f"' +\
        ' Connections="%s" BallSize="30" StickRadius="10" Components="%s"/>\n'
    bond_template = '\t\t\t<Bond ID="%d" Connects="%d,%d"/>\n'
    if bonds is None:
        bonds = c.get_bonds()

    #写文件头
    f = open(file_name, 'w')
    f.write('''<?xml version="1.0" encoding="latin1"?>
<!DOCTYPE XSD []>\n''')
    f.write('<XSD Version="7.0" WrittenBy="CGA">\n')
    f.write('\t<AtomisticTreeRoot ID="1" NumProperties="%d" NumChildren="1">\n' %len(properties))
    for p in properties: #写property属性
        f.write(property_template %(p[0],p[1],p[2]))
    f.write('\t\t<Molecule ID="2" NumChildren="%d">\n' %(len(bonds) + c.get_size()))

    #写原子
    for i in range(c.get_size()):
        connection = []
        for idx,bond in enumerate(bonds):
            if bond[0] == i or bond[1] == i:
                connection.append(str(idx+c.get_size()+3))
        connection = ','.join(connection)
        f.write(atom_template %(i+3, c.atoms[i].elem, i+1, i+1,
                                c.atoms[i].x, c.atoms[i].y, c.atoms[i].z,
                                connection, c.atoms[i].elem))

    #写连接
    for i in range(len(bonds)):
        f.write(bond_template %(i+c.get_size()+3, bonds[i][0]+3, bonds[i][1]+3))
    #写结尾
    f.write('''\t\t</Molecule>
\t</AtomisticTreeRoot>
</XSD>''')
    f.close()

def read_xsd(file_name):
    '''读xsd文件
    参数file_name: 字符串类型,要读取的文件名
    返回: 团簇，邻接矩阵'''
    if not os.path.isfile(file_name):
        print('no such file',file_name)
        return None, None

    #解析xsd
    atom = []
    bonds = []
    parse(file_name, _XsdParser(atom, bonds))

    #读取元素,坐标
    id_ = []
    atoms = []
    for a in atom:
        if 'XYZ' not in a: #删除原子后，连接还在
            continue
        coo = a['XYZ'].split(',')
        try:
            atoms.append(Atom(float(coo[0]), float(coo[1]), float(coo[2]), a['Components']))
        except:
            raise utils.FileFormatErr
        id_.append(int(a['ID']))

    #读取连接
    adj = np.zeros((len(atoms), len(atoms)))
    for b in bonds:
        endpoint = b['Connects'].split(',')
        try: #删除原子后，连接还在
            start = id_.index(int(endpoint[0]))
            end = id_.index(int(endpoint[1]))
        except:
            continue
        adj[start,end] = 1
        adj[end,start] = 1

    return Cluster(atoms), adj

def read_outmol(file_name):
    '''从outmol中读取结构.仅支持团簇
    未修改
    返回：团簇'''
    if not os.path.isfile(file_name):
        print('no such file', file_name)
        return None
    step = 4096
    f = open(file_name)
    f.seek(0, 2)
    bottom_pos = f.tell()

    flag = 0
    while True: #分段往前读
        f.seek(-min(step+80, bottom_pos), 1)
        while f.tell() < bottom_pos:
            line = f.readline()
            if line.find('Final Coordinates') != -1 or \
                line.find('Input Coordinates') != -1: #后者表示未算完
                flag = 1
                f.readline(),f.readline()
            elif line.startswith('$coordinates'): #单点能
                flag = 2
            if flag != 0:
                atoms = []
                while True:
                    data = f.readline().split()
                    if len(data) < 4:
                        break
                    try:
                        atoms.append(Atom(float(data[-3]),float(data[-2]),float(data[-1]), data[-4]))
                    except:
                        raise utils.FileFormatErr
                if flag == 2:
                    Cluster(atoms).zoom(0.529177249)
                return
        if f.tell() < step:
            break
        else:
            f.seek(-step, 1)
            bottom_pos = f.tell()
    return Cluster(atoms)

def outmol_structures(file_name):
    '''未修改
    返回：团簇列表'''
    result = []
    f = open(file_name)
    while True:
        line = f.readline()
        if not line:
            break
        elif line.find('Input Coordinates') != -1:
            f.readline(), f.readline()
            c = Cluster()
            while True:
                line = f.readline().split()
                if len(line) != 5:
                    break
                c.add_atom(Point(float(line[-3]),float(line[-2]),float(line[-1])),line[1])
        elif line.startswith('opt=='):
            line = line.split()
            try:
                int(line[1])
            except:
                pass
            else:
                c.set_energy(float(line[2]))
                result.append(c)
    return result
    
###############################################################################
#                 VASP
###############################################################################
def write_vasp(c, f, direct=False, fix=[]):
    '''c：结构
    f：文件对象
    direct：分数坐标还是笛卡尔坐标
    fix：那些原子需要固定'''
    elements = c.get_elements_count()
    elements_sort = list(elements.keys())
    elements_sort.sort(key=lambda e: get_element_id(e))

    #写标题
    f.write('_'.join([e+str(elements[e]) for e in elements_sort]) +  '\n')

    #写晶格常数
    f.write('1.0\n')
    if hasattr(c, 'cell'):
        cell = c.cell
    else:
        cell = np.diag(np.sum(c.coord_matrix(axis=0)) + 15)
        c2 = c.deepcopy()
        c2.move(c2.get_center() +
               Point(np.linalg.norm(c2.cell[0]),
                     np.linalg.norm(c2.cell[1]),
                     np.linalg.norm(c2.cell[2])))
    for v in cell:
        for x in v:
            f.write(' %15.6f'%x)
        f.write('\n')
        
    #写元素
    f.write(' '.join(elements_sort) + '\n')
    #写原子数
    f.write(' '.join([str(elements[e]) for e in elements_sort]) +  '\n')

    if fix:
        f.write('Selective dynamics\n')
    #写坐标
    coord = c2.coord_matrix()
    if direct:
        f.write('Direct\n')
        #坐标需要转换
        coord = np.linalg.solve(cell.T, coord.T).T #笛卡尔坐标转分数坐标
    else:
        f.write('Cartesian\n')
    for e in elements_sort:
        for i in range(c2.get_size()):
            if c2.atoms[i].elem == e:
                f.write('%f\t%f\t%f' %tuple(coord[i]))
            flag = ''
            if fix: #固定坐标
                flag += 'T' if fix[i]%2 else 'F'
                flag += 'T' if fix[i]//2%2 else 'F'
                flag += 'T' if fix[i]//4 else 'F'
                f.write(' '+flag)                        
            f.write('\n')

def write_poscar(c, folder_name='', direct=False, fix=[]):
    '''写vasp的坐标文件
    参数folder_name: 写到哪个文件夹里,默认当前文件夹'''
    f = open(os.path.join(folder_name,'POSCAR'), 'w')
    write_vasp(c, f, direct, fix)

def write_contcar(c, folder_name='', direct=False, fix=[]):
    '''写vasp的坐标文件
    参数folder_name: 写到哪个文件夹里,默认当前文件夹'''
    f = open(os.path.join(folder_name,'CONTCAR'), 'w')
    write_vasp(c, f, direct, fix)

def read_vasp_stru(file_name):
    '''从POSCAR或CONTCAR文件中提取结构
    返回：晶体或团簇。自动判断是晶体还是团簇'''
    if not os.path.isfile(file_name):
        print('no such file', file_name)
        return None, None
    f = open(file_name)
    
    f.readline() #title
    try:
        scale = float(f.readline())
        cell = np.array([[float(x) for x in f.readline().split()],
                         [float(x) for x in f.readline().split()],
                         [float(x) for x in f.readline().split()]])
        cell *= scale
        elem_types = f.readline().split()
        atom_nums = [int(num) for num in f.readline().split()]
    except:
        raise utils.FileFormatErr()
        
    fraction = False #是否分数坐标
    while True:
        line = f.readline()
        if line.startswith('D'): #Direct
            fraction = True
            break
        elif line.startswith('C'): #Cartesion
            break
            
    coord = []
    for n in atom_nums:
        for _ in range(n):
            line = f.readline().split()
            try:
                coord.append([float(line[0]), float(line[1]),float(line[2])])
            except:
                raise utils.FileFormatErr()
    coord = np.array(coord)
    if fraction: #分数坐标转笛卡尔坐标
        coord = np.dot(coord, cell)
        
    atoms = []
    i = 0
    for n,e in zip(atoms, elem_types):
        atoms.append(Atom(*tuple(coord[i]), e))
        
    c = Crystal(atoms, cell)
    if c.is_cluster():
        return Cluster(atoms)
    else:
        return c

def read_contcar(folder_name=''):
    '''从CONTCAR文件中提取坐标
    参数folder_name: 从哪个文件夹里读contcar,默认为当前文件夹
    返回：团簇或晶体'''
    return read_vasp_stru(os.path.join(folder_name,'CONTCAR'))

###############################################################################
#                 Gaussian
###############################################################################
def write_gaussian_stru(c, file_name, params={}):
    '''写gaussian输入文件，后缀为gif
    需要传入并行的核数
    仅支持团簇'''
    f = open(file_name, 'w')
    if 'core_num' in params:
        f.write('%%nprocshared=%d\n'%params['core_num'])
    f.write('%nosave\n%mem=48GB\n')
    f.write('%%chk=%s.chk\n'%os.path.split(file_name)[1][:-4])

    param_line = '#p'
    if 'task' in params: #opt, freq,
        param_line += params['task']
    if 'basis' in params: #基组
        param_line += params['basis']
    f.write(param_line+'\n')

    f.write('\ntemp\n\n')
    charge = params['charge'] if 'charge' in params else 0
    spin = 2 if c.get_spin() else 1
    f.write('%d %d \n' %(charge, spin)) #charge spin

    for a in c.atoms:
        f.write('%s %15.6f%15.6f%15.6f\n' %(a.elem, a.x, a.y, a.z))
    f.write('\n\n\n')
    f.close()
        
def read_gaussian_stru(file_name):
    '''读gaussian结构文件（输入文件），后缀为gif
    仅支持团簇
    返回：团簇'''
    f = open(file_name)
    while True:
        line = f.readline()
        if line == '':
            break
        if line.startswith('#p'):
            for _ in range(5):
                f.readline()
            break
    atoms = []
    while True:
        line = f.readline().split()
        if len(line) != 4:
            break
        atoms.append(Atom(float(line[1]), float(line[2]), float(line[3]), line[0]))
    return Cluster(atoms)

def read_gaussian_out_stru(file_name):
    '''读取gaussian输出文件中的结构
    仅支持团簇
    返回：团簇'''
    if not os.path.isfile(file_name):
        print('no such file', file_name)
        return None
    f = open(file_name)
    f.seek(0, 2)
    end_pos = f.tell()
    step = 4096

    while True: #分段往前读
        start_pos = max(end_pos-step, 0)
        f.seek(start_pos, 0)
        while f.tell() < end_pos:
            line = f.readline()
            if line.find('Standard orientation') != -1:
                f.readline(),f.readline(),f.readline(),f.readline()
                atoms = []
                while True:
                    line = f.readline().split()
                    if len(line) != 6:
                        return Cluster(atoms)
                    try:
                        atoms.append(Atom(float(line[-3]), float(line[-2]), float(line[-1]), get_element(int(line[1]))))
                    except:
                        raise utils.FileFormatErr(line)
        end_pos = start_pos
    return None

###############################################################################
#                 Lammps
###############################################################################
def write_lammps_stru(c, file_name=''):
    f = open(file_name, 'w')
    #写标题
    elements_count = c.get_elements_count()
    f.write(file_name +  ' cluster, written by CGA.\n')
    f.write('%d atoms\n'%c.get_size())
    f.write('%d atom types\n' %len(elements_count))

    x_max = max([a.x for a in c.atoms]) + 8
    x_min = min([a.x for a in c.atoms]) - 8
    y_max = max([a.y for a in c.atoms]) + 8
    y_min = min([a.y for a in c.atoms]) - 8
    z_max = max([a.z for a in c.atoms]) + 8
    z_min = min([a.z for a in c.atoms]) - 8
    f.write('%f %f xlo xhi\n %f %f ylo yhi\n%f %f zlo zhi\n' \
            %(x_min, x_max, y_min, y_max, z_min, z_max))
    f.write('\nAtoms\n\n')

    for i,e in enumerate(elements_count.keys()):
        for j,a in enumerate(c.atoms):
            if a.elem == e:
                f.write('%d %d '%(j+1,i+1) + str(a.coordinate()) + '\n')
        
def read_lammps_stru(file_name):
    '''文件名固定为dump.dat
    返回：团簇'''
    if not os.path.isfile(file_name):
        print('no such file',file_name)
        return None
    else:
        f = open(file_name)
    
    result = None
    while True:
        try:
            out = read_xyz_obj(f)
        except:
            return result
        else:
            result = out
    return result

###############################################################################
#                 其他格式
###############################################################################
def read_cif(file_name):
    '''读取Crystallographic Information File格式文件，该格式仅支持晶体
    返回：晶体'''
    f = open(file_name)
    c = Crystal()
    while True:
        line = f.readline()
        if line.startswith('_symmetry_space_group_name_H-M'):
            c.group = line.split()[1][1:-1]
            break
            
    while True:
        line = f.readline()
        if line.startswith('_cell_length_a'):
            va = float(line.split()[1])
            break
    vb = float(f.readline().split()[1])
    vc = float(f.readline().split()[1])
    alpha = float(f.readline().split()[1])
    beta = float(f.readline().split()[1])
    gamma = float(f.readline().split()[1])
    c.cell = param6_cell(va,vb,vc,alpha,beta,gamma)
    
    while True:
        line = f.readline()
        if line.startswith('_atom_site_occupancy'):
            while True:
                line = f.readline().split()
                if line[0] == 'loop_':
                    break
                c.add_atom(Atom(float(line[2])*va, float(line[3])*vb, float(line[4])*vc, line[1]))
            break
    
    return c

def write_cell(c, file_name):
    '''castep输入文件'''
    sample = open('castep.cell').read().split('%ENDBLOCK POSITIONS_FRAC')[1]
    f = open(file_name, 'w')
    f.write('%BLOCK LATTICE_CART\n')
    f.write('%f 0.0 0.0\n 0.0 %f 0.0\n0.0 0.0 %f\n' %(c.a,c.b,c.c))
    f.write('%ENDBLOCK LATTICE_CART\n')
    f.write('%BLOCK POSITIONS_FRAC\n')
    for p in zip(c.atoms):
        f.write('%3s%21f%21f%21f\n' %(p.elem,p.x,p.y,p.z))
    f.write('%ENDBLOCK POSITIONS_FRAC\n')
    f.write(sample)

def read_castep_stru(file_name):
    '''返回：晶体'''
    c = Crystal()
    _, f = utils.rfind_file(file_name, contain_word='Lattice parameters')
    line = f.readline()
    va = float(line.split()[2])
    alpha = float(line.split()[-1])
    line = f.readline()
    vb = float(line.split()[2])
    beta = float(line.split()[-1])
    line = f.readline()
    vc = float(line.split()[2])
    gamma = float(line.split()[-1])
    c.cell = param6_cell(va,vb,vc,alpha,beta,gamma)

    for _ in range(11):
        f.readline()
    while True:
        line = f.readline().split()
        if line[0][:3] == 'xxx':
            break
        c.add_atom(Atom(float(line[-3]), float(line[-3]), float(line[-3]), line[1]))
        
    return c

def write_abacus_stru(c, folder_name='', relax_list=None):
    '''写ABACUS结构文件
    参数folder_name: 写到哪个文件夹里,默认当前文件夹'''
    head = open('STRU').read().split('LATTICE_VECTORS')[0]
    f = open(os.path.join(folder_name,'STRU'), 'w')
    #写标题
    f.write(head)

    a = max([m.x for m in c.atoms]) - min([m.x for m in c.atoms])
    b = max([m.y for m in c.atoms]) - min([m.y for m in c.atoms])
    c = max([m.z for m in c.atoms]) - min([m.z for m in c.atoms])
    const = (max(a,b,c) + 10) / 1.88972687777

    #写晶格常数
    f.write('LATTICE_VECTORS\n')
    f.write('%f\t0\t0\n0\t%f\t0\n0\t0\t%f\n\n'
            %(const, const, const))

    #写坐标
    f.write('ATOMIC_POSITIONS\nDirect\n\n')
    elements = c.get_elements_count()
    for e in elements:
        f.write(e+'\n')
        f.write('0\n')
        f.write(str(elements[e])+'\n')
        for i in range(c.get_size()):
            if c.atoms[i].elem == e:
                f.write(str(c.atoms[i]/const + Point(0.5,0.5,0.5)))
                if relax_list is None or i in relax_list:
                    f.write(' 1 1 1\n')
                else:
                    f.write(' 0 0 0\n')

def read_abacus_stru(file_name):
    '''从ABACUS的STRU_ION_D文件中提取坐标
    返回：团簇'''
    if not os.path.isfile(file_name):
        print('no such file', file_name)
        return
    f = open(file_name)
    while True:
        line = f.readline()
        if line.find('LATTICE_CONSTANT') != -1:
            crystal_length = float(f.readline())
            break
    while True:
        line = f.readline()
        if line.find('LATTICE_VECTORS') != -1:
            break
    a = float(f.readline().split()[0]) * crystal_length
    b = float(f.readline().split()[1]) * crystal_length
    c = float(f.readline().split()[2]) * crystal_length

    atoms = []
    while True:
        line = f.readline()
        if line == '':
            break
        if len(line.split()) > 0 and get_element_id(line.split()[0]) > 0:
            e = line.split()[0]
            f.readline()
            n = int(f.readline().split()[0])
            for _ in range(n):
                line = f.readline().split()
                atoms.append(Atom(float(line[0])*a, float(line[1])*b,float(line[2])*c, e))
    return Cluster(atoms)


def read(file_name):
    '''读取任意格式文件'''
    name, ext = os.path.splitext(file_name)
    if ext == '':
        path,name = os.path.split(file_name)
        if name == 'CONTCAR' or name == 'POSCAR':
            return read_vasp_stru(file_name)
        else:
            print('do not support this file:',file_name)
        return None
    ext = ext[1:]
    if ext == 'xyz':
        return read_xyz(file_name)
    elif ext == 'car':
        return read_car(file_name)
    elif ext == 'gjf' or ext == 'LOG': #暂时没有gauss晶体的需求
        return read_gaussian_stru(file_name)
    elif ext == 'cif':
        return read_cif(file_name)
    elif ext == 'xsd': #xsd文件暂时只支持团簇
        return read_xsd(file_name)[0]
    elif ext == 'outmol': #outmol文件暂时只支持团簇
        return read_outmol(file_name)
    elif ext == 'out':  # 添加对 .qe 文件的支持
        return read_qe(file_name)
    elif name.endswith('dump') and ext=='dat':
        return read_lammps_stru(file_name)
    else:
        print('do not support this file:',file_name)
        return None
