import os
import sys
import tempfile
import shutil

from shutil import rmtree,copyfile
import configparser
from configparser import ConfigParser
import re
import time
import subprocess
import platform
import psutil
from signal import SIGTERM
from periodictable import elements
import random
import string
import json

import numpy as np
from . import utils
from .point import Point
from .io import *
from .cluster import Cluster
try:
    import scipy.optimize as opt
    import torch
    from ..gcn.gcn import GCN,GCNF
    import warnings
    warnings.simplefilter("ignore") #关闭经验势计算中的警告
except:
    #print('warning: 没有scipy或torch模块，不使用经验势优化可忽略')
    pass


_accelrys_path = ''
_ms_version = ''



def set_accelrys_path():
    global _accelrys_path
    if platform.system() == 'Windows':
        _accelrys_path = Abinitio.cfg.get('optimize', 'accelrys_path_windows')
    else:
        _accelrys_path = Abinitio.cfg.get('optimize', 'accelrys_path_linux')

    global _ms_version
    _ms_version = re.compile(r'\d+\.?\d*').findall(_accelrys_path.split('etc')[0])[-1]
    if '.' in _ms_version:
        _ms_version = int(float(_ms_version)*10)
    else:
        _ms_version = int(_ms_version)


# Utility Functions
def load_pseudo_files(pseudo_dir):
    pseudo_files = {}
    try:
        for file in os.listdir(pseudo_dir):
            if file.endswith((".UPF", ".upf")):
                element = file.split(".")[0]
                pseudo_files[element] = file
        print(f"[信息] 共加载 {len(pseudo_files)} 个赝势文件.")
    except Exception as e:
        print(f"[错误] 无法加载赝势文件: {e}")
    return pseudo_files

def get_mass(element):
    try:
        return getattr(elements, element).mass
    except AttributeError:
        return 0.0

def calculate_lattice_from_positions(atomic_positions, offset=0):
    positions = np.array([list(map(float, pos[1:])) for pos in atomic_positions])
    max_pos = np.max(positions, axis=0)
    min_pos = np.min(positions, axis=0)
    return max_pos - min_pos + offset
    

def parse_molecule_from_config(config_path):
    config = ConfigParser()
    config.read(config_path)
    elements = config["cluster"]["element"].split()
    atom_counts = list(map(int, config["cluster"]["atom_num"].split()))
    return [(elements[i], atom_counts[i]) for i in range(len(elements))]


def generate_atoms(parsed_elements):
    """
    根据解析的分子式生成原子列表，包含原子的 symbol 和初始坐标。
    :param parsed_elements: 元素和数量的列表，例如 [("Ca", 1), ("Ni", 2)]
    :return: 原子列表，每个原子包含 symbol 和初始坐标 (x, y, z)
    """
    atoms = []
    offset = 1.0  # 坐标偏移量
    for symbol, count in parsed_elements:
        for i in range(count):
            x = offset * i
            y = offset * i
            z = offset * i
            atoms.append({"symbol": symbol, "x": x, "y": y, "z": z})
    print(f"[调试] 生成的 atoms: {atoms}")
    return atoms


def write_qe_input(cluster, file_name, pseudo_dir, pseudo_mapping, config_file):
    # 读取配置文件中的参数
    if isinstance(config_file, str):
        # 如果 config_file 是字符串（路径），则读取文件
        with open(config_file, "r") as f:
            config = json.load(f)
    else:
        raise ValueError(f"Expected config_file to be a string path, but got {type(config_file)}")
    
    # 从配置文件中读取参数，如果没有则使用默认值
    forc_conv_thr = config.get('forc_conv_thr', 1e-3)
    nstep = config.get('nstep', 30)
    ecutwfc = config.get('ecutwfc', 40.0)
    ecutrho = config.get('ecutrho', 240.0)
    electron_maxstep = config.get('electron_maxstep', 30)
    conv_thr = config.get('conv_thr', '1e-7')
    degauss = config.get('degauss', 0.001)
    nspin = config.get('nspin', 2)
    smearing = config.get('smearing', 'gaussian')
    mixing_beta = config.get('mixing_beta', 0.3)

    element_counts = {}
    for atom in cluster.atoms:
        element_counts[atom.elem] = element_counts.get(atom.elem, 0) + 1
    prefix = "".join([f"{k}{v}" for k, v in element_counts.items()])

    atomic_positions = [(atom.elem, atom.x, atom.y, atom.z) for atom in cluster.atoms]
    lattice_constants = calculate_lattice_from_positions(atomic_positions)

    # 创建一个新的输入文件
    with open(file_name, "w") as out:
        # 写入 CONTROL Section
        out.write("&control\n")
        out.write(f"   prefix = '{prefix}'\n")
        out.write("   calculation = 'relax'\n")
        out.write("   restart_mode = 'from_scratch'\n")
        out.write("   wf_collect = .false.\n")
        out.write("   tprnfor = .true.\n")
        out.write("   outdir = './scratch'\n")
        out.write(f"   pseudo_dir = '{pseudo_dir}'\n")
        out.write(f"   forc_conv_thr = {forc_conv_thr}\n")
        out.write(f"   nstep = {nstep}\n")
        out.write("/\n\n")

        # 写入 SYSTEM Section
        out.write("&system\n")
        out.write("   ibrav = 0\n")
        out.write(f"   nat = {len(atomic_positions)}\n")
        out.write(f"   ntyp = {len(set([p[0] for p in atomic_positions]))}\n")
        out.write(f"   ecutwfc = {ecutwfc}\n")
        out.write(f"   ecutrho = {ecutrho}\n")
        out.write("   occupations = 'smearing'\n")
        out.write(f"   smearing = '{smearing}'\n")
        out.write(f"   degauss = {degauss}\n")
        out.write(f"   nspin = {nspin}\n")
        out.write(f"   starting_magnetization(1) = 0.001\n")
        out.write("/\n\n")

        # 写入 ELECTRONS Section
        out.write("&electrons\n")
        out.write("   scf_must_converge = .false.\n")
        out.write(f"   electron_maxstep = {electron_maxstep}\n")
        out.write(f"   conv_thr = {conv_thr}\n")
        out.write(f"   mixing_mode = 'plain'\n")
        out.write(f"   mixing_beta = {mixing_beta}\n")
        out.write("   diagonalization = 'david'\n")
        out.write("/\n\n")

        # 写入 IONS Section
        out.write("&ions\n")
        out.write("   ion_dynamics = 'bfgs'\n")
        out.write("/\n\n")

        # 写入 K_POINTS Section
        out.write("K_POINTS Gamma\n\n")

        # ATOMIC_SPECIES Section
        out.write("ATOMIC_SPECIES\n")
        for elem in set([p[0] for p in atomic_positions]):
            pseudo_file = pseudo_mapping.get(elem)
            if not pseudo_file:
                raise ValueError(f"[ERROR] Missing pseudo file for {elem}.")
            out.write(f"  {elem} {get_mass(elem):.2f} {pseudo_file}\n")
        out.write("\n")

        # CELL_PARAMETERS Section
        out.write("CELL_PARAMETERS (angstrom)\n")
        if config.get("CELL_PARAMETERS") == "offset":
            offset = config.get("offset", 10)
            out.write(f"  {lattice_constants[0] + offset:.6f} 0.0 0.0\n")
            out.write(f"  0.0 {lattice_constants[1] + offset:.6f} 0.0\n")
            out.write(f"  0.0 0.0 {lattice_constants[2] + offset:.6f}\n")
        else:
            lattice_value = config.get("CELL_PARAMETERS", 10)
            out.write(f"  {lattice_value:.6f} 0.0 0.0\n")
            out.write(f"  0.0 {lattice_value:.6f} 0.0\n")
            out.write(f"  0.0 0.0 {lattice_value:.6f}\n")
        out.write("\n")

        # ATOMIC_POSITIONS Section
        out.write("ATOMIC_POSITIONS (angstrom)\n")
        for atom in atomic_positions:
            out.write(f"  {atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n")

        print(f"[信息] 输入文件: {file_name}")




class Abinitio(object):
    optimize_name = '' #第一性原理软件名称
    abinitio_cmd = ''
    write_stru = 'write_car' #写输入结构函数名
    wait_run_time = 600
    wait_start_time = 100
    wait_stru_time = 30
    wait_energy_time = 10
    core_num = 1
    keep_temp_folder = False
    temp_path = '.'

    @classmethod
    def set_path(cls, folder_name):
        '''设置临时文件夹的位置'''
        cls.temp_path = folder_name

    @classmethod
    def read_cmd(cls):
        pass

    @classmethod
    def config(cls):
        cfg = ConfigParser()
        try:
            if os.path.isfile('config.ini'):
                cfg.read('config.ini')
            else:
                #cfg.read(os.path.dirname(inspect.stack()[0][1])+'/config.ini')
                cfg.read(os.path.join(os.path.split(__file__)[0],'config.ini'))
        except:
            print('no config.ini')
            os._exit(0)
        #读取并行计算的核数
        cls.core_num = cfg.getint('optimize', 'core_num')
        cls.cfg = cfg
        #第一性原理计算方法
        method = cfg.get('optimize', 'method')
        set_optimize_method(method)
        #是否保留临时文件夹
        cls.keep_temp_folder = cfg.getboolean('optimize', 'keep_temp_folder')

    def __init__(self, cluster, folder_name=None, **kw):
        self.cluster = cluster
        self.folder_name = folder_name #在哪个文件夹里计算，不指定则用临时文件夹
        self.use_temp_folder = False
        self.write_stru_ext_param = None
        self.process = None

    def __call__(self):
        if self.abinitio_cmd == '':
            self.read_cmd()
        self.cluster.set_energy(0.)
        self.mkdir()
        if self.write_stru_ext_param is None:
            eval(self.write_stru)(self.cluster, self.get_in_stru_name())
        else:
            eval(self.write_stru)(self.cluster, self.get_in_stru_name(), self.write_stru_ext_param)
        self.prepare_input()
        if hasattr(self, 'fix_coord') and self.fix_coord is not None:
            self.write_fix_coord()
        self.call_abinitio()
        utils.write_log('job start at: ' + time.ctime(time.time()))
        utils.write_log('folder name: ' + self.folder_name)
        self._wait()
        utils.write_log('job end at: ' + time.ctime(time.time()))
        self.get_stru()
        self.get_energy()
        if self.use_temp_folder and not self.keep_temp_folder:
            self.rmdir()

    def mkdir(self):
        '''创建文件夹，如果不指定文件夹名，则使用临时文件夹'''
        if self.folder_name is None:
            # 如果未指定文件夹名，使用临时文件夹
            self.folder_name = tempfile.mktemp(prefix=self.optimize_name+'_', dir=self.temp_path)
            self.use_temp_folder = True

        # 确保文件夹存在
        os.makedirs(self.folder_name, exist_ok=True)


    def rmdir(self):
        '''删除临时文件夹'''
        if self.use_temp_folder:
            try:
                rmtree(self.folder_name)
            except:
                utils.write_log('del temp folder "%s" failure.' %self.folder_name)

    def get_in_stru_name(self):
        '''获取写入结构的文件名'''
        pass

    def get_out_stru_name(self):
        '''获取计算后的结构文件名'''
        pass

    def get_out_file_name(self):
        '''获取能量文件名'''
        pass

    def prepare_input(self):
        '''准备input文件'''
        pass

    def write_fix_coord(self):
        '''优化时固定哪些坐标'''
        pass

    def call_abinitio(self):
        '''调用第一性原理软件'''
        pass

    def _wait(self):
        '''等待第一性原理计算结束'''
        file_name = self.get_out_file_name()
        #先监视一段时间，若还没有产生输出文件，则有错。
        for _ in range(int(self.wait_start_time)*100):
            time.sleep(0.01)
            if os.path.isfile(file_name):
                break
        else:
            utils.write_log('can not find %s, maybe task start failure.' %file_name)
            return
        #继续监视，wait_time时间后输出文件无更新则认为失败
        lastSize = os.path.getsize(file_name)
        count = 0
        while self.process.poll() == None:
            time.sleep(0.01)
            count += 1
            if count > self.wait_run_time / 0.01:
                utils.write_log('out file not update in %d seconds.'%self.wait_run_time)
                self.kill()
                return
            size = os.path.getsize(file_name)
            if size > lastSize:
                lastSize = size
                count = 0
        if self.process.poll() != 0:
            utils.write_log('task terminate abnormally, return value: %d' %self.process.poll())

    def get_stru(self):
        '''获取计算后的结构'''
        #可能传输结果需要时间，这里等待一会
        for _ in range(int(self.wait_stru_time)*100):#获取结构
            if os.path.isfile(self.get_out_stru_name()):
                break
            time.sleep(0.01)
        else:
            utils.write_log('can not find structure file.')
        self.cluster.deepcopy(read(self.get_out_stru_name()))
        if hasattr(self.cluster, 'after_read'):
            self.cluster.after_read()

    def get_energy(self):
        '''获取计算后的能量'''
        #可能传输结果需要时间，这里等待一会
        for _ in range(int(self.wait_energy_time)*100):
            if os.path.isfile(self.get_out_file_name()):
                break
            time.sleep(0.01)
        else:
            utils.write_log('can not find energy file.')
            return
        #能量所在文件出现，但结果未传完，这里等待一会
        for _ in range(int(self.wait_energy_time)*100):
            self.cluster.set_energy(self.read_energy(self.get_out_file_name()))
            if self.cluster.get_energy() != 0.:
                break
            time.sleep(0.01)
        else:
            utils.write_log('can not get energy value.')
        print(f"[DEBUG] 准备在abinitio中返回能量值")

    @staticmethod
    def read_energy(file_name):
        pass

    @staticmethod
    def check_symmetry_on(file_name=''):
        return True

    def kill(self):
        if platform.system() == 'Windows':
            parent = psutil.Process(self.process.pid)
            children = parent.children(recursive=True)
            children.append(parent)
            for p in children:
                p.send_signal(SIGTERM)
        elif platform.system() == 'Linux':
            os.killpg(self.process.pid, SIGTERM)

class QE(Abinitio):
    optimize_name = 'qe'
    
    # 读取 config.ini 配置文件
    config_file = os.path.join(os.getcwd(), "search", "config.ini")
    config = configparser.ConfigParser()
    if os.path.exists(config_file):
        config.read(config_file)
    else:
        raise FileNotFoundError(f"[错误] 未找到配置文件: {config_file}")

    # 从配置文件中获取参数
    abinitio_cmd = config.get("optimize", "qe_dir", fallback="/home/cast/software/q-e-qe-7.2/bin/pw.x")

    def __init__(self, cluster, folder_name=None, **kw):
        super(QE, self).__init__(cluster, folder_name)
        if not folder_name:
            raise ValueError("必须提供文件夹路径用于QE计算。")

        self.folder_name = os.path.abspath(folder_name)
        parent_dir = os.path.dirname(self.folder_name)
        config_file = os.path.join(parent_dir, "config.ini")

        # 读取 config.ini 配置文件
        config = configparser.ConfigParser()
        if os.path.exists(config_file):
            config.read(config_file)
        else:
            raise FileNotFoundError(f"[错误] 未找到配置文件: {config_file}")

        # 从配置文件中获取参数
        self.pseudo_dir = config.get("optimize", "qe_pseudo_dir", fallback="/home/cast/users/htj/pseudo")
        self.core_num = config.getint("optimize", "core_num", fallback=32)

    def prepare_input(self):
        """准备 QE 输入文件"""
        # 加载伪势文件
        pseudo_mapping = load_pseudo_files(self.pseudo_dir)
        # 获取输入结构文件名
        in_stru_name = self.get_in_stru_name()
                # 获取上一级目录的路径
        parent_dir = os.path.dirname(self.folder_name)
        
        # 获取 qein.json 配置文件路径（仅用于计算参数，不包含路径信息）
        config_file = os.path.join(parent_dir, "qein.json")

        # 调用 write_qe_input，并传递 config_file 路径
        write_qe_input(self.cluster, in_stru_name, self.pseudo_dir, pseudo_mapping, config_file)


    def get_in_stru_name(self):
        """返回输入文件的路径"""
        return os.path.join(self.folder_name, "qe.in")
    
    def call_abinitio(self):
        """运行 QE 并生成 .qe 文件"""
        input_file = self.get_in_stru_name()
        output_file = self.get_out_stru_name()
        print(f"[信息] 输出文件: {output_file}")
        command = f"mpirun -np {self.core_num} {self.abinitio_cmd} < {input_file} > {output_file}"
        self.process = subprocess.Popen(
            command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid,
            stderr=subprocess.PIPE
        )

        try:
            sys.stdout.write("[信息] 正在等待 QE 进程完成...")  
            sys.stdout.flush()  # 立刻输出到终端
            self._wait()
            sys.stdout.write("\r" + " " * 50 + "\r")  # 用空格清空行，并回到行首
            sys.stdout.flush()

        except KeyboardInterrupt:
            self.terminate_process(self.process.pid)
            sys.exit(1)

        # 获取能量值，如果电子步未收敛，直接返回
        energy = self.get_energy()
        if energy == 0:
            return 0

        result_file = self.parse_qe_results(output_file)
        if result_file:
            return result_file
        else:
            return 0



    def _wait(self):
        """等待 QE 计算完成"""
        output_file = self.get_out_stru_name()
        if self.process is not None:
            while self.process.poll() is None:
                time.sleep(1)
        else:
            raise RuntimeError("未启动 QE 进程。请确保 call_abinitio 正确调用。")

        if os.path.getsize(output_file) == 0:
            raise ValueError(f"输出文件为空: {output_file}")

    def terminate_process(self, pid):
        """终止进程及其所有子进程"""
        parent = psutil.Process(pid)
        for child in parent.children(recursive=True):
            child.terminate()
        parent.terminate()
        for child in parent.children(recursive=True):
            child.wait()
        parent.wait()

    def get_energy(self):
        """从 QE 输出中提取总能量"""
        output_file = self.get_out_stru_name()
        energy = self.read_energy(output_file)

        if energy is None:
            raise ValueError(f"未能在输出文件中找到能量值: {output_file}")

        # 只在这里警告一次
        if self.process.returncode == 3:
            print("[警告] 离子步未收敛，读取最后的能量值")  

        self.cluster.set_energy(energy)
        return energy

    def read_energy(self, output_file):
        """从 QE 输出文件中读取能量值"""
        energy = None

        # 检查进程是否正常执行
        if self.process.returncode == 2:
            print("[错误] 电子步未收敛")
            return 0
        elif self.process.returncode != 0 and self.process.returncode != 3:
            print(f"[警告] QE 进程未正常执行，返回值为 {self.process.returncode}")
            return 0

        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()

            # 如果 returncode == 3，从下往上查找 "energy new ="
            if self.process.returncode == 3:
                for line in reversed(lines):
                    if "     energy             new  =" in line:
                        try:
                            energy_str = line.split('=')[-1].strip().split()[0]
                            energy = float(energy_str)
                            return energy
                        except ValueError:
                            return 0  # 解析失败时返回 0
                return 0  # 未找到时返回 0

            # 正常查找 "Final energy" (适用于 returncode == 0)
            for line in reversed(lines):
                if "Final energy" in line:
                    try:
                        energy_str = line.split('=')[-1].strip().split()[0]
                        energy = float(energy_str)
                        return energy
                    except ValueError:
                        return 0

            return 0  # 仍未找到能量，返回 0

        except Exception as e:
            print(f"[错误] 读取 {output_file} 失败: {e}")
            return 0


    def parse_qe_results(self, output_file):
        result_file = os.path.join(self.folder_name, "qe.qe")
        energy = self.get_energy()
        positions = []
        extracting_positions = False
        
        try:
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # 反向读取文件，找到 ATOMIC_POSITIONS (angstrom)
            for i in range(len(lines) - 1, -1, -1):
                line = lines[i].strip()
                if "ATOMIC_POSITIONS (angstrom)" in line:
                    extracting_positions = True
                    start_index = i + 1
                    break
            
            if not extracting_positions:
                print("[警告] 未找到 ATOMIC_POSITIONS (angstrom)")
                return None
            
            # 从找到的位置开始向下读取坐标
            for line in lines[start_index:]:
                line = line.strip()
                
                if line == "" or "End final coordinates" in line:
                    break
                
                data = line.split()
                if len(data) >= 4:
                    element = data[0]
                    try:
                        x, y, z = map(float, data[1:4])  # 将坐标转为浮点数
                        positions.append((element, x, y, z))
                    except ValueError:
                        continue  # 如果无法转换为浮点数，则跳过该行
            
            # 如果找到有效的坐标，则保存到文件
            if positions:
                with open(result_file, 'w') as f:
                    f.write(f"# Final energy: {energy:.10f} Ry\n")
                    f.write(f"# Atomic positions (angstrom):\n")
                    for pos in positions:
                        f.write(f"{pos[0]}  {pos[1]:.10f}  {pos[2]:.10f}  {pos[3]:.10f}\n")
                return result_file
            else:
                print("[警告] 未找到有效的坐标")
                return None
        
        except Exception as e:
            print(f"[错误] 解析 {output_file} 时发生错误: {e}")
            return None

    def get_out_stru_name(self):
        """返回当前迭代的输出文件路径（.out文件）"""
        return os.path.join(self.folder_name, "qe.out")
    


class Dmol3(Abinitio):
    optimize_name = 'dmol3'
    write_stru = 'write_car' #写输入结构函数名

    def __init__(self, cluster, folder_name=None, **kw):
        '''必选参数：folder_name:在哪个文件夹里计算，不指定则用临时文件夹
        可选参数：task_name：文件夹中car和input文件的名字
        input_name：提供的input文件的名字
        fix_coord：优化时坐标的固定情况。列表中每个值代表对应一个原子，每位为1表示固定，例如5表示xz固定'''
        super(Dmol3,self).__init__(cluster, folder_name)
        if 'task_name' in kw:
            self.task_name = kw['task_name']
        else:
            self.task_name = 'dmol'
        if 'input_name' in kw:
            self.input_name = kw['input_name']
        else:
            self.input_name = 'optimize'
        self.fix_coord = kw['fix_coord'] if 'fix_coord' in kw else None

    @classmethod
    def read_cmd(cls):
        set_accelrys_path()
        global _accelrys_path
        cls.abinitio_cmd = _accelrys_path + '/etc/DMol3/bin/RunDMol3'
        if platform.system() == 'Windows':
            cls.abinitio_cmd += '.bat'
        else:
            cls.abinitio_cmd += '.sh'
        if not os.path.isfile(cls.abinitio_cmd):
            print('Warning, dmol path is error:', cls.abinitio_cmd)
        cls.abinitio_cmd = '"' + cls.abinitio_cmd + '"'

    def get_in_stru_name(self):
        return os.path.join(self.folder_name,self.task_name+'.car')

    def get_out_stru_name(self):
        return os.path.join(self.folder_name, self.task_name+'.car')

    def get_out_file_name(self):
        return os.path.join(self.folder_name, self.task_name+'.outmol')

    def prepare_input(self):
        if not os.path.isfile(self.input_name+'.input'):
            if self.input_name == 'energy' and os.path.isfile('optimize.input'):
                self.optimize2energy()
            else:
                print('no input file.')
                raise
        copyfile(self.input_name+'.input', os.path.join(self.folder_name, self.task_name+'.input'))

    def call_abinitio(self):
        command = '%s -np %d %s' %(self.abinitio_cmd, self.core_num, self.task_name)
        if platform.system() == 'Windows':
            self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name)
        elif platform.system() == 'Linux':
            self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid)
        else:
            print('do not support this operating system.')
            raise

    def write_fix_coord(self):
        '''根据优化的input文件写优化但固定原子的input文件
        self.fix_coord: 优化时不固定的原子列表'''
        optimize_input = open('optimize.input').read()
        #关闭对称性
        symmetry_on = re.compile('Symmetry +on').findall(optimize_input)
        if symmetry_on != []:
            optimize_input = optimize_input.replace(symmetry_on[0],'Symmetry                      off')
        f = open(os.path.join(self.folder_name,'%s.input'%self.task_name), 'w')
        f.write(optimize_input)
        f.write('\n#Cartesian constraints\nOpt_fixed\n')
        for i in range(len(self.fix_coord)):
            if self.fix_coord[i] == 0:
                continue
            f.write(str(i+1)+' ')
            if self.fix_coord[i]&4:
                f.write('X')
            if self.fix_coord[i]&2:
                f.write('Y')
            if self.fix_coord[i]&1:
                f.write('Z')
            f.write('\n')
        f.close()

    @staticmethod
    def read_energy(file_name):
        '''从outmol文件中提取能量
        参数file_name: 字符串类型,要读取的outmol文件名
        返回: 浮点型,读取到的能量,若为0则出错.'''
        energy = 0.
        line, f = utils.rfind_file(file_name, 'opt==')
        if line is None: #单点能
            line, f = utils.rfind_file(file_name, 'Ef')
            if line is None:
                return 0.
            else:
                energy = float(line[2:22])
        else: #几何优化
            next_line = f.readline()
            if next_line.startswith('opt=='): #只优化了一代
                line = f.readline()
                if line.startswith('opt=='):
                    energy = float(line[12:29])
            else:
                energy = float(line[12:29])

        while True:
            line = f.readline()
            if line == '':
                break
            if line.startswith('Message: DMol3 job finished successfully'):
                return energy
        #作业异常，返回0
        return 0.

    @staticmethod
    def homo_lumo(file_name):
        '''从outmol文件中读取HOMO_LUMO gap
        单位: eV
        适用于MS7.0'''
        line, f = utils.rfind_file(file_name, 'Energy of Highest Occupied Molecular Orbital')
        if line is None:
            return [0., 0., 0.]
        else:
            homo = float(line.split('eV')[0].split()[-1])
            next_line = f.readline()
            lumo = float(next_line.split('eV')[0].split()[-1])
            return homo, lumo, lumo-homo

    @staticmethod
    def force_coord_energy(file_name):
        '''从outmol文件中提取每代的力、坐标和能量'''
        result = []
        f = open(file_name)
        c = F = None
        while True:
            line = f.readline()
            if line == '':
                break
            if line.startswith('df') and line.find('DERIVATIVES') != -1: #读力
                f.readline()
                F = Cluster()
                try:
                    while True:
                        line = f.readline().split()
                        if line[1] == 'binding':
                            break
                        F.add_atom(Point(float(line[-3]), float(line[-2]), float(line[-1])), line[1])
                    F.zoom(0.5291772083) #au转埃
                except:
                    print('read force error in file',file_name)
                    F = None
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
                if line[1].isdigit():
                    if c is not None and F is not None:
                        c.set_energy(float(line[2]))
                        result.append((c,F))
                    c = F = None
        return result

    @staticmethod
    def bond_order(c, file_name):
        '''从Dmol的outmol文件中读取键级,返回矩阵'''
        f = open(file_name)
        while True:
            line = f.readline()
            if line == '' or line.find('Mayer bond orders') != -1:
                break
        if line.find('Mayer bond orders') == -1:
            print(file_name,'error.')
        bond_orders = np.zeros((c.get_size(),c.get_size()))
        while True:
            line = f.readline().split()
            if len(line) != 5:
                break
            bond_orders[int(line[1])-1, int(line[3])-1] = float(line[4])
            bond_orders[int(line[3])-1, int(line[1])-1] = float(line[4])
        return bond_orders

    @staticmethod
    def optimize2energy(optimize_name='optimize.input',energy_name='energy.input'):
        '''将几何优化的input转换为单点能的input'''
        energy_input = open(optimize_name).read()
        calculate = re.compile('Calculate +optimize').findall(energy_input)
        if calculate != []:
            energy_input = energy_input.replace(calculate[0],'Calculate                     energy')
        f = open(energy_name, 'w')
        f.write(energy_input)
        f.close()

    @staticmethod
    def energy2optimize(optimize_name='optimize.input',energy_name='energy.input'):
        '''将单点能的input转换为几何优化的input'''
        optimize_input = open(energy_name).read()
        energy = re.compile('Calculate +energy').findall(optimize_input)
        if energy != []:
            optimize_input = optimize_input.replace(energy[0],'Calculate                     optimize')
        f = open(optimize_name, 'w')
        f.write(optimize_input)
        f.close()

    @staticmethod
    def check_symmetry_on(input_name='optimize.input'):
        '''检查input文件中的对称性是否打开'''
        return re.compile('Symmetry +on').findall(open(input_name).read()) != []


class Vasp(Abinitio):
    optimize_name = 'vasp'
    write_stru = 'write_poscar'
    wait_run_time = 1800

    @classmethod
    def read_cmd(cls):
        cls.abinitio_cmd = cls.cfg.get('abinitio', 'vasp_path')

    def __init__(self, cluster, folder_name=None, **kw):
        '''必选参数：folder_name:在哪个文件夹里计算，不指定则用临时文件夹
        可选参数：task_name：文件夹中car和input文件的名字
        input_name：提供的input文件的名字
        fix_coord：优化时坐标的固定情况。列表中每个值代表对应一个原子，每位为1表示固定，例如5表示xz固定'''
        super(Vasp,self).__init__(cluster, folder_name)
        self.fix_coord = kw['fix_coord'] if 'fix_coord' in kw else None

    def get_in_stru_name(self):
        return self.folder_name

    def get_out_stru_name(self):
        return os.path.join(self.folder_name, 'CONTCAR')

    def get_out_file_name(self):
        return os.path.join(self.folder_name, 'OUTCAR')

    def prepare_input(self):
        '''复制输入文件INCAR、KPOINTS、POTCAR'''
        if not os.path.isfile('INCAR'):
            print('no INCAR file.')
            return
        copyfile('INCAR', os.path.join(self.folder_name, 'INCAR'))
        if not os.path.isfile('KPOINTS'):
            print('no KPOINTS file.')
            return
        copyfile('KPOINTS', os.path.join(self.folder_name, 'KPOINTS'))
        if not os.path.isfile('POTCAR'):
            print('no POTCAR file.')
            return
        copyfile('POTCAR', os.path.join(self.folder_name, 'POTCAR'))

    def call_abinitio(self):
        #调用vasp
        self.process = subprocess.Popen('mpirun -np %d %s > vasp.out' %(self.core_num, self.abinitio_cmd), shell=True, cwd=self.folder_name, stdout=subprocess.PIPE, close_fds=True, preexec_fn=os.setsid)

    def write_fix_coord(self):
        '''在POSCAR文件中写优化时要固定的坐标，F为固定'''
        self.modify_incar(os.path.join(self.folder_name, 'INCAR'), 'ISYM', 0) #不用对称性
        lines = open(os.path.join(self.folder_name, 'POSCAR')).readlines()
        f = open(os.path.join(self.folder_name, 'POSCAR'), 'w')
        for i in range(7):
            f.write(lines[i])
        f.write('Selective dynamics\n')
        if lines[7].find('Selective') == -1:
            f.write(lines[7])
            lines = lines[8:]
        else:
            f.write(lines[8])
            lines = lines[9:]
        for i in range(len(self.fix_coord)):
            f.write('\t'.join(lines[i].split()[:3]))
            if self.fix_coord[i] & 4:
                f.write('\tF')
            else:
                f.write('\tT')
            if self.fix_coord[i] & 2:
                f.write('\tF')
            else:
                f.write('\tT')
            if self.fix_coord[i] & 1:
                f.write('\tF')
            else:
                f.write('\tT')
            f.write('\n')

    @staticmethod
    def read_energy(file_name=''):
        '''从outcar文件中提取能量
        参数file_name: 字符串类型,要读取的outcar文件夹名
        返回: 浮点型,读取到的能量,若为0则出错.'''
        line, _ = utils.rfind_file(file_name,  '  energy  without entropy')
        if line is None:
            return 0.
        else:
            return float(line[26:43])

    @staticmethod
    def modify_incar(file_name, key, value):
        '''key为字符串,value不确定'''
        incar = open(file_name).read()
        key_split = re.compile('%s +='%key).split(incar,1)
        if len(key_split) == 1:
            print('INCAR file no %s'%key)
            return
        incar = key_split[0] + '%s = '%key+str(value)+'\n' + key_split[1].split('\n',1)[1]
        f = open(file_name, 'w')
        f.write(incar)
        f.close()

    @staticmethod
    def get_incar_value(file_name, key):
        '''key为字符串,返回value'''
        incar = open(file_name).read()
        key_split = re.compile('%s +='%key).split(incar,1)
        if len(key_split) == 1:
            return None
        return key_split[1].split('\n',1)[0].strip()

    @staticmethod
    def optimize2energy(file_name):
        '''将IBRION的值改为-1'''
        Vasp.modify_incar(file_name, 'IBRION', -1)

    @staticmethod
    def check_symmetry_on(file_name='INCAR'):
        '''检查INCAR文件中的对称性是否打开'''
        return Vasp.get_incar_value(file_name, 'ISYM') == '2'


class Gaussian(Abinitio):
    optimize_name = 'gaussian'
    write_stru = 'write_gaussian_stru' #写输入结构函数名
    wait_run_time = 1800
    in_ext = '.com' #gjf或com

    @classmethod
    def read_cmd(cls):
        cls.abinitio_cmd = cls.cfg.get('abinitio', 'gaussian_path')

    @staticmethod
    def set_parameter(cls, params):
        for k,v in params.items():
            cls.write_stru_ext_param[k] = v

    def __init__(self, cluster, folder_name=None, task_name='temp', **kw):
        super(Gaussian,self).__init__(cluster, folder_name)
        self.task_name = task_name
        if 'core_num' not in kw:
            kw['core_num'] = self.core_num
        if 'charge' not in kw:
            kw['charge'] = 0
        if 'basis' not in kw: #注意：默认basis为tpssh/6-311+g(d)
            kw['basis'] = 'tpssh/6-311+g(d)'
        self.write_stru_ext_param = kw

    def get_in_stru_name(self):
        return os.path.join(self.folder_name, self.task_name+Gaussian.in_ext)

    def get_out_stru_name(self):
        return os.path.join(self.folder_name, self.task_name+'.LOG')

    def get_out_file_name(self):
        return os.path.join(self.folder_name, self.task_name+'.LOG')

    def call_abinitio(self):
        command = '%s %s%s %s.LOG' %(self.abinitio_cmd, self.task_name, self.in_ext, self.task_name)
        self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid)

    @staticmethod
    def read_energy(file_name):
        '''从gaussian的输出文件.log中提取能量
        参数file_name: 字符串类型,要读取的log文件
        返回: 浮点型,读取到的能量,若为0则出错.'''
        if not os.path.isfile(file_name):
            print('no such file', file_name)
            return
        step = 10000
        f = open(file_name)
        f.seek(0, 2)
        end_pos = f.tell()

        while True: #分段往前读
            start_pos = max(end_pos-step, 0)
            f.seek(start_pos, 0)
            data = f.read(step+80).split('\n')
            data = ''.join([line.lstrip() for line in data])
            if data.find('\CCSD(T)=') != -1:
                energy_line = data.split('\CCSD(T)=')
            else:
                energy_line = data.split('\HF=')
            if len(energy_line) > 1:
                return float(energy_line[-1].split('\\')[0])
            elif data.find('SCF Done') != -1:
                return float(data.split('SCF Done')[1].split('=')[1].split()[0])
            if start_pos == 0:
                break
            else:
                end_pos = start_pos

        f.close()
        line, f = utils.rfind_file(file_name, ' SCF Done')
        if line is None:
            return 0.
        else:
            return float(line.split('=')[1].split()[0])

    @staticmethod
    def zpe(file_name=''):
        '''从gaussian的输出文件.log中提取零点能修正
        参数file_name: 字符串类型,要读取的log文件
        返回: 浮点型,读取到的能量,若为0则出错.'''
        if not os.path.isfile(file_name):
            print('no such file', file_name)
            return
        f = open(file_name)
        for _ in range(100):
            line = f.readline()
            if line.find('freq') != -1:
                break
        else:
            return None

        line, f = utils.rfind_file(file_name, ' Zero-point correction')
        zpc = float(line.split()[-2])
        for _ in range(10):
            line = f.readline()
            if line.startswith(' Sum of electronic and zero-point Energies'):
                return zpc, float(line.split()[-1])

        return None

    @staticmethod
    def homo_lumo(file_name):
        '''从gaussian的输出文件.log中提取HOMO,LUMO和gap
        返回: 浮点型,读取到的能量,若为0则出错.'''
        homo = lumo = 0.
        if not os.path.isfile(file_name):
            print('no such file', file_name)
            return
        step = 10000
        f = open(file_name)
        f.seek(0, 2)
        bottom_pos = f.tell()

        while True: #分段往前读
            f.seek(-min(step+80, bottom_pos), 1)
            while f.tell() < bottom_pos:
                line = f.readline()
                if line.startswith(' Alpha  occ. eigenvalues'):
                    homo = float(line[28:].split()[-1])
                    break
            else:
                if f.tell() < step:
                    break
                else:
                    f.seek(-step, 1)
                    bottom_pos = f.tell()
                    continue

            while True:
                line = f.readline()
                if line.startswith(' Alpha  occ. eigenvalues'):
                    homo = float(line[28:].split()[-1])
                elif line.startswith(' Alpha virt. eigenvalues'):
                    lumo = float(line[28:].split()[0])
                    if homo == 0.:
                        print('gauss_homo_lumo error.')
                        break
                    return homo,lumo,lumo-homo
                else:
                    print('gauss_homo_lumo error.')
                    break

            if homo != 0.:
                if line.startswith(' Alpha virt. eigenvalues'):
                    lumo = float(line[28:].split()[0])
                    return homo,lumo,lumo-homo
                else:
                    print('extract homo lumo error.')
                    return [0.,0.,0.]

        return [0.,0.,0.]

    @staticmethod
    def bond_order(c, file_name):
        '''从Gaussian的LOG文件中读取键级,返回矩阵'''
        f = open(file_name)
        while True:
            line = f.readline()
            if line == '':
                break
            if not line.startswith(' Wiberg bond index matrix in the NAO basis'):
                continue
            bond_orders = np.zeros((c.get_size(),c.get_size()))
            n_column = 0
            while True:
                line = f.readline()
                line = f.readline()
                line = f.readline()
                if not line.startswith('     ----'):
                    break
                for i in range(c.get_size()):
                    line = f.readline().split()
                    for j in range(len(line)-2):
                        bond_orders[i,j+n_column] = float(line[j+2])
                n_column += len(line)-2
        return bond_orders


class Abacus(Abinitio):
    optimize_name = 'abacus'
    write_stru = 'write_abacus_stru' #写输入结构函数名

    @classmethod
    def read_cmd(cls):
        cls.abinitio_cmd = cls.cfg.get('abinitio', 'abacus_path')

    def __init__(self, cluster, folder_name=None, **kw):
        super(Abacus,self).__init__(cluster, folder_name)
        #优化时不固定的原子列表
        self.relax_list = kw['relax_list'] if 'relax_list' in kw else None
        self.write_stru_ext_param = self.relax_list

    def get_in_stru_name(self):
        return self.folder_name

    def get_out_stru_name(self):
        return os.path.join(self.folder_name,'OUT.ABACUS/STRU_ION_D')

    def get_out_file_name(self):
        return os.path.join(self.folder_name, 'OUT.ABACUS/running_relax.log')

    def prepare_input(self):
        #复制input文件
        copyfile('INPUT', os.path.join(self.folder_name,'INPUT'))
        copyfile('KPT', os.path.join(self.folder_name,'KPT'))
        f = open('STRU')
        while True:
            line = f.readline()
            if line == '':
                break
            if line.startswith('ATOMIC_SPECIES'):
                upf_name = f.readline().split()[2]
                copyfile(upf_name, os.path.join(self.folder_name, upf_name))
            if line.startswith('NUMERICAL_ORBITAL'):
                orbital_name = f.readline().split()[0]
                copyfile(orbital_name, os.path.join(self.folder_name, orbital_name))

    def call_abinitio(self):
        command = 'mpiexec -np %d %s > abacus_log.txt' %(self.core_num, self.abinitio_cmd)
        self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid)

    @staticmethod
    def read_energy(file_name):
        '''从abacus的输出文件running_relax.log中提取能量
        参数file_name: 字符串类型,要读取的running_relax.log所在文件夹名,默认当前文件夹
        返回: 浮点型,读取到的能量,若为0则出错.'''
        line, _ = utils.rfind_file(file_name, ' !FINAL_ETOT_IS')
        if line is None:
            return 0.
        else:
            return float(line.split('FINAL_ETOT_IS')[1].split('eV')[0])


class Lammps(Abinitio):
    optimize_name = 'lammps'
    write_stru = 'write_lammps_stru' #写输入结构函数名
    wait_run_time = 5
    wait_start_time = 1
    wait_stru_time = 1
    wait_energy_time = 1

    @classmethod
    def read_cmd(cls):
        if platform.system() == 'Windows':
            cls.abinitio_cmd = cls.cfg.get('abinitio', 'lammps_windows')
        else:
            cls.abinitio_cmd = cls.cfg.get('abinitio', 'lammps_linux')

    def __init__(self, cluster, folder_name=None, **kw):
        super(Lammps,self).__init__(cluster, folder_name)
        if 'input_name' in kw:
            self.input_name = kw['input_name']
        else:
            self.input_name = 'optimize'

    def get_in_stru_name(self):
        return os.path.join(self.folder_name,'pre_opt.dat')

    def get_out_stru_name(self):
        return os.path.join(self.folder_name,'dump.dat')

    def get_out_file_name(self):
        return os.path.join(self.folder_name,'log.lammps')

    def prepare_input(self):
        copyfile(self.input_name+'.input', os.path.join(self.folder_name,'optimize.input'))

    def call_abinitio(self):
        command = 'mpirun -np %d %s -in %s.input > log.lammps' %(self.core_num, self.abinitio_cmd, self.input_name)
        if platform.system() == 'Windows':
            command = '%s -in optimize.input > log.lammps' %self.abinitio_cmd
            self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name)
        elif platform.system() == 'Linux':
            self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid)
        else:
            print('do not support this operating system.')
            raise

    @staticmethod
    def read_energy(file_name):
        _, f = utils.rfind_file(file_name, contain_word='Energy initial')
        return float(f.readline().split()[2])


class Clean(Abinitio):
    optimize_name = 'clean'
    write_stru = 'write_xsd' #写输入结构函数名

    @classmethod
    def read_cmd(cls):
        set_accelrys_path()
        global _accelrys_path
        global _ms_version
        if _ms_version >= 70:
            cls.abinitio_cmd = _accelrys_path + r'/etc/Scripting/bin/RunMatScript'
        else:
            cls.abinitio_cmd = _accelrys_path + r'/etc/MatServer/bin/RunMatServer'
        if platform.system() == 'Windows':
            cls.abinitio_cmd += '.bat'
        else:
            cls.abinitio_cmd += '.sh'
        if not os.path.isfile(cls.abinitio_cmd):
            print('Warning, RunMatServer path is error:', cls.abinitio_cmd)
        cls.abinitio_cmd = '"' + cls.abinitio_cmd + '"'
        return cls.abinitio_cmd

    def __init__(self, cluster, folder_name=None, **kw):
        '''额外参数：task_name,bonds'''
        super(Clean,self).__init__(cluster, folder_name)
        self.bonds = kw['bonds'] if 'bonds' in kw else self.cluster.bond_adj()
        self.task_name = kw['task_name'] if 'task_name' in kw else 'clean'
        self.write_stru_ext_param = self.bonds

    def get_in_stru_name(self):
        return os.path.join(self.folder_name,self.task_name+'.xsd')

    def get_out_stru_name(self):
        global _ms_version
        if _ms_version >= 70:
            return os.path.join(self.folder_name, 'clean_Files/Documents/%s.xsd' %self.task_name)
        else:
            return os.path.join(self.folder_name, 'clean Files/Documents/%s.xsd' %self.task_name)

    def call_abinitio(self):
        #写入perl脚本
        f = open(os.path.join(self.folder_name, self.task_name+'.pl'), 'w')
        f.write('''#!perl
use strict;
use MaterialsScript qw(:all);
Documents->Import("%s.xsd");
my $doc = $Documents{"%s.xsd"};
$doc->Clean();''' %(self.task_name, self.task_name))
        f.close()
        #运行脚本
        command = self.abinitio_cmd + ' ' + self.task_name
        self.process = subprocess.call(command, shell=True, cwd=self.folder_name, stdout=subprocess.PIPE)

    def get_energy(self):
        pass

    def _wait(self):
        pass


class Symmetry(Abinitio):
    optimize_name = 'symmetry'
    write_stru = 'write_xsd' #写输入结构函数名

    @classmethod
    def read_cmd(cls):
        '''同clean'''
        cls.abinitio_cmd = Clean.read_cmd()

    def __init__(self, cluster, folder_name=None, **kw):
        '''额外参数：task_name,tolerance'''
        super(Symmetry,self).__init__(cluster, folder_name)
        self.tolerance = kw['tolerance'] if 'tolerance' in kw else 0.1
        self.task_name = kw['task_name'] if 'task_name' in kw else 'sym'

    def get_in_stru_name(self):
        return os.path.join(self.folder_name,self.task_name+'.xsd')

    def get_out_stru_name(self):
        return os.path.join('%s/%s_Files/Documents/%s.xsd' %(self.folder_name,self.task_name,self.task_name))

    def get_out_file_name(self):
        global _ms_version
        if _ms_version >= 70:
            return os.path.join(self.folder_name, '%s.pl.out' %self.task_name)
        else:
            return os.path.join(self.folder_name, 'stdout.txt')

    def call_abinitio(self):
        #写入perl脚本
        f = open(os.path.join(self.folder_name, self.task_name+'.pl'), 'w')
        f.write('''#!perl
use strict;
use MaterialsScript qw(:all);
Documents->Import("%s.xsd");
my $doc = $Documents{"%s.xsd"};
Tools->Symmetry->ChangeSettings([PositionTolerance => %f]);
my $result = Tools->Symmetry->FindSymmetry->Find($doc);
printf "%%s", $result->PointGroupSchoenfliesName;
''' %(self.task_name, self.task_name, self.tolerance))
        f.close()
        command = self.abinitio_cmd + ' ' + self.task_name
        self.process = subprocess.call(command, shell=True, cwd=self.folder_name, stdout=subprocess.PIPE)

    def _wait(self):
        pass

    def get_energy(self):
        '''获取对称性，保存在self.sym'''
        #可能传输结果需要时间，这里等待10s
        for _ in range(1000):
            if os.path.isfile(self.get_out_file_name()):
                break
            time.sleep(0.01)
        else:
            utils.write_log('can not find energy file.')
            return
        self.sym = open(self.get_out_file_name()).read().rstrip()


class Castep(Abinitio):
    optimize_name = 'castep'
    write_stru = 'write_cell' #写输入结构函数名
    wait_run_time = 1800
    in_ext = '.cell'

    @classmethod
    def read_cmd(cls):
        set_accelrys_path()
        global _accelrys_path
        cls.abinitio_cmd = _accelrys_path + '/etc/CASTEP/bin/RunCASTEP'
        if platform.system() == 'Windows':
            cls.abinitio_cmd += '.bat'
        else:
            cls.abinitio_cmd += '.sh'
        if not os.path.isfile(cls.abinitio_cmd):
            print('Warning, dmol path is error:', cls.abinitio_cmd)
        cls.abinitio_cmd = '"' + cls.abinitio_cmd + '"'

    def __init__(self, cluster, folder_name=None, **kw):
        super(Castep,self).__init__(cluster, folder_name)
        if 'task_name' in kw:
            self.task_name = kw['task_name']
        else:
            self.task_name = 'castep'
        if 'input_name' in kw:
            self.input_name = kw['input_name']
        else:
            self.input_name = 'castep'

    def get_in_stru_name(self):
        return os.path.join(self.folder_name, self.task_name+Castep.in_ext)

    def get_out_stru_name(self):
        return os.path.join(self.folder_name, self.task_name+'.castep')

    def get_out_file_name(self):
        return os.path.join(self.folder_name, self.task_name+'.castep')

    def prepare_input(self):
        if not os.path.isfile(self.input_name+'.param'):
            copyfile(self.input_name+'.param', os.path.join(self.folder_name, self.task_name+'.param'))

    def call_abinitio(self):
        command = '%s %s%s %s.LOG' %(self.abinitio_cmd, self.task_name, self.in_ext, self.task_name)
        self.process = subprocess.Popen(command, shell=True, cwd=self.folder_name, close_fds=True, preexec_fn=os.setsid)

    @staticmethod
    def read_energy(file_name):
        '''从castep的输出文件.castep中提取能量
        参数file_name: 字符串类型,要读取的castep文件
        返回: 浮点型,读取到的能量,若为0则出错.'''
        line, _ = utils.rfind_file(file_name, start_word='Final energy')
        if line is None:
            return 0.
        else:
            return float(line[15:32])


class Dmol3Pbs(Dmol3):
    '''未调试'''
    optimize_name = 'dmol_pbs'

    @classmethod
    def read_cmd(cls):
        cls.abinitio_cmd = cls.cfg.get('abinitio', 'dmol3_pbs')

    def __init__(self, **kw):
        Dmol3.__init__(self, **kw)

    def call_abinitio(self):
        self.process = subprocess.Popen(Dmol3Pbs.abinitio_cmd, shell=True, cwd=self.folder_name,
        stdout=subprocess.PIPE, close_fds=True, preexec_fn=os.setsid)
        print('submit print is:\n',self.process.stdout)


    def _wait(self):
        '''等待第一性原理计算结束'''
        file_name = self.get_out_file_name()
        if not file_name:
            print(f"[WARNING] 输出文件路径为 None，跳过等待步骤。")
            return
        # 监控文件生成
        for _ in range(int(self.wait_start_time * 100)):
            if os.path.isfile(file_name):
                break
            time.sleep(0.01)
        else:
            print(f"[ERROR] 输出文件 {file_name} 未生成，可能任务启动失败。")
            return

    @staticmethod
    def get_free_core():
        '''获取集群空闲的核数'''
        p = subprocess.Popen('qstat -an | grep node', shell=True, stdout=subprocess.PIPE)
        out = p.stdout.readlines()
        print(out)

    @staticmethod
    def get_pbs_no(s):
        '''根据subprocess.Popen返回的结果字符串获取pbs作业号'''
        return s[0].split()[2] #适用于xu

    @staticmethod
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
        
        
class Empirical_Potentail:
    '''经验势优化基类'''
    optimize_name = 'ep'
        
    def __init__(self, c, n_iter=1000, **kw):
        '''c：要优化的团簇
        n_iter:迭代次数'''
        self.cluster = c
        self.n_iter = n_iter

    @classmethod
    def energy_func(cls):
        pass
        
    @classmethod
    def force_func(cls):
        pass
        
    def optim(self, r):
        '''用lbfgs优化'''
        def fp(x):
            return -self.force_func(x.reshape(-1,3)).flatten()
        
        return opt.fmin_bfgs(self.energy_func, r, fp, maxiter=self.n_iter, disp=0).reshape(-1,3)

    def __call__(self):
        r = self.cluster.coord_matrix()
        r = self.optim(r)
        for i in range(self.cluster.get_size()):
            self.cluster.atoms[i].x = r[i,0]
            self.cluster.atoms[i].y = r[i,1]
            self.cluster.atoms[i].z = r[i,2]
        if hasattr(self.cluster, 'after_read'):
            self.cluster.after_read()
        self.cluster.set_energy(self.energy_func(r))
        
        
class LJ_Potentail(Empirical_Potentail):
    '''Lennard-Jones势优化'''
    optimize_name = 'lj'
        
    def __init__(self, c, n_iter=1000, **kw):
        super(LJ_Potentail,self).__init__(c, n_iter)

    @classmethod
    def energy_func(cls, coord, dist=None):
        '''V=4ε[(σ/r)^(-12) - (σ/r)^(-6)]
        coord：坐标。维数(n,3)或(n*3,)
        dist：距离矩阵(n*n)'''    
        epsilon = 2
        coord = coord.reshape(-1,3)
        if dist is None:
            r = coord[np.newaxis,:,:] - coord[:,np.newaxis,:]
            dist = np.linalg.norm(r, axis=-1)
        e = np.power(dist, -12) - np.power(dist,-6)
        e[np.diag_indices_from(e)] = 0
        return epsilon*np.sum(e)
        
    @classmethod
    def force_func(cls, coord, dist=None):
        '''coord：坐标(n*3)
        dist：距离矩阵(n*n)'''
        epsilon = 24
        r = coord[np.newaxis,:,:] -coord[:,np.newaxis,:]
        if dist is None:
            dist = np.linalg.norm(r, axis=-1)
        fc = -2*np.power(dist, -14) + np.power(dist,-8)
        fc[np.diag_indices_from(fc)] = 0
        return epsilon * (fc[:,:,np.newaxis]*r).sum(axis=1)
        
        
class GCN_Potentail(Empirical_Potentail):
    '''图卷积拟合势优化'''
    optimize_name = 'gcn'
    
    if 'torch' in sys.modules:
        device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
        _path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
        net_e = GCN(36,28,80,120,72)
        net_e.load_state_dict(torch.load(os.path.join(_path, 'gcn/energy.pth'), 
                              map_location=lambda storage, loc: storage)) #能量预测网络
        net_e.eval()
        net_f = GCNF(36,28,80,120)
        net_f.load_state_dict(torch.load(os.path.join(_path, 'gcn/force.pth'), 
                              map_location=lambda storage, loc: storage)) #力预测网络
        net_f.eval()
        
    def __init__(self, c, n_iter=1000, **kw):
        super(GCN_Potentail, self).__init__(c, n_iter)

    @classmethod
    def energy_func(cls, x):
        x = torch.from_numpy(x).view(1,-1,3).float()
        r = x.unsqueeze(1) - x.unsqueeze(2)
        d = torch.norm(r, dim=-1)
        energy = GCN_Potentail.net_e(d).item()
        return -energy #预测的是平均结合能，返回负值才能计算最小值
        
    @classmethod
    def force_func(cls, x):
        '''x：坐标矩阵。维数(n,3)或(n*3,)'''
        x = torch.from_numpy(x).view(1,-1,3).float()
        r = x.unsqueeze(1) - x.unsqueeze(2)
        d = torch.norm(r, dim=-1)
        return GCN_Potentail.net_f(x, d).detach().numpy()
    
    def optim(self, r):
        def fp(x):
            return self.force_func(x).flatten() #负负得正
        
        return opt.fmin_bfgs(self.energy_func, r, fprime=fp, 
                             gtol=1e-4,  
                             maxiter=500, disp=0).reshape(-1,3)

###############################################################################
#                 第一性原理计算类 定义结束
###############################################################################

_optimize_type = Abinitio


def set_optimize_method(method):
    method_list = {
        'qe': QE,
        'dmol3': Dmol3,
        'vasp': Vasp,
        'gaussian': Gaussian,
        'abacus': Abacus,
        'lammps': Lammps,
        'clean': Clean,
        'symmetry': Symmetry,
        'dmol3_pbs': Dmol3Pbs,
        'castep': Castep,
        'lj': LJ_Potentail,
        'gcn': GCN_Potentail
    }
    # 检查方法是否存在
    if method not in method_list:
        raise ValueError(f"[ERROR] 未知的优化方法: {method}")

    # 检查 abinitio_cmd 是否正确设置
    if not hasattr(method_list[method], 'abinitio_cmd') or not method_list[method].abinitio_cmd:
        method_list[method].read_cmd()  # 调用 read_cmd 方法重新加载

    # 检查命令路径是否正确
    cmd = method_list[method].abinitio_cmd
    if cmd[0] == '"':
        cmd = cmd[1:-1]
    if not os.path.isfile(cmd):
        raise FileNotFoundError(f"[ERROR] 命令 {cmd} 不存在，请检查配置或安装情况")

    global _optimize_type
    _optimize_type = method_list[method]



def get_optimize_name():
    return _optimize_type.optimize_name


def get_optimize_method():
    return _optimize_type


def optimize(cluster, folder_name=None, **kw):
    """
    Main optimization function.
    """
    if folder_name is None:
        folder_name = os.path.abspath("./")
    
    # Generate a unique folder for this optimization
    base_name = f"qe_{''.join(random.choices(string.ascii_lowercase + string.digits, k=8))}"
    folder_name = os.path.join(folder_name, base_name)
    print(f"\n")
    
    # Call the optimization method with simulation mode
    optimize_method = _optimize_type(cluster, folder_name=folder_name, simulate=False, simulated_output='/home/cast/users/htj/cqc/qe.out', **kw)
    optimize_method()




    
Abinitio.set_path(os.path.abspath("./"))
Abinitio.config() #第一性原理计算配置
