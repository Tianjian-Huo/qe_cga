import os
import shutil
import json
import re

element_ecut_map = {
    "H": 20, "He": 25, "Li": 10, "Be": 10, "B": 10, "C": 12, "N": 20, "O": 17.5, "F": 17.5, "Ne": 20,
    "Na": 17.5, "Mg": 10, "Al": 12, "Si": 10, "P": 15, "S": 15, "Cl": 17.5, "Ar": 20, "K": 15, "Ca": 12,
    "Sc": 10, "Ti": 12, "V": 10, "Cr": 25, "Mn": 12, "Fe": 12, "Co": 12, "Ni": 20, "Cu": 15, "Zn": 12,
    "Ga": 12, "Ge": 15, "As": 10, "Se": 10, "Br": 10, "Kr": 15, "Rb": 10, "Sr": 10, "Y": 10, "Zr": 12,
    "Nb": 10, "Mo": 10, "Tc": 10, "Ru": 10, "Rh": 10, "Pd": 17.5, "Ag": 12, "Cd": 12, "In": 12, "Sn": 12,
    "Sb": 10, "Te": 10, "I": 10, "Xe": 20, "Cs": 17.5, "Ba": 10, "Hf": 10, "Ta": 10, "W": 10, "Re": 10,
    "Os": 10, "Lr": 10, "Pt": 15, "Au": 15, "Hg": 15, "Pb": 10, "Bi": 10, "Po": 10, "At": -1, "Rn": 15
}

def setup_structure():
    structure_name = input("请输入原子结构名称，例如 'Ca1Ba2': ").strip()
    base_dir = os.getcwd()
    structure_path = os.path.join(base_dir, structure_name)
    search_path = os.path.join(structure_path, "search")

    os.makedirs(search_path, exist_ok=True)
    print(f"已创建目录: {structure_path} 和 {search_path}")

    cmds_dir = os.path.join(base_dir, "cmds")
    required_files = ["config.ini", "optimize.input"]

    for file in required_files:
        src = os.path.join(cmds_dir, file)
        dst = os.path.join(search_path, file)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"已复制文件: {file} 到 {search_path}")
        else:
            print(f"[WARNING] 源文件 {file} 不存在: {src}")

    config_path = os.path.join(search_path, "config.ini")
    modify_config_file(config_path, structure_name)
    copy_additional_files(structure_path, cmds_dir)
    create_qein_json(search_path, structure_name)

def modify_config_file(config_path, structure_name):
    element_counts = {}
    pattern = r"([A-Za-z]+)(\d*)"
    matches = re.findall(pattern, structure_name)

    for match in matches:
        element, num_str = match
        count = int(num_str) if num_str else 1
        element_counts[element] = count

    total_atoms = sum(element_counts.values())
    pop_size = 6 if total_atoms <= 6 else 16
    max_gen = 100 if total_atoms <= 6 else 250

    with open(config_path, "r") as file:
        lines = file.readlines()

    with open(config_path, "w") as file:
        for line in lines:
            if line.startswith("element"):
                file.write(f"element = {' '.join(element_counts.keys())}\n")
            elif line.startswith("atom_num"):
                file.write(f"atom_num = {' '.join(map(str, element_counts.values()))}\n")
            elif line.startswith("pop_size"):
                file.write(f"pop_size = {pop_size}\n")
            elif line.startswith("max_gen"):
                file.write(f"max_gen = {max_gen}\n")
            else:
                file.write(line)
    print(f"[SUCCESS] config.ini 已更新: element, atom_num, pop_size, max_gen")

def copy_additional_files(structure_path, cmds_dir):
    extra_files = ["cmds_cmd.py", "cga.pbs", "config.ini"]
    for file in extra_files:
        src = os.path.join(cmds_dir, file)
        dst = os.path.join(structure_path, file)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"已复制文件: {file} 到 {structure_path}")
        else:
            print(f"[WARNING] 源文件 {file} 不存在: {src}")

def create_qein_json(search_path, structure_name):
    element_counts = {}
    pattern = r"([A-Za-z]+)(\d*)"
    matches = re.findall(pattern, structure_name)
    
    for match in matches:
        element, num_str = match
        count = int(num_str) if num_str else 1
        element_counts[element] = count

    max_ecut = max([element_ecut_map.get(el, 20) for el in element_counts.keys()])
    ecutwfc = max_ecut * 2
    ecutrho = ecutwfc * 6

    qein_content = {
        "forc_conv_thr": 1e-3,
        "nstep": 30,
        "ecutwfc": ecutwfc,
        "ecutrho": ecutrho,
        "electron_maxstep": 30,
        "conv_thr": "1e-6",
        "degauss": 0.02,
        "smearing": "gaussian",
        "mixing_beta": 0.5,
        "CELL_PARAMETERS": "offset",
        "offset": 10.0,
    }
    
    qein_path = os.path.join(search_path, "qein.json")
    
    with open(qein_path, "w") as json_file:
        json.dump(qein_content, json_file, indent=4)

    print(f"[INFO] 创建 qein.json 文件并写入内容: {qein_path}")

if __name__ == "__main__":
    setup_structure()
