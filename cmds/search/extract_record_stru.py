import sys
sys.path.append('../src')
sys.path.append('../cga')
from cluster import Cluster
from cga_analysis import Cga_Analysis


folder_list = ['B56_C2','B56_C3','B56_planar','B56-C1-ok','B56-C2-2',
               'B56-C2-n','B56-C2-n-ok','B56-C3-2','B56-C3-n','B56-C3-n-ok','B56-Ci-n',
               'B56-Cs-1-ok','B56-Cs-2','B56-Cs-n']

total_stru = []
for folder in folder_list:
    cgaa = Cga_Analysis()
    cgaa.rebuild(r'E:\cluster\B\search\56'+'\\'+folder)
    for p in cgaa.pop:
        total_stru.append(p[0])
    for g in cgaa.iter_info:
        total_stru.append(g.opt_cluster)
    print('找到%d个结构'%len(cgaa.iter_info))
print('共%d个结构'%len(total_stru))
        
non_iso_list = []
for c in total_stru:
    c.calc_fingerprint()
    for c2 in non_iso_list:
        if c.isomorphic(c2):
            break
    else:
        non_iso_list.append(c)
print('共%d个不同构结构'%len(non_iso_list))
        
for c in non_iso_list:
    c.append_arc('B56.arc')