import sys
sys.path.append('../..')
import os
from cmds.src.cluster import Cluster


class Generation:
    def __init__(self):
        self.method = None
        self.father = None
        self.mother = None
        self.pre_cluster = None
        self.opt_cluster = None
        self.replace = None
        self.gen = None
        
class Cga_Analysis:
    def __init__(self):
        self.pop = []
        self.iter_info = []
        
    def read_init(self, f):
        f.readline()
        while True:            
            line = f.readline()
            if not line.startswith('init'):
                break
            try:
                c = Cluster()
                line = c.read_file_obj(f)
            except:
                print('文件格式错误')
                raise
            else:
                c.gen = 0
                self.pop.append([c])
        return line
                
    def read_iter(self, f, line): 
        while line != '':
            g = Generation()
            try:
                assert line.startswith('gen')
                g.gen = int(line.split()[1])
                
                line = f.readline()
                assert line.startswith('method')
                g.method = line.split('=')[1]
                
                line = f.readline()
                assert line.startswith('father')
                g.father = line.split()[1]
                
                line = f.readline()
                if line.startswith('mother'):
                    g.mother = line.split()[1]
                    line = f.readline()
                    
                assert line.startswith('son')
                c = Cluster()
                line = c.read_file_obj(f)
                g.pre_cluster = c
                assert line.startswith('optimized')
                line = c.read_file_obj(f)
                c.gen = g.gen
                g.opt_cluster = c
                
                if line.startswith('replace'):
                    g.replace = int(line.split()[-1])
                    self.pop[g.replace].append(g.opt_cluster)
            except:
                while True:
                    line = f.readline()
                    if line == '':
                        return
                    elif line.startswith('gen'):
                        break
                continue
            else:
                self.iter_info.append(g)
                f.readline()
                line = f.readline()            
            
    def rebuild(self, folder_name):
        f = open(os.path.join(folder_name,'record.txt'))
        line = self.read_init(f)
        self.read_iter(f, line)

    def to_recover(self):
        f = open('recover.txt', 'w')
        for i,p in enumerate(self.pop):
            f.write('pop %d\n' %i)
            f.write(str(self.p[-1]))
            f.write('\n')
            
    def get_structure(self, idx, gen):
        '''idx：种群编号
        gen：代数'''
        result = None
        for c in self.pop[idx]:
            if c.gen <= gen:
                result = c
            else:
                break
        return result
    
    def get_pop(self, gen):
        '''这一代之前的种群'''
        return [self.get_structure(i, gen) for i in range(len(self.pop))]
    
    
if __name__ == '__main__':
    cgaa = Cga_Analysis()
    cgaa.rebuild(r'E:\cluster\B\search\56\B56_C2')
    c = cgaa.get_structure(4,3148)
    c.write_xyz('b.xyz')
    print(c.gen)