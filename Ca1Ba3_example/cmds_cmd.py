'''可以把该文件放到你任何位置
运行方式为：python cmds_cmd.py 参数1 参数2 参数3
参数1为要执行什么操作，候选值为cga,sort,extract
参数2为在哪个文件夹下运行
参数3为额外参数'''
import sys
sys.path.append('..')
from cmds import CGA,sort_pop,extract_result
import os

if len(sys.argv) == 1:
    print('请输入参数')
    sys.exit(0)
    
if sys.argv[1] == 'cga':
    if len(sys.argv) > 2:
        os.chdir(sys.argv[2])
    ga = CGA()
    ga()
elif sys.argv[1] == 'sort':
    if len(sys.argv) > 3:
        sort_pop(sys.argv[2], sys.argv[3]) #参数3为写结构的文件名前缀
    else:
        sort_pop(sys.argv[2]) #参数2为recover所在路径
elif sys.argv[1] == 'extract':
    extract_result(sys.argv[2]) #参数2为要提取的文件夹所在的文件夹