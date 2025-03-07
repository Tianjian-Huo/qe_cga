# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cga_ui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog,QTableWidgetItem
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from copy import deepcopy
from vasp_set import Ui_VASP_Set
from cga_set import Ui_CGA_Set
from dmol3_set import Ui_DMol3_set
import sys
sys.path.append('..')
from src.atom import get_color
from src.cluster import Cluster
from src.genetic_algorithm import GA


#通过继承FigureCanvas类，使得该类既是一个PyQt5的Qwidget，又是一个matplotlib的FigureCanvas
#这是连接pyqt5与matplotlib的关键
class EnergyFigure(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        #创建一个Figure，注意：该Figure为matplotlib下的figure，不是matplotlib.pyplot下面的figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.move(350,30)
        
    def plot(self):
        self.axes = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[0.1,0.1,1.,1.]) #否则左边的energy标签会被挡住
        
        lines = open('energy.txt').readlines()
        step = []
        energy = []
        for line in lines:
            col = line.split()
            if len(col) != 5:
                continue
            step.append(int(col[0][:-1]))
            energy.append(float(col[2]))
        self.axes.plot(step,energy,'b-',linewidth=1,ls='-')
        for s,e in zip(step,energy):
            self.axes.scatter(s,e,marker='o',color='r')
        self.axes.set_xlabel('step')
        self.axes.set_ylabel('energy')
        
        self.axes.set_title(u'energy curve')
        self.move(350,30)
        self.draw()
        
        
class ClusterFigure(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.move(870,30)
        
    def plot(self, c):      
        self.axes = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[0.1,0.1,1.,1.]) #否则左边的energy标签会被挡住
        
        adj = c.adjacent()
        ax = Axes3D(self.fig)
        ax.patch.set_facecolor('black') #背景色
        ax.set_axis_off() #不显示坐标轴
        
        for i in xrange(c.get_size()):
            ax.scatter(c.atoms[i].x, c.atoms[i].y, c.atoms[i].z, 
                       s=400, c=get_color(c.atoms[i].elem)) #画原子
            for j in xrange(i): #画键
                if adj[i,j] == 0:
                    continue
                x = [c.atoms[i].x, c.atoms[j].x]
                y = [c.atoms[i].y, c.atoms[j].y]
                z = [c.atoms[i].z, c.atoms[j].z]
                ax.plot(x, y, z, lw=4, c='#7f7f7f')
        
        #设置坐标轴显示范围
        left_bound = 0.8*min(min([p.x for p in c.atoms]),
                         min([p.y for p in c.atoms]),
                         min([p.z for p in c.atoms]))
        right_bound = 0.8*max(max([p.x for p in c.atoms]),
                         max([p.y for p in c.atoms]),
                         max([p.z for p in c.atoms]))
        ax.set_xlim(left_bound, right_bound)
        ax.set_ylim(left_bound, right_bound)
        ax.set_zlim(left_bound, right_bound)
        
        self.draw()


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1401, 789)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("logo.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pop_info = QtWidgets.QTableWidget(self.centralwidget) #改为QTableWidget
        self.pop_info.setGeometry(QtCore.QRect(10, 30, 320, 400))
        self.pop_info.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.pop_info.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.pop_info.setSortingEnabled(True)
        self.pop_info.setObjectName("pop_info")
        self.pop_info.horizontalHeader().setSortIndicatorShown(True)
        self.pop_info.verticalHeader().setVisible(False)
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 500, 291, 31))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(self.horizontalLayoutWidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.prefix = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.prefix.setObjectName("prefix")
        self.horizontalLayout.addWidget(self.prefix)
        self.save_2 = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.save_2.setObjectName("save_2")
        self.horizontalLayout.addWidget(self.save_2)
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(330, 440, 1000, 340))
        self.textBrowser.setObjectName("textBrowser")
        #MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 901, 23))
        self.menubar.setObjectName("menubar")
        self.menu_abinitio = QtWidgets.QMenu(self.menubar)
        self.menu_abinitio.setObjectName("menu_abinitio")
        self.menuCGA = QtWidgets.QMenu(self.menubar)
        self.menuCGA.setObjectName("menuCGA")
        self.menu_help = QtWidgets.QMenu(self.menubar)
        self.menu_help.setObjectName("menu_help")
        self.menu_file = QtWidgets.QMenu(self.menubar)
        self.menu_file.setObjectName("menu_file")
        self.menu = QtWidgets.QMenu(self.menubar)
        self.menu.setObjectName("menu")
        #MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        #MainWindow.setStatusBar(self.statusbar)
        self.dmol3_set = QtWidgets.QAction(MainWindow)
        self.dmol3_set.setCheckable(True)
        self.dmol3_set.setObjectName("dmol3_set")
        self.vasp_set = QtWidgets.QAction(MainWindow)
        self.vasp_set.setCheckable(True)
        self.vasp_set.setObjectName("vasp_set")
        self.gaussian_set = QtWidgets.QAction(MainWindow)
        self.gaussian_set.setCheckable(True)
        self.gaussian_set.setObjectName("gaussian_set")
        self.lammps_set = QtWidgets.QAction(MainWindow)
        self.lammps_set.setCheckable(True)
        self.lammps_set.setObjectName("lammps_set")
        self.abacus_set = QtWidgets.QAction(MainWindow)
        self.abacus_set.setCheckable(True)
        self.abacus_set.setObjectName("abacus_set")
        self.new_cga = QtWidgets.QAction(MainWindow)
        self.new_cga.setObjectName("new_cga")
        self.continue_cga = QtWidgets.QAction(MainWindow)
        self.continue_cga.setObjectName("continue_cga")
        self.action_CGA = QtWidgets.QAction(MainWindow)
        self.action_CGA.setObjectName("action_CGA")
        self.open_folder = QtWidgets.QAction(MainWindow)
        self.open_folder.setObjectName("open_folder")
        self.open_file = QtWidgets.QAction(MainWindow)
        self.open_file.setObjectName("open_file")
        self.run_abinitio = QtWidgets.QAction(MainWindow)
        self.run_abinitio.setObjectName("run_abinitio")
        self.save = QtWidgets.QAction(MainWindow)
        self.save.setObjectName("save")
        self.menu_abinitio.addAction(self.dmol3_set)
        self.menu_abinitio.addAction(self.vasp_set)
        self.menu_abinitio.addAction(self.gaussian_set)
        self.menu_abinitio.addAction(self.lammps_set)
        self.menu_abinitio.addAction(self.abacus_set)
        self.menu_abinitio.addSeparator()
        self.menu_abinitio.addAction(self.run_abinitio)
        self.menuCGA.addAction(self.new_cga)
        self.menuCGA.addAction(self.continue_cga)
        self.menu_help.addAction(self.action_CGA)
        self.menu_file.addAction(self.open_folder)
        self.menu_file.addAction(self.open_file)
        self.menu_file.addAction(self.save)
        self.menubar.addAction(self.menu_file.menuAction())
        self.menubar.addAction(self.menu_abinitio.menuAction())
        self.menubar.addAction(self.menuCGA.menuAction())
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menu_help.menuAction())
        
        self.open_folder.triggered.connect(self.open_folder_msg)
        self.open_file.triggered.connect(self.open_file_msg)
        self.dmol3_set.triggered.connect(self.dmol3_set_msg)
        self.vasp_set.triggered.connect(self.vasp_set_msg)
        self.gaussian_set.triggered.connect(self.gaussian_set_msg)
        self.lammps_set.triggered.connect(self.lammps_set_msg)
        self.abacus_set.triggered.connect(self.abacus_set_msg)
        self.run_abinitio.triggered.connect(self.run_abinitio_msg)
        self.new_cga.triggered.connect(self.new_cga_msg)
        self.continue_cga.triggered.connect(self.continue_cga_msg)

        self.retranslateUi(MainWindow)
        self.save_2.clicked.connect(MainWindow.save_pop)
        self.pop_info.cellDoubleClicked.connect(MainWindow.show_pop)
        #self.connect(self.pop_info, QtCore.SIGNAL("itemDoubleClicked(*item)"),self.show_pop)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
        self.extra_init()


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "CGA"))
        self.label.setText(_translate("MainWindow", "文件名前缀："))
        self.save_2.setText(_translate("MainWindow", "排序并保存种群结构"))
        self.menu_abinitio.setTitle(_translate("MainWindow", "第一性原理软件"))
        self.menuCGA.setTitle(_translate("MainWindow", "CGA"))
        self.menu_help.setTitle(_translate("MainWindow", "帮助"))
        self.menu_file.setTitle(_translate("MainWindow", "文件"))
        self.menu.setTitle(_translate("MainWindow", "分析"))
        self.dmol3_set.setText(_translate("MainWindow", "DMol3"))
        self.vasp_set.setText(_translate("MainWindow", "VASP"))
        self.gaussian_set.setText(_translate("MainWindow", "Gaussian"))
        self.lammps_set.setText(_translate("MainWindow", "Lammps"))
        self.abacus_set.setText(_translate("MainWindow", "Abacus"))
        self.new_cga.setText(_translate("MainWindow", "新搜索"))
        self.continue_cga.setText(_translate("MainWindow", "继续搜索"))
        self.action_CGA.setText(_translate("MainWindow", "关于CGA"))
        self.open_folder.setText(_translate("MainWindow", "打开文件夹"))
        self.open_file.setText(_translate("MainWindow", "打开文件"))
        self.run_abinitio.setText(_translate("MainWindow", "运行"))
        self.save.setText(_translate("MainWindow", "保存"))
        
        
    def extra_init(self):
        self.energy_figure = EnergyFigure(self)
        self.cluster_figure = ClusterFigure(self)
        self.pop_table_set()
        
        
    def open_folder_msg(self):
        work_folder = QFileDialog.getExistingDirectory(self,"选择文件夹","·")
        os.chdir(work_folder)
        if os.path.isfile('recover.txt'):
            self.energy_figure.plot()
            self.ga = GA()
            self.ga.read_recover()
            self.pop_display(self.ga.pop)
        
        
    def open_file_msg(self):
        self.file_names = QFileDialog.getOpenFileNames(self,"打开文件","C:\Users\Administrator\Desktop","Txt files(*.*)")[0]
        print self.file_names
        self.cluster = Cluster()
        self.cluster.read(self.file_names[0])
        self.cluster_figure.plot(self.cluster)
        
        
    def dmol3_set_msg(self):
        self.dmol3_set.setChecked(True)
        self.vasp_set.setChecked(False)
        self.gaussian_set.setChecked(False)
        self.lammps_set.setChecked(False)
        self.abacus_set.setChecked(False)
        dmol3_set_dialog = Ui_DMol3_set()
        dmol3_set_dialog.exec_()
        
        
    def vasp_set_msg(self):
        self.dmol3_set.setChecked(False)
        self.vasp_set.setChecked(True)
        self.gaussian_set.setChecked(False)
        self.lammps_set.setChecked(False)
        self.abacus_set.setChecked(False)
        vasp_set_dialog = Ui_VASP_Set()
        #vasp_dialog.show()
        vasp_set_dialog.exec_()
        
        
    def gaussian_set_msg(self):
        self.dmol3_set.setChecked(False)
        self.vasp_set.setChecked(False)
        self.gaussian_set.setChecked(True)
        self.lammps_set.setChecked(False)
        self.abacus_set.setChecked(False)
        
        
    def lammps_set_msg(self):
        self.dmol3_set.setChecked(False)
        self.vasp_set.setChecked(False)
        self.gaussian_set.setChecked(False)
        self.lammps_set.setChecked(True)
        self.abacus_set.setChecked(False)
        
        
    def abacus_set_msg(self):
        self.dmol3_set.setChecked(False)
        self.vasp_set.setChecked(False)
        self.gaussian_set.setChecked(False)
        self.lammps_set.setChecked(False)
        self.abacus_set.setChecked(True)
        
        
    def run_abinitio_msg(self):
        print 'run abinitio.'
        pass
        
        
    def new_cga_msg(self):
        cga_set_dialog = Ui_CGA_Set()
        #cga_set_dialog.setupUi(self)
        cga_set_dialog.exec_()
        if cga_set_dialog.run_cga:
            self.continue_cga_msg()
        
        
    def continue_cga_msg(self):
        self.ga = GA()
        self.ga.init_pop()
        self.ga.iterator()
        
        
    def pop_table_set(self):
        self.pop_info.setColumnCount(4)
        self.pop_info.setRowCount(1)
        self.pop_info.setHorizontalHeaderLabels([u'编号',u'能量',u'转动惯量1',u'转动惯量2'])
        self.pop_info.setColumnWidth(0,40)
        self.pop_info.setColumnWidth(1,100)
        self.pop_info.setColumnWidth(2,70)
        self.pop_info.setColumnWidth(3,80)
        
        
    def pop_display(self, pop):
        self.pop_info.setRowCount(len(pop))
        for i in range(len(pop)):
            self.pop_info.setItem(i,0,QTableWidgetItem('%02d'%i))
            self.pop_info.setItem(i,1,QTableWidgetItem('%.6f'%pop[i].get_energy()))
            self.pop_info.setItem(i,2,QTableWidgetItem('%.2f'%pop[i].get_inertia()[0]))
            self.pop_info.setItem(i,3,QTableWidgetItem('%.2f'%pop[i].get_inertia()[1]))
            
            
    def show_pop(self, row, column):
        idx = int(self.pop_info.item(row,0).text())
        self.cluster_figure.plot(self.ga.pop[idx]) 
    
    
    def save_pop(self):
        pop = deepcopy(self.ga.pop)
        prefix = self.prefix.text()
        for i in range(len(pop)):
            idx = np.argmin([c.get_energy() for c in pop])
            pop[idx].write_car('%s_%02d.car' %(prefix,i))
            print pop[idx].get_energy()
            del pop[idx]




class MyWindow(QtWidgets.QWidget,Ui_MainWindow):
    def __init__(self): 
        super(MyWindow,self).__init__() 
        self.setupUi(self)
        
    def closeEvent(self, event):
        print 'close'



if __name__=="__main__": 
    import sys
    if '-f' in sys.argv:
        idx = sys.argv.index('-f') + 1
        if idx == len(sys.argv) or not os.path.isdir(sys.argv[idx]):
            print 'command paramters error.'
            os._exit(0)
        else:
            os.chdir(sys.argv[idx])
    if '-s' in sys.argv:
        ga = GA()
        ga.init_pop()
        ga.iterator()
    elif '-b' in sys.argv:
        pass
    else:
        app = QtWidgets.QApplication(sys.argv) 
        myshow = MyWindow()
        myshow.show()
        sys.exit(app.exec_())
