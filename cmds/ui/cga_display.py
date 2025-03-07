# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cga_display.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


#通过继承FigureCanvas类，使得该类既是一个PyQt5的Qwidget，又是一个matplotlib的FigureCanvas
#这是连接pyqt5与matplotlib的关键
class EnergyFigure(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        #创建一个Figure，注意：该Figure为matplotlib下的figure，不是matplotlib.pyplot下面的figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        
        self.axes = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[0.1,0.1,1.,1.]) #否则左边的energy标签会被挡住
        self.plot()
        FigureCanvas.updateGeometry(self)
        
    def plot(self):
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
        self.draw()
        #self.move(350,0)
        

class Ui_CGA_Display(object):
    def setupUi(self, CGA_Display):
        CGA_Display.setObjectName("CGA_Display")
        CGA_Display.resize(760, 613)
        self.pop_info = QtWidgets.QTableView(CGA_Display)
        self.pop_info.setGeometry(QtCore.QRect(10, 290, 311, 301))
        self.pop_info.setObjectName("pop_info")
        self.pop_info.horizontalHeader().setSortIndicatorShown(True)
        self.pop_info.verticalHeader().setVisible(False)
        self.textBrowser = QtWidgets.QTextBrowser(CGA_Display)
        self.textBrowser.setGeometry(QtCore.QRect(10, 10, 311, 231))
        self.textBrowser.setObjectName("textBrowser")
        self.horizontalLayoutWidget = QtWidgets.QWidget(CGA_Display)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 250, 291, 31))
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
        self.save = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.save.setObjectName("save")
        self.horizontalLayout.addWidget(self.save)

        self.retranslateUi(CGA_Display)
        QtCore.QMetaObject.connectSlotsByName(CGA_Display)

    def retranslateUi(self, CGA_Display):
        _translate = QtCore.QCoreApplication.translate
        CGA_Display.setWindowTitle(_translate("CGA_Display", "CGA结果"))
        self.label.setText(_translate("CGA_Display", "文件名前缀："))
        self.save.setText(_translate("CGA_Display", "排序并保存种群结构"))

        #初始化表格
        self.pop_table_set()
        #self.plot_energy()
        
        
    def pop_table_set(self):
        #添加表头：
        self.model = QtGui.QStandardItemModel(self.pop_info)
        #设置表格属性：
        self.model.setRowCount(2)
        self.model.setColumnCount(4)
        #设置表头
        self.model.setHeaderData(0,QtCore.Qt.Horizontal,u'编号')
        self.model.setHeaderData(1,QtCore.Qt.Horizontal,u'能量')
        self.model.setHeaderData(2,QtCore.Qt.Horizontal,u'转动惯量1')
        self.model.setHeaderData(3,QtCore.Qt.Horizontal,u'转动惯量2')
        self.pop_info.setModel(self.model)
        #设置列宽
        self.pop_info.setColumnWidth(0,60)
        self.pop_info.setColumnWidth(1,100)
        self.pop_info.setColumnWidth(2,80)
        self.pop_info.setColumnWidth(3,80)
        
    def pop_display(self, pop):
        for i in len(pop):
            self.model.setItem(i,0,'%02d'%i)
            self.model.setItem(i,1,'%f'%pop[i].get_energy())
            self.model.setItem(i,2,'%.2f'%pop[i].get_inertia()[0])
            self.model.setItem(i,3,'%.2f'%pop[i].get_inertia()[1])
            
    def plot_energy(self):
        print 'here'
        self.energy_figure = EnergyFigure(self)
        self.energy_figure.plot()  # 画图
        self.gridLayoutWidget = QtWidgets.QWidget()
        self.gridLayoutWidget.setGeometry(QtCore.QRect(180, 10, 1100, 500))  # 定义gridLayout控件的大小和位置，4个数字分别为左边坐标，上边坐标，长，宽
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.graphicview = QtWidgets.QGraphicsView(self.gridLayoutWidget)  # 第一步，创建一个QGraphicsView，注意同样以gridLayoutWidget为参
        graphicscene = QtWidgets.QGraphicsScene()  # 第三步，创建一个QGraphicsScene，因为加载的图形（FigureCanvas）不能直接放到graphicview控件中，必须先放到graphicScene，然后再把graphicscene放到graphicview中
        graphicscene.addWidget(self.energy_figure)  # 第四步，把图形放到QGraphicsScene中，注意：图形是作为一个QWidget放到QGraphicsScene中的
        self.graphicview.setScene(graphicscene) # 第五步，把QGraphicsScene放入QGraphicsView
        self.graphicview.show()  # 最后，调用show方法呈现图形！Voila!!
    
    
#a = Ui_CGA_Display()
#a.plot_energy()    


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = EnergyFigure()
    w.show()
    sys.exit(app.exec_())

