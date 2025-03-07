# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dmol3_set.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QFileDialog
import shutil
import os
import platform
from ConfigParser import ConfigParser
import sys
sys.path.append('..')
from src.abinitio import Abinitio

class Ui_DMol3_set(QDialog):
    def __init__(self, parent=None):
        super(Ui_DMol3_set,self).__init__(parent)
        self.setupUi(self)
        self.setFixedSize(self.width(), self.height())
        
    def setupUi(self, dmol3_set):
        dmol3_set.setObjectName("dmol3_set")
        dmol3_set.resize(406, 92)
        dmol3_set.setMinimumSize(QtCore.QSize(0, 0))
        dmol3_set.setSizeIncrement(QtCore.QSize(0, 0))
        dmol3_set.setModal(True)
        self.horizontalLayoutWidget = QtWidgets.QWidget(dmol3_set)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 50, 378, 31))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setSpacing(50)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(self.horizontalLayoutWidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.core_num = QtWidgets.QSpinBox(self.horizontalLayoutWidget)
        self.core_num.setMinimum(1)
        self.core_num.setProperty("value", 1)
        self.core_num.setObjectName("core_num")
        self.horizontalLayout.addWidget(self.core_num)
        self.pushButton = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout.addWidget(self.pushButton)
        self.pushButton_2 = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.horizontalLayout.addWidget(self.pushButton_2)
        self.label_2 = QtWidgets.QLabel(dmol3_set)
        self.label_2.setGeometry(QtCore.QRect(10, 20, 131, 16))
        self.label_2.setObjectName("label_2")
        self.ms_path = QtWidgets.QLineEdit(dmol3_set)
        self.ms_path.setGeometry(QtCore.QRect(150, 20, 241, 20))
        self.ms_path.setObjectName("ms_path")

        self.retranslateUi(dmol3_set)
        self.pushButton.clicked.connect(dmol3_set.import_input)
        self.pushButton_2.clicked.connect(dmol3_set.ok)
        QtCore.QMetaObject.connectSlotsByName(dmol3_set)

    def retranslateUi(self, dmol3_set):
        _translate = QtCore.QCoreApplication.translate
        dmol3_set.setWindowTitle(_translate("dmol3_set", "DMol3导入"))
        self.label.setText(_translate("dmol3_set", "核数："))
        self.pushButton.setText(_translate("dmol3_set", "导入input"))
        self.pushButton_2.setText(_translate("dmol3_set", "确定"))
        self.label_2.setText(_translate("dmol3_set", "Materials Studio路径："))        
        
    def import_input(self):
        file_name = QFileDialog.getOpenFileName(self,"打开文件","C:\Users\Administrator\Desktop","input file(*.input)")[0]
        shutil.copyfile(file_name, 'optimize.input')
    
    def ok(self):
        cfg = ConfigParser()
        if os.path.isfile('config.ini'):
            cfg.read('config.ini')
        if not os.path.isfile('config.ini') or 'abinitio' not in cfg.sections():
            cfg.add_section('abinitio')
        cfg.set('abinitio', 'core_num', self.core_num.text())
        if platform.system() == 'Windows':
            cfg.set('abinitio', 'accelrys_path_windows', self.ms_path.text())
        else:
            cfg.set('abinitio', 'accelrys_path_linux', self.ms_path.text())
        cfg.set('abinitio', 'method', 'dmol3')
        f = open('config.ini', 'w')
        cfg.write(f)
        Abinitio.config()
        self.close()

