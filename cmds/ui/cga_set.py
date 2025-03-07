# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cga_set.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QDialog
from ConfigParser import ConfigParser,RawConfigParser
import os
import sys
sys.path.append('..')
from src.atom import get_element_id


class Ui_CGA_Set(QDialog):
    def __init__(self, parent=None):
        super(Ui_CGA_Set,self).__init__(parent)
        self.setupUi(self)
        
    def setupUi(self, CGA_Set):
        CGA_Set.setObjectName("CGA_Set")
        CGA_Set.resize(579, 309)
        CGA_Set.setModal(True)
        self.groupBox = QtWidgets.QGroupBox(CGA_Set)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 291, 201))
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 20, 166, 41))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_cluster_type = QtWidgets.QLabel(self.horizontalLayoutWidget)
        self.label_cluster_type.setObjectName("label_cluster_type")
        self.horizontalLayout.addWidget(self.label_cluster_type)
        self.cluster_type_set = QtWidgets.QComboBox(self.horizontalLayoutWidget)
        self.cluster_type_set.setObjectName("cluster_type_set")
        self.cluster_type_set.addItem("")
        self.cluster_type_set.addItem("")
        self.cluster_type_set.addItem("")
        self.cluster_type_set.addItem("")
        self.horizontalLayout.addWidget(self.cluster_type_set)
        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(10, 70, 270, 41))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_element1 = QtWidgets.QLabel(self.horizontalLayoutWidget_2)
        self.label_element1.setObjectName("label_element1")
        self.horizontalLayout_6.addWidget(self.label_element1)
        self.element1 = QtWidgets.QLineEdit(self.horizontalLayoutWidget_2)
        self.element1.setObjectName("element1")
        self.horizontalLayout_6.addWidget(self.element1)
        self.label_atom_num1 = QtWidgets.QLabel(self.horizontalLayoutWidget_2)
        self.label_atom_num1.setObjectName("label_atom_num1")
        self.horizontalLayout_6.addWidget(self.label_atom_num1)
        self.atom_num1 = QtWidgets.QLineEdit(self.horizontalLayoutWidget_2)
        self.atom_num1.setObjectName("atom_num1")
        self.horizontalLayout_6.addWidget(self.atom_num1)
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(10, 120, 270, 31))
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_element2 = QtWidgets.QLabel(self.horizontalLayoutWidget_3)
        self.label_element2.setObjectName("label_element2")
        self.horizontalLayout_7.addWidget(self.label_element2)
        self.element2 = QtWidgets.QLineEdit(self.horizontalLayoutWidget_3)
        self.element2.setObjectName("element2")
        self.horizontalLayout_7.addWidget(self.element2)
        self.label_atom_num2 = QtWidgets.QLabel(self.horizontalLayoutWidget_3)
        self.label_atom_num2.setObjectName("label_atom_num2")
        self.horizontalLayout_7.addWidget(self.label_atom_num2)
        self.atom_num2 = QtWidgets.QLineEdit(self.horizontalLayoutWidget_3)
        self.atom_num2.setObjectName("atom_num2")
        self.horizontalLayout_7.addWidget(self.atom_num2)
        self.horizontalLayoutWidget_4 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_4.setGeometry(QtCore.QRect(10, 160, 158, 31))
        self.horizontalLayoutWidget_4.setObjectName("horizontalLayoutWidget_4")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_4)
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_9 = QtWidgets.QLabel(self.horizontalLayoutWidget_4)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_9.addWidget(self.label_9)
        self.symmetry_set = QtWidgets.QComboBox(self.horizontalLayoutWidget_4)
        self.symmetry_set.setObjectName("symmetry_set")
        self.symmetry_set.addItem("")
        self.symmetry_set.addItem("")
        self.symmetry_set.addItem("")
        self.symmetry_set.addItem("")
        self.symmetry_set.addItem("")
        self.horizontalLayout_9.addWidget(self.symmetry_set)
        self.groupBox_2 = QtWidgets.QGroupBox(CGA_Set)
        self.groupBox_2.setGeometry(QtCore.QRect(320, 10, 251, 291))
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayoutWidget_5 = QtWidgets.QWidget(self.groupBox_2)
        self.horizontalLayoutWidget_5.setGeometry(QtCore.QRect(10, 20, 231, 31))
        self.horizontalLayoutWidget_5.setObjectName("horizontalLayoutWidget_5")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_5)
        self.horizontalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label_pop_size = QtWidgets.QLabel(self.horizontalLayoutWidget_5)
        self.label_pop_size.setObjectName("label_pop_size")
        self.horizontalLayout_10.addWidget(self.label_pop_size)
        self.popular_size = QtWidgets.QLineEdit(self.horizontalLayoutWidget_5)
        self.popular_size.setObjectName("popular_size")
        self.horizontalLayout_10.addWidget(self.popular_size)
        self.label_iter_num = QtWidgets.QLabel(self.horizontalLayoutWidget_5)
        self.label_iter_num.setObjectName("label_iter_num")
        self.horizontalLayout_10.addWidget(self.label_iter_num)
        self.iterator_num = QtWidgets.QLineEdit(self.horizontalLayoutWidget_5)
        self.iterator_num.setObjectName("iterator_num")
        self.horizontalLayout_10.addWidget(self.iterator_num)
        self.label_12 = QtWidgets.QLabel(self.groupBox_2)
        self.label_12.setGeometry(QtCore.QRect(20, 60, 54, 12))
        self.label_12.setObjectName("label_12")
        self.label_13 = QtWidgets.QLabel(self.groupBox_2)
        self.label_13.setGeometry(QtCore.QRect(120, 60, 54, 12))
        self.label_13.setObjectName("label_13")
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 80, 160, 206))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_method4 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method4.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method4.setObjectName("label_method4")
        self.gridLayout.addWidget(self.label_method4, 3, 0, 1, 1)
        self.label_method3 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method3.setObjectName("label_method3")
        self.gridLayout.addWidget(self.label_method3, 2, 0, 1, 1)
        self.label_method1 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method1.setObjectName("label_method1")
        self.gridLayout.addWidget(self.label_method1, 0, 0, 1, 1)
        self.perturbation_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.perturbation_prob.setObjectName("perturbation_prob")
        self.gridLayout.addWidget(self.perturbation_prob, 1, 1, 1, 1)
        self.mating_principle_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.mating_principle_prob.setObjectName("mating_principle_prob")
        self.gridLayout.addWidget(self.mating_principle_prob, 3, 1, 1, 1)
        self.mating_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.mating_prob.setObjectName("mating_prob")
        self.gridLayout.addWidget(self.mating_prob, 0, 1, 1, 1)
        self.mutation_one_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.mutation_one_prob.setObjectName("mutation_one_prob")
        self.gridLayout.addWidget(self.mutation_one_prob, 2, 1, 1, 1)
        self.label_method6 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method6.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method6.setObjectName("label_method6")
        self.gridLayout.addWidget(self.label_method6, 5, 0, 1, 1)
        self.label_method5 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method5.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method5.setObjectName("label_method5")
        self.gridLayout.addWidget(self.label_method5, 4, 0, 1, 1)
        self.label_method2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method2.setObjectName("label_method2")
        self.gridLayout.addWidget(self.label_method2, 1, 0, 1, 1)
        self.exchange_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.exchange_prob.setObjectName("exchange_prob")
        self.gridLayout.addWidget(self.exchange_prob, 6, 1, 1, 1)
        self.label_method7 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method7.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method7.setObjectName("label_method7")
        self.gridLayout.addWidget(self.label_method7, 6, 0, 1, 1)
        self.change_naxis_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.change_naxis_prob.setObjectName("change_naxis_prob")
        self.gridLayout.addWidget(self.change_naxis_prob, 5, 1, 1, 1)
        self.horizontalLayout_16 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.mating_radial_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.mating_radial_prob.setObjectName("mating_radial_prob")
        self.horizontalLayout_16.addWidget(self.mating_radial_prob)
        self.gridLayout.addLayout(self.horizontalLayout_16, 4, 1, 1, 1)
        self.label_method8 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_method8.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_method8.setObjectName("label_method8")
        self.gridLayout.addWidget(self.label_method8, 7, 0, 1, 1)
        self.adjust_H_prob = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.adjust_H_prob.setObjectName("adjust_H_prob")
        self.gridLayout.addWidget(self.adjust_H_prob, 7, 1, 1, 1)
        self.horizontalLayoutWidget_7 = QtWidgets.QWidget(CGA_Set)
        self.horizontalLayoutWidget_7.setGeometry(QtCore.QRect(70, 270, 160, 31))
        self.horizontalLayoutWidget_7.setObjectName("horizontalLayoutWidget_7")
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_7)
        self.horizontalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.saveButton = QtWidgets.QPushButton(self.horizontalLayoutWidget_7)
        self.saveButton.setObjectName("saveButton")
        self.horizontalLayout_12.addWidget(self.saveButton)
        self.runButton = QtWidgets.QPushButton(self.horizontalLayoutWidget_7)
        self.runButton.setObjectName("runButton")
        self.horizontalLayout_12.addWidget(self.runButton)

        self.retranslateUi(CGA_Set)
        self.cluster_type_set.currentIndexChanged['QString'].connect(CGA_Set.change_cluster_type)
        self.symmetry_set.currentIndexChanged['QString'].connect(CGA_Set.change_symmetry)
        self.saveButton.clicked.connect(CGA_Set.save)
        self.runButton.clicked.connect(CGA_Set.run)
        QtCore.QMetaObject.connectSlotsByName(CGA_Set)

        self.extra_init(CGA_Set)

    def retranslateUi(self, CGA_Set):
        _translate = QtCore.QCoreApplication.translate
        CGA_Set.setWindowTitle(_translate("CGA_Set", "CGA设置"))
        self.groupBox.setTitle(_translate("CGA_Set", "团簇体系设置"))
        self.label_cluster_type.setText(_translate("CGA_Set", "团簇类型："))
        self.cluster_type_set.setItemText(0, _translate("CGA_Set", "单质团簇"))
        self.cluster_type_set.setItemText(1, _translate("CGA_Set", "二元合金"))
        self.cluster_type_set.setItemText(2, _translate("CGA_Set", "纯分子团簇"))
        self.cluster_type_set.setItemText(3, _translate("CGA_Set", "掺杂分子团簇"))
        self.label_element1.setText(_translate("CGA_Set", "元素1："))
        self.label_atom_num1.setText(_translate("CGA_Set", "原子数1："))
        self.label_element2.setText(_translate("CGA_Set", "元素2："))
        self.label_atom_num2.setText(_translate("CGA_Set", "原子数2："))
        self.label_9.setText(_translate("CGA_Set", "对称性："))
        self.symmetry_set.setItemText(0, _translate("CGA_Set", "C1"))
        self.symmetry_set.setItemText(1, _translate("CGA_Set", "Cs"))
        self.symmetry_set.setItemText(2, _translate("CGA_Set", "C2"))
        self.symmetry_set.setItemText(3, _translate("CGA_Set", "C3"))
        self.symmetry_set.setItemText(4, _translate("CGA_Set", "C5"))
        self.groupBox_2.setTitle(_translate("CGA_Set", "遗传算法设置"))
        self.label_pop_size.setText(_translate("CGA_Set", "种群大小："))
        self.popular_size.setText(_translate("CGA_Set", "16"))
        self.label_iter_num.setText(_translate("CGA_Set", "迭代次数："))
        self.iterator_num.setText(_translate("CGA_Set", "1000"))
        self.label_12.setText(_translate("CGA_Set", "遗传方法"))
        self.label_13.setText(_translate("CGA_Set", "概率"))
        self.label_method4.setText(_translate("CGA_Set", "主方向交叉"))
        self.label_method3.setText(_translate("CGA_Set", "移动一个原子"))
        self.label_method1.setText(_translate("CGA_Set", "交叉 "))
        self.perturbation_prob.setText(_translate("CGA_Set", "0"))
        self.mating_principle_prob.setText(_translate("CGA_Set", "0"))
        self.mating_prob.setText(_translate("CGA_Set", "0"))
        self.mutation_one_prob.setText(_translate("CGA_Set", "0"))
        self.label_method6.setText(_translate("CGA_Set", "变异对称轴"))
        self.label_method5.setText(_translate("CGA_Set", "径向交叉"))
        self.label_method2.setText(_translate("CGA_Set", "微扰变异"))
        self.exchange_prob.setText(_translate("CGA_Set", "0"))
        self.label_method7.setText(_translate("CGA_Set", "交换元素"))
        self.change_naxis_prob.setText(_translate("CGA_Set", "0"))
        self.mating_radial_prob.setText(_translate("CGA_Set", "0"))
        self.label_method8.setText(_translate("CGA_Set", "调整氢键"))
        self.adjust_H_prob.setText(_translate("CGA_Set", "0"))
        self.saveButton.setText(_translate("CGA_Set", "保存"))
        self.runButton.setText(_translate("CGA_Set", "运行"))


    def extra_init(self, CGA_Set):
        '''其它初始化工作'''
        CGA_Set.setWindowFlags(QtCore.Qt.WindowCloseButtonHint)
        self.symmetry = self.symmetry_set.currentText()
        self.change_cluster_type()
        '''self.method_list = [self.mating_prob,self.perturbation_prob,self.mutation_one_prob,
            self.mating_principle_prob,self.mating_radial_prob,self.change_naxis_prob,
            self.exchange_prob, self.adjust_H_prob]
        self.read_config()'''

    def read_config(self):
        '''读取配置文件'''
        config = ConfigParser()
        config.read('config.ini')

        self.symmetry = config.get('cluster', 'symmetry')
        elements = config.get('cluster', 'element').split()
        if len(elements) == 1:
            if get_element_id(elements[0]) != -1:
                self.cluster_type = u'单质团簇'
            else:
                self.cluster_type = u'纯分子团簇'
        elif len(elements) == 2:
            if get_element_id(elements[0]) != -1:
                self.cluster_type = u'二元合金'
            else:
                self.cluster_type = u'掺杂分子团簇'
        else:
            pass

        """#读取遗传算法参数: 种群大小,迭代次数
        self.popular_size = config.getint('genetic algorithm', 'pop_size')
        self.iterator_num = config.getint('genetic algorithm', 'max_gen')
        #读取产生后代方式及概率
        self.method_str = config.get('genetic algorithm', 'breed_method').split(',')
        self.probability = config.get('genetic algorithm', 'breed_probability').split(',')
        assert len(self.method_str) == len(self.probability)
        self.method = []
        for i in range(len(self.method_str)):
            self.method_str[i] = self.method_str[i].strip()
            self.probability[i] = float(self.probability[i])
            if not hasattr(self.class_name, self.method_str[i]):
                print 'no function %s in' %self.method_str[i], self.class_name
            self.method.append(getattr(self.class_name, self.method_str[i]))
        for i in range(len(self.method_str)):
            if self.method_str == 'mating':
                pass"""

    def set_method_prob(self):
        '''根据self.cluster_type和self.symmetry的值来设置遗传方法是否可用和概率的默认值'''
        if self.cluster_type == u'单质团簇':
            self.exchange_prob.setEnabled(False)
            self.adjust_H_prob.setEnabled(False)
            self.mutation_one_prob.setEnabled(True)
            if self.symmetry == 'C1':
                self.change_naxis_prob.setEnabled(False)
                self.mating_radial_prob.setEnabled(False)
                self.mating_principle_prob.setEnabled(True)
                self.mating_prob.setText('0.4')
                self.perturbation_prob.setText('0.1')
                self.mutation_one_prob.setText('0.1')
                self.mating_principle_prob.setText('0.4')
            else:
                self.change_naxis_prob.setEnabled(True)
                self.mating_radial_prob.setEnabled(False)
                self.mating_principle_prob.setEnabled(False)
                self.mating_prob.setText('0.7')
                self.perturbation_prob.setText('0.1')
                self.mutation_one_prob.setText('0.1')
                self.change_naxis_prob.setText('0.1')
        elif self.cluster_type == u'二元合金':
            self.exchange_prob.setEnabled(True)
            self.adjust_H_prob.setEnabled(False)
            self.mutation_one_prob.setEnabled(True)
            if self.symmetry == 'C1':
                self.change_naxis_prob.setEnabled(False)
                self.mating_radial_prob.setEnabled(False)
                self.mating_principle_prob.setEnabled(True)
                self.mating_prob.setText('0.4')
                self.perturbation_prob.setText('0.1')
                self.mutation_one_prob.setText('0.1')
                self.mating_principle_prob.setText('0.3')
                self.exchange_prob.setText('0.1')
            else:
                self.change_naxis_prob.setEnabled(True)
                self.mating_radial_prob.setEnabled(False)
                self.mating_principle_prob.setEnabled(False)
                self.mating_prob.setText('0.6')
                self.perturbation_prob.setText('0.1')
                self.mutation_one_prob.setText('0.1')
                self.exchange_prob.setText('0.1')
                self.change_naxis_prob.setText('0.1')
        elif self.cluster_type == u'纯分子团簇' or self.cluster_type == u'掺杂分子团簇':
            self.exchange_prob.setEnabled(False)
            self.mutation_one_prob.setEnabled(False)
            self.mating_radial_prob.setEnabled(False)
            self.change_naxis_prob.setEnabled(False)
            self.change_naxis_prob.setEnabled(False)
            self.mating_radial_prob.setEnabled(False)
            self.mating_principle_prob.setEnabled(False)
            self.adjust_H_prob.setEnabled(True)
            self.mating_prob.setText('0.7')
            self.perturbation_prob.setText('0.1')
            self.adjust_H_prob.setText('0.2')

    def change_cluster_type(self):
        '''根据选择的团簇类型，改变团簇体系参数的状态'''
        self.cluster_type = self.cluster_type_set.currentText()
        if self.cluster_type == u'单质团簇':
            self.label_element1.setText(u'元素1')
            self.label_atom_num1.setText(u'原子数1')
            self.element2.setEnabled(False)
            self.atom_num2.setEnabled(False)
            self.symmetry_set.setEnabled(True)
        elif self.cluster_type == u'二元合金':
            self.label_element1.setText(u'元素1')
            self.label_atom_num1.setText(u'原子数1')
            self.label_element2.setText(u'元素2')
            self.label_atom_num2.setText(u'原子数2')
            self.element2.setEnabled(True)
            self.atom_num2.setEnabled(True)
            self.symmetry_set.setEnabled(True)
        elif self.cluster_type == u'纯分子团簇':
            self.label_element1.setText(u'分子结构')
            self.label_atom_num1.setText(u'分子数')
            self.element2.setEnabled(False)
            self.atom_num2.setEnabled(False)
            self.symmetry_set.setEnabled(False)
        elif self.cluster_type == u'掺杂分子团簇':
            self.label_element1.setText(u'主体分子结构')
            self.label_atom_num1.setText(u'主体分子数')
            self.label_element2.setText(u'客体分子结构')
            self.label_atom_num2.setText(u'客体分子数')
            self.element2.setEnabled(True)
            self.atom_num2.setEnabled(True)
            self.symmetry_set.setEnabled(False)
        self.set_method_prob()

    def change_symmetry(self):
        self.symmetry = self.symmetry_set.currentText()
        self.set_method_prob()
        
    def get_method_prob(self):
        '''获取选择的遗传方法和概率
        保存在self.method_list和self.prob_list'''
        self.method_list = []
        self.prob_list = []
        if float(self.mating_prob.text())!=0.:
            self.method_list.append('mating')
            self.prob_list.append(self.mating_prob.text())
        if float(self.perturbation_prob.text())!=0.:
            self.method_list.append('perturbation')
            self.prob_list.append(self.perturbation_prob.text())
        if self.mutation_one_prob.isEnabled() and float(self.mutation_one_prob.text())!=0.:
            self.method_list.append('mutation_one_better')
            self.prob_list.append(self.mutation_one_prob.text())
        if self.mating_principle_prob.isEnabled() and float(self.mating_principle_prob.text())!=0.:
            self.method_list.append('mating_principal_direction')
            self.prob_list.append(self.mating_principle_prob.text())
        if self.mating_radial_prob.isEnabled() and float(self.mating_radial_prob.text())!=0.:
            self.method_list.append('mating_radial')
            self.prob_list.append(self.mating_radial_prob.text())
        if self.change_naxis_prob.isEnabled() and float(self.change_naxis_prob.text())!=0.:
            self.method_list.append('change_naxis')
            self.prob_list.append(self.change_naxis_prob.text())
        if self.exchange_prob.isEnabled() and float(self.exchange_prob.text())!=0.:
            self.method_list.append('exchange')
            self.prob_list.append(self.exchange_prob.text())
        if self.adjust_H_prob.isEnabled() and float(self.adjust_H_prob.text())!=0.:
            self.method_list.append('adjust_H')
            self.prob_list.append(self.adjust_H_prob.text())
            
    def check(self):
        '''检查填写的参数是否正确'''
        if get_element_id(self.element1.text())==-1 and not os.path.isfile(self.element1.text()):
            QtWidgets.QMessageBox.critical(self, '', u'元素（分子）1错误')
            return False       
        if self.element2.isEnabled():
            if get_element_id(self.element2.text())==-1 and not os.path.isfile(self.element2.text()):
                QtWidgets.QMessageBox.critical(self, '', u'元素（分子）2错误')
                return False
        try:
            if int(self.atom_num1.text()) <= 0:
                raise
        except:
            QtWidgets.QMessageBox.critical(self, '', u'原子（分子）数1错误')
            return False
        if self.atom_num2.isEnabled():
            try:
                if int(self.atom_num2.text()) <= 0:
                    raise
            except:
                QtWidgets.QMessageBox.critical(self, '', u'原子（分子）数2错误')
                return False        
        try:
            self.get_method_prob()
            assert 0.9 < sum([float(prob) for prob in self.prob_list]) < 1.1
        except:
            QtWidgets.QMessageBox.critical(self, '', u'概率错误')
            return False
        return True

    def save(self):
        '''保存参数到config.ini文件'''
        if not self.check():
            return
        cfg = RawConfigParser()
        cfg.read('config.ini')
        
        #写cluster节
        if not cfg.has_section('cluster'):
            cfg.add_section('cluster')
        elem = self.element1.text()
        if self.cluster_type == u'二元合金' or self.cluster_type == u'掺杂分子团簇':
            elem += ' ' + self.element2.text()
        cfg.set('cluster', 'element', elem)         
        atom_num = self.atom_num1.text()
        if self.cluster_type == u'二元合金' or self.cluster_type == u'掺杂分子团簇':
            atom_num += ' ' + self.atom_num2.text()
        cfg.set('cluster', 'atom_num', atom_num)
        cfg.set('cluster', 'symmetry', self.symmetry_set.currentText())
        
        #写genetic algorithm节
        if not cfg.has_section('genetic algorithm'):
            cfg.add_section('genetic algorithm')
        cfg.set('genetic algorithm', 'pop_size', self.popular_size.text())
        cfg.set('genetic algorithm', 'max_gen', self.iterator_num.text())
        cfg.set('genetic algorithm', 'breed_method', ','.join(self.method_list))
        cfg.set('genetic algorithm', 'breed_probability', ','.join(self.prob_list))
        
        cfg.write(open('config.ini','r+'))
        self.close()
        self.run_cga = False #告诉主窗口是否运行CGA

    def run(self):
        self.save()
        self.close()
        self.run_cga = True #告诉主窗口是否运行CGA

