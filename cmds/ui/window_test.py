from PyQt5 import QtWidgets 
from cga_set import Ui_CGA_Set
from vasp_set import Ui_VASP_Set 
from cga import Ui_MainWindow
from cga_display import Ui_CGA_Display

class MyWindow(QtWidgets.QWidget,Ui_CGA_Display):
    def __init__(self): 
        super(MyWindow,self).__init__() 
        self.setupUi(self)
        
    def closeEvent(self, event):
        print 'close'

if __name__=="__main__": 
    import sys 
    app = QtWidgets.QApplication(sys.argv) 
    myshow = MyWindow() 
    myshow.show() 
    sys.exit(app.exec_())