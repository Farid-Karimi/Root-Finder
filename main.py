from math import *
from PyQt5.QtWidgets import QFileDialog, QLineEdit, QInputDialog, QWidget, QApplication
from PyQt5 import QtCore, QtGui, QtWidgets
import time
import sys
import sympy as sym
from sympy import *
from scipy.misc import derivative

class Ui_MainWindow(object):

    # -------------------------------Root finding implementations------------------------------------

    # Bisection method implementation
    def bisection(self, f1, a, b, N, e):
        print("bisection function")
        if self.f(f1, a)*self.f(f1, b) >= 0:
            self.restb.setText("fails")
            print("Bisection method fails.")
            # methodfail()
            # self.MethodError.setText("Method fails with the given input")
            # self.MethodError.show()
            return None
        m_old = 0.0
        self.errortb.insertPlainText("--\n")
        for n in range(1, N+1):
            if n == 1:
                self.reserror.insertPlainText(str(n)+"\n")
            m_n = (a + b)/2
            f_m_n = self.f(f1, m_n)
            if n > 1:
                self.reserror.insertPlainText(str(n)+"\n")
                error = abs(m_n - m_old) / abs(m_n)
                self.errortb.insertPlainText("%.8f" % error + "\n")
                if error <= e:
                    break
            if self.f(f1, a)*f_m_n < 0:
                b = m_n
            elif self.f(f1, b)*f_m_n < 0:
                a = m_n
            elif f_m_n == 0:
                print("Found exact solution.")
                return m_n
            else:
                print("Bisection method fails.")
                return None
            m_old = m_n
        return (a + b)/2

    # Implementing False Position Method
    def falsePosition(self, f1, a, b, N, e):
        print("false position function")
        if self.f(f1, a) * self.f(f1, b) >= 0:
            self.restb.setText("fails")
            print("False position method fails.")
            return None
        step = 1
        m_old = 0.0
        m_n = 0.0
        print('\n\n False Position Method')
        self.errortb.insertPlainText("--\n")
        for n in range(1, N+1):
            m_n = (a * self.f(f1, b) - b * self.f(f1, a)) / (self.f(f1, b) - self.f(f1, a))
            if n == 1:
                self.reserror.insertPlainText(str(n)+"\n")
            if n > 1:
                self.reserror.insertPlainText(str(n)+"\n")
                error = abs(m_n - m_old) / abs(m_n)
                self.errortb.insertPlainText("%.8f" % error + "\n")
                if error <= e:
                    break
            print('Iteration-%d, x2 = %0.6f and f(x2) = %0.6f' % (step, m_n, self. f(f1, m_n)))
            if self.f(f1, a) * self.f(f1, m_n) < 0:
                b = m_n
            elif self.f(f1, b) * self.f(f1, m_n) < 0:
                a = m_n
            elif self.f(f1, m_n) == 0:
                return m_n
            step = step + 1
            m_old = m_n
        print('\nRequired root is: %0.8f' % m_n)
        return m_n

    # Implementing Newton Raphson Method
    def newton(self, f1, xi, N, e):
        print("newton function")
        Df = self.diff_func(f1, xi)
        print(Df)
        m_old = 0.0
        self.errortb.insertPlainText("--\n")
        for n in range(1, N+1):
            if n == 1:
                self.reserror.insertPlainText(str(n)+"\n")
            fxi = self.f(f1, xi)
            Dfxi = self.diff_func(f1, xi)
            if Dfxi == 0:
                self.restb.setText("zero derivative")
                print('Zero derivative. No solution found.')
                return
            xi = float(xi - (fxi / Dfxi))
            if n > 1:
                self.reserror.insertPlainText(str(n)+"\n")
                error = abs((xi - m_old)) / abs(xi)
                self.errortb.insertPlainText("%.8f" % error + "\n")
                if error <= e:
                    print('Found solution after', n, 'iterations.')
                    return xi
            m_old = xi
        return xi

    # Implementing Secant Method
    def secant(self, f1, x0, x1, N, e):
        print("secant function")
        self.errortb.insertPlainText("--\n")
        for n in range(1, N+1):
            if n == 1:
                self.reserror.insertPlainText(str(n)+"\n")
            if self.f(f1, x0) == self.f(f1, x1):
                self.restb.setText("fail")
                print('Divide by zero error!')
                self.restb.setText("Divide by 0")
                break
            x2 = float(x0 - (x1 - x0) * self.f(f1, x0) / (self.f(f1, x1) - self.f(f1, x0)))
            if n > 1:
                self.reserror.insertPlainText(str(n)+"\n")
                error = abs((x2 - x1)) / abs(x2)
                self.errortb.insertPlainText("%.8f" % error + "\n")
                if error <= e:
                    print('Found solution after', n, 'iterations.')
                    return x2
            x0 = x1
            x1 = x2
        return None

    # Implementing Fixed Point Iteration Method
    def fixedPointIteration(self, g, x0, N, e):
        print("fixed point function")
        f1 = self.diff_func(g, x0)
        if abs(f1) >= 1:
            self.restb.setText("diverge")
            print('diverage')
            return
        x1=0.0
        self.errortb.insertPlainText("--\n")
        for i in range(1, N+1):
            x1 = self.f(g, x0)
            if i == 1:
                self.reserror.insertPlainText(str(i)+"\n")

            if i > 1:
                self.reserror.insertPlainText(str(i)+"\n")
                error = abs(float(x1-x0)/x1)
                self.errortb.insertPlainText("%.8f" % error + "\n")
                if error < e:
                    print(i)
                    return x1
            x0 = x1
        print(i)
        return x1

    def Calc(self, m):
        self.EnterRequirementL.hide()
        self.X0L.hide()
        self.X0TB.hide()
        self.X1L.hide()
        self.X1TB.hide()
        self.EnterIntervalL.hide()
        self.NLabel.hide()
        self.NtextBox.hide()
        self.ErrorLabel.hide()
        self.ErrorTextBox.hide()
        self.CalculateButton.hide()
        self.clearButton.hide()
        self.minusButton.hide()
        self.mulButton.hide()
        self.divButton.hide()
        self.sqrtButton.hide()
        self.leftBracket.hide()
        self.rightBracket.hide()
        self.cosButton.hide()
        self.sinButton.hide()
        self.powerButton.hide()
        self.plusButton.hide()
        self.button7.hide()
        self.button4.hide()
        self.button1.hide()
        self.expButton.hide()
        self.button8.hide()
        self.button5.hide()
        self.button2.hide()
        self.button0.hide()
        self.button9.hide()
        self.button6.hide()
        self.button3.hide()
        self.buttonx.hide()
        self.buttonDot.hide()
        self.FGxTextBox.hide()
        self.FxL.hide()
        self.buttonFile.hide()
        try:
            self.reserror.show()
            self.reslb.show()
            self.errorlb.show()
            self.errortb.show()
            start = time.time()
            n = 0
            er = 0.0
            if self.NtextBox.text() == "":
                n = 50
            else:
                n = int(self.NtextBox.text())
            if self.ErrorTextBox.text() == "":
                er = 0.00001
            else:
                er = float(self.ErrorTextBox.text())
            if self.comboBoxMethods.currentIndex() == 0:
                res = self.bisection(m, float(self.X0TB.text()), float(self.X1TB.text()), n, er)
                print(res)
                if res == None:
                    self.resl.show()
                    self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! Wrong Interval</span></p></body></html>")
                    self.errortb.hide()
                    self.reserror.hide()
                    return
            elif self.comboBoxMethods.currentIndex() == 1:
                res = self.falsePosition(m, float(self.X0TB.text()), float(self.X1TB.text()), n, er)
                print(res)
                if res == None:
                    self.resl.show()
                    self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! Wrong Interval</span></p></body></html>")
                    self.errortb.hide()
                    self.reserror.hide()
                    return
            elif self.comboBoxMethods.currentIndex() == 2:
                res = self.fixedPointIteration(m, float(self.X0TB.text()), n, er)
                print(res)
                if res == None:
                    self.resl.show()
                    self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! Function Diverge</span></p></body></html>")
                    self.errortb.hide()
                    self.reserror.hide()
                    return
            elif self.comboBoxMethods.currentIndex() == 3:
                res = self.newton(m, float(self.X0TB.text()), n, er)
                print(res)
                if res == None:
                    self.resl.show()
                    self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! Zero derivative.</span></p></body></html>")
                    self.errortb.hide()
                    self.reserror.hide()
                    return
            else:
                res = self.secant(m, float(self.X0TB.text()), float(self.X1TB.text()), n, er)
                print(res)
                if res == None:
                    self.resl.show()
                    self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! Divide by zero occurred</span></p></body></html>")
                    self.errortb.hide()
                    self.reserror.hide()
                    return
            end = time.time()
            print(end-start)
            self.timerlb.setText("<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">Time: </span></p></body></html>")
            self.timertb.setText(" %.8f" %float(end-start) + " s")

            self.timertb.show()
            self.timerlb.show()
            self.restb.show()
            self.resl.show()
            self.reslb.setText("<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">N</span></p></body></html>")
            self.errorlb.setText("<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">E</span></p></body></html>")
            self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">Result: </span></p></body></html>")
            self.restb.setText(str(res))
        except:
            self.resl.show()
            self.resl.setText("<html><head/><body><p><span style=\" font-size:11pt; color:red; font-weight:600;\">ERROR! in inputs Try Again</span></p></body></html>")
    def clear(self):
        self.FGxTextBox.clear()
    def tryFunc(self, m):
        self.FGxTextBox.insertPlainText(m)

    # math functions
    def cos(self):
        self.FGxTextBox.insertPlainText("cos(")
    def sin(self):
        self.FGxTextBox.insertPlainText("sin(")
    def exp(self):
        self.FGxTextBox.insertPlainText("exp(")
    def sqrt(self):
        self.FGxTextBox.insertPlainText("sqrt(")
    def power(self):
        self.FGxTextBox.insertPlainText("**")
    def f(self, func, x):
        return eval(func)
    def diff_func(self, x, x0):
        h = 1e-5
        return (self.f(x, x0+h)-self.f(x, x0-h))/(2*h)

    #----------------------------File Handling-----------------------------

    # function to get file
    def getfile(self):
        response = QFileDialog.getOpenFileName(filter="Text files (*.txt)")
        # filename = QFileDialog.getOpenFileName()
        # path = filename[0]
        path = response[0]
        print(path)
        f = open(path, 'r')
        with f:
            data = f.read()
            self.FGxTextBox.setText(data)

    # -------------------------------GUI------------------------------------

    # needs to be refactored and checked
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(871, 646)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.comboBoxMethods = QtWidgets.QComboBox(self.centralwidget)
        self.comboBoxMethods.setGeometry(QtCore.QRect(220, 130, 171, 31))
        self.comboBoxMethods.setObjectName("comboBoxMethods")
        self.comboBoxMethods.addItem("")
        self.comboBoxMethods.addItem("")
        self.comboBoxMethods.addItem("")
        self.comboBoxMethods.addItem("")
        self.comboBoxMethods.addItem("")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(250, 30, 321, 81))
        self.label.setObjectName("label")
        self.ChooseMethodLabel = QtWidgets.QLabel(self.centralwidget)
        self.ChooseMethodLabel.setGeometry(QtCore.QRect(20, 130, 191, 31))
        self.ChooseMethodLabel.setObjectName("ChooseMethodLabel")
        self.confirmMethodButton = QtWidgets.QPushButton(self.centralwidget)
        self.confirmMethodButton.setGeometry(QtCore.QRect(450, 130, 161, 31))
        self.confirmMethodButton.setObjectName("confirmMethodButton")
        self.EnterRequirementL = QtWidgets.QLabel(self.centralwidget)
        self.EnterRequirementL.setGeometry(QtCore.QRect(20, 200, 191, 31))
        self.EnterRequirementL.setObjectName("EnterRequirementL")
        self.FxL = QtWidgets.QLabel(self.centralwidget)
        self.FxL.setGeometry(QtCore.QRect(20, 235, 61, 41))
        self.FxL.setObjectName("FxL")
        self.FGxTextBox = QtWidgets.QTextEdit(self.centralwidget)
        self.FGxTextBox.setGeometry(QtCore.QRect(70, 240, 271, 31))
        self.FGxTextBox.setObjectName("FGxTextBox")
        self.button9 = QtWidgets.QPushButton(self.centralwidget)
        self.button9.setGeometry(QtCore.QRect(20, 300, 51, 28))
        self.button9.setObjectName("buttonx")
        self.buttonx = QtWidgets.QPushButton(self.centralwidget)
        self.buttonx.setGeometry(QtCore.QRect(140, 420, 51, 28))
        self.buttonx.setObjectName("buttonx")
        self.buttonx.setText("x")

        # add file button
        self.buttonFile = QtWidgets.QPushButton(self.centralwidget)
        self.buttonFile.setGeometry(QtCore.QRect(150, 470, 130, 28))
        self.buttonFile.setObjectName("buttonFile")
        self.buttonFile.setText("Open File")

        self.button8 = QtWidgets.QPushButton(self.centralwidget)
        self.button8.setGeometry(QtCore.QRect(80, 300, 51, 28))
        self.button8.setObjectName("button8")
        self.button7 = QtWidgets.QPushButton(self.centralwidget)
        self.button7.setGeometry(QtCore.QRect(140, 300, 51, 28))
        self.button7.setObjectName("button7")
        self.button6 = QtWidgets.QPushButton(self.centralwidget)
        self.button6.setGeometry(QtCore.QRect(20, 330, 51, 28))
        self.button6.setObjectName("button6")
        self.button5 = QtWidgets.QPushButton(self.centralwidget)
        self.button5.setGeometry(QtCore.QRect(80, 330, 51, 28))
        self.button5.setObjectName("button5")
        self.button4 = QtWidgets.QPushButton(self.centralwidget)
        self.button4.setGeometry(QtCore.QRect(140, 330, 51, 28))
        self.button4.setObjectName("button4")
        self.button3 = QtWidgets.QPushButton(self.centralwidget)
        self.button3.setGeometry(QtCore.QRect(20, 360, 51, 28))
        self.button3.setObjectName("button3")
        self.button2 = QtWidgets.QPushButton(self.centralwidget)
        self.button2.setGeometry(QtCore.QRect(80, 360, 51, 28))
        self.button2.setObjectName("button2")
        self.button0 = QtWidgets.QPushButton(self.centralwidget)
        self.button0.setGeometry(QtCore.QRect(80, 390, 51, 28))
        self.button0.setObjectName("button0")
        self.clearButton = QtWidgets.QPushButton(self.centralwidget)
        self.clearButton.setGeometry(QtCore.QRect(290, 300, 93, 28))
        self.clearButton.setObjectName("clearButton")
        self.minusButton = QtWidgets.QPushButton(self.centralwidget)
        self.minusButton.setGeometry(QtCore.QRect(290, 330, 93, 28))
        self.minusButton.setObjectName("minusButton")
        self.mulButton = QtWidgets.QPushButton(self.centralwidget)
        self.mulButton.setGeometry(QtCore.QRect(290, 360, 93, 28))
        self.mulButton.setObjectName("mulButton")
        self.divButton = QtWidgets.QPushButton(self.centralwidget)
        self.divButton.setGeometry(QtCore.QRect(290, 390, 93, 28))
        self.divButton.setObjectName("divButton")
        self.button1 = QtWidgets.QPushButton(self.centralwidget)
        self.button1.setGeometry(QtCore.QRect(140, 360, 51, 28))
        self.button1.setObjectName("button1")
        self.buttonDot = QtWidgets.QPushButton(self.centralwidget)
        self.buttonDot.setGeometry(QtCore.QRect(20, 390, 51, 28))
        self.buttonDot.setObjectName("buttonDot")
        self.plusButton = QtWidgets.QPushButton(self.centralwidget)
        self.plusButton.setGeometry(QtCore.QRect(200, 300, 81, 28))
        self.plusButton.setObjectName("plusButton")
        self.powerButton = QtWidgets.QPushButton(self.centralwidget)
        self.powerButton.setGeometry(QtCore.QRect(200, 330, 81, 28))
        self.powerButton.setObjectName("powerButton")
        self.sinButton = QtWidgets.QPushButton(self.centralwidget)
        self.sinButton.setGeometry(QtCore.QRect(200, 360, 81, 28))
        self.sinButton.setObjectName("sinButton")
        self.cosButton = QtWidgets.QPushButton(self.centralwidget)
        self.cosButton.setGeometry(QtCore.QRect(200, 390, 81, 28))
        self.cosButton.setObjectName("cosButton")
        self.expButton = QtWidgets.QPushButton(self.centralwidget)
        self.expButton.setGeometry(QtCore.QRect(140, 390, 51, 28))
        self.expButton.setObjectName("expButton")
        self.rightBracket = QtWidgets.QPushButton(self.centralwidget)
        self.rightBracket.setGeometry(QtCore.QRect(20, 420, 51, 28))
        self.rightBracket.setObjectName("rightBracket")
        self.leftBracket = QtWidgets.QPushButton(self.centralwidget)
        self.leftBracket.setGeometry(QtCore.QRect(80, 420, 51, 28))
        self.leftBracket.setObjectName("leftBracket")
        self.sqrtButton = QtWidgets.QPushButton(self.centralwidget)
        self.sqrtButton.setGeometry(QtCore.QRect(200, 420, 81, 28))
        self.sqrtButton.setObjectName("sqrtButton")
        self.X0L = QtWidgets.QLabel(self.centralwidget)
        self.X0L.setGeometry(QtCore.QRect(420, 240, 31, 31))
        self.X0L.setObjectName("X0L")
        self.X1L = QtWidgets.QLabel(self.centralwidget)
        self.X1L.setGeometry(QtCore.QRect(580, 240, 31, 31))
        self.X1L.setObjectName("X1L")
        self.X0TB = QtWidgets.QLineEdit(self.centralwidget)
        self.X0TB.setGeometry(QtCore.QRect(460, 240, 91, 31))
        self.X0TB.setObjectName("X0TB")
        self.X1TB = QtWidgets.QLineEdit(self.centralwidget)
        self.X1TB.setGeometry(QtCore.QRect(620, 240, 91, 31))
        self.X1TB.setObjectName("X1TB")
        self.EnterIntervalL = QtWidgets.QLabel(self.centralwidget)
        self.EnterIntervalL.setGeometry(QtCore.QRect(500, 200, 151, 31))
        self.EnterIntervalL.setObjectName("EnterIntervalL")
        self.NLabel = QtWidgets.QLabel(self.centralwidget)
        self.NLabel.setGeometry(QtCore.QRect(410, 300, 141, 41))
        self.NLabel.setObjectName("NLabel")
        self.NtextBox = QtWidgets.QLineEdit(self.centralwidget)
        self.NtextBox.setGeometry(QtCore.QRect(560, 300, 110, 31))
        self.NtextBox.setObjectName("NtextBox")
        self.ErrorLabel = QtWidgets.QLabel(self.centralwidget)
        self.ErrorLabel.setGeometry(QtCore.QRect(410, 360, 131, 41))
        self.ErrorLabel.setObjectName("ErrorLabel")
        self.ErrorTextBox = QtWidgets.QLineEdit(self.centralwidget)
        self.ErrorTextBox.setGeometry(QtCore.QRect(560, 360, 110, 31))
        self.ErrorTextBox.setObjectName("ErrorTextBox")
        self.CalculateButton = QtWidgets.QPushButton(self.centralwidget)
        self.CalculateButton.setGeometry(QtCore.QRect(560, 420, 111, 41))
        self.CalculateButton.setObjectName("CalculateButton")
        self.resl = QtWidgets.QLabel(self.centralwidget)
        self.resl.setGeometry(QtCore.QRect(600, 350, 400, 200))
        self.resl.setObjectName("resl")
        self.restb = QtWidgets.QLabel(self.centralwidget)
        self.restb.setGeometry(QtCore.QRect(670, 350, 300, 200))
        self.restb.setObjectName("restb")
        self.restb.setStyleSheet("background-color: transparent;")
        self.reslb = QtWidgets.QLabel(self.centralwidget)
        self.reslb.setGeometry(QtCore.QRect(60, 200, 80, 150))
        self.reslb.setObjectName("reslb")
        self.reserror = QtWidgets.QTextEdit(self.centralwidget)
        self.reserror.setGeometry(QtCore.QRect(100, 250, 150, 200))
        self.reserror.setObjectName("reserror")
        self.errorlb = QtWidgets.QLabel(self.centralwidget)
        self.errorlb.setGeometry(QtCore.QRect(300, 200, 80, 150))
        self.errorlb.setObjectName("errorlb")
        self.errortb = QtWidgets.QTextEdit(self.centralwidget)
        self.errortb.setGeometry(QtCore.QRect(350, 250, 150, 200))
        self.errortb.setObjectName("errortb")
        self.timerlb = QtWidgets.QLabel(self.centralwidget)
        self.timerlb.setGeometry(QtCore.QRect(600, 400, 200, 150))
        self.timerlb.setObjectName("timerlb")
        self.timerlb.setStyleSheet("background-color: transparent;")
        self.timertb = QtWidgets.QLabel(self.centralwidget)
        self.timertb.setGeometry(QtCore.QRect(650, 400, 200, 150))
        self.timertb.setObjectName("timertb")
        self.timertb.setStyleSheet("background-color: transparent;")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1004, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        #hide buttons
        self.EnterRequirementL.hide()
        self.FxL.hide()
        self.FGxTextBox.hide()
        self.X0L.hide()
        self.X0TB.hide()
        self.X1L.hide()
        self.X1TB.hide()
        self.EnterIntervalL.hide()
        self.NLabel.hide()
        self.NtextBox.hide()
        self.ErrorLabel.hide()
        self.ErrorTextBox.hide()
        self.CalculateButton.hide()
        self.clearButton.hide()
        self.minusButton.hide()
        self.mulButton.hide()
        self.divButton.hide()
        self.sqrtButton.hide()
        self.leftBracket.hide()
        self.rightBracket.hide()
        self.cosButton.hide()
        self.sinButton.hide()
        self.powerButton.hide()
        self.plusButton.hide()
        self.button7.hide()
        self.button4.hide()
        self.button1.hide()
        self.expButton.hide()
        self.button8.hide()
        self.button5.hide()
        self.button2.hide()
        self.button0.hide()
        self.button9.hide()
        self.button6.hide()
        self.button3.hide()
        self.buttonDot.hide()
        self.buttonx.hide()
        self.reserror.hide()
        self.errortb.hide()
        self.errorlb.hide()
        self.buttonFile.hide()
        #Call function showReq to show hidden buttons
        self.confirmMethodButton.clicked.connect(self.showReq)
        self.buttonx.clicked.connect(lambda: self.tryFunc(self.buttonx.text()))
        self.button0.clicked.connect(lambda: self.tryFunc(self.button0.text()))
        self.button1.clicked.connect(lambda: self.tryFunc(self.button1.text()))
        self.button2.clicked.connect(lambda: self.tryFunc(self.button2.text()))
        self.button3.clicked.connect(lambda: self.tryFunc(self.button3.text()))
        self.button4.clicked.connect(lambda: self.tryFunc(self.button4.text()))
        self.button5.clicked.connect(lambda: self.tryFunc(self.button5.text()))
        self.button6.clicked.connect(lambda: self.tryFunc(self.button6.text()))
        self.button7.clicked.connect(lambda: self.tryFunc(self.button7.text()))
        self.button8.clicked.connect(lambda: self.tryFunc(self.button8.text()))
        self.button9.clicked.connect(lambda: self.tryFunc(self.button9.text()))
        self.buttonDot.clicked.connect(lambda: self.tryFunc(self.buttonDot.text()))
        self.cosButton.clicked.connect(lambda: self.cos())
        self.sinButton.clicked.connect(lambda: self.sin())
        self.powerButton.clicked.connect(lambda: self.power())
        self.plusButton.clicked.connect(lambda: self.tryFunc(self.plusButton.text()))
        self.expButton.clicked.connect(lambda: self.exp())
        self.clearButton.clicked.connect(self.clear)
        self.minusButton.clicked.connect(lambda: self.tryFunc(self.minusButton.text()))
        self.mulButton.clicked.connect(lambda: self.tryFunc(self.mulButton.text()))
        self.divButton.clicked.connect(lambda: self.tryFunc(self.divButton.text()))
        self.sqrtButton.clicked.connect(lambda: self.sqrt())
        self.leftBracket.clicked.connect(lambda: self.tryFunc(self.leftBracket.text()))
        self.rightBracket.clicked.connect(lambda: self.tryFunc(self.rightBracket.text()))
        self.CalculateButton.clicked.connect(lambda: self.Calc(self.FGxTextBox.toPlainText()))
        # added line file function
        self.buttonFile.clicked.connect(self.getfile)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        self.set_dark_mode(MainWindow)

    def showReq(self):
        elements_to_clear = [
            self.FGxTextBox,
            self.X0TB,
            self.X1TB,
            self.NtextBox,
            self.ErrorTextBox
        ]
        for element in elements_to_clear:
            element.clear()

        chosen_method = self.comboBoxMethods.currentIndex()
        if chosen_method != 2 or chosen_method != 3:
            self.X1L.show()
            self.X1TB.show()

        if chosen_method == 2:
            self.FxL.setText("G(x)")
        else:
            self.FxL.setText("F(x)")

        if chosen_method == 3 or chosen_method == 2:
            self.X1L.hide()
            self.X1TB.hide()

        elements_to_hide = [
            self.timerlb,
            self.timertb,
            self.resl,
            self.restb,
            self.reslb,
            self.errorlb,
            self.errortb,
            self.reserror
        ]
        for element in elements_to_hide:
            element.hide()
            element.clear()

        elements_to_show = [
            self.EnterRequirementL,
            self.FxL,
            self.FGxTextBox,
            self.X0L,
            self.X0TB,
            self.EnterIntervalL,
            self.NLabel,
            self.NtextBox,
            self.ErrorLabel,
            self.ErrorTextBox,
            self.CalculateButton,
            self.clearButton,
            self.minusButton,
            self.mulButton,
            self.divButton,
            self.sqrtButton,
            self.leftBracket,
            self.rightBracket,
            self.cosButton,
            self.sinButton,
            self.powerButton,
            self.plusButton,
            self.button7,
            self.button4,
            self.button1,
            self.expButton,
            self.button8,
            self.button5,
            self.button2,
            self.button0,
            self.button9,
            self.button6,
            self.button3,
            self.buttonDot,
            self.buttonx,
            self.buttonFile
        ]
        for element in elements_to_show:
            element.show()

    def retranslateUi (self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.comboBoxMethods.setItemText(0, _translate("MainWindow", "Bisection"))
        self.comboBoxMethods.setItemText(1, _translate("MainWindow", "False Position"))
        self.comboBoxMethods.setItemText(2, _translate("MainWindow", "Fixed Point Iteration"))
        self.comboBoxMethods.setItemText(3, _translate("MainWindow", "Newton Raphson"))
        self.comboBoxMethods.setItemText(4, _translate("MainWindow", "Secant"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:22pt; font-weight:600; font-style:italic;\">Root Calculator</span></p></body></html>"))
        self.ChooseMethodLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">Choose Method</span></p></body></html>"))
        self.confirmMethodButton.setText(_translate("MainWindow", "Confirm Method"))
        self.EnterRequirementL.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">Enter Requirements</span></p></body></html>"))
        self.FxL.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">F(x):</span></p></body></html>"))
        self.button9.setText(_translate("MainWindow", "9"))
        self.button8.setText(_translate("MainWindow", "8"))
        self.button7.setText(_translate("MainWindow", "7"))
        self.button6.setText(_translate("MainWindow", "6"))
        self.button5.setText(_translate("MainWindow", "5"))
        self.button4.setText(_translate("MainWindow", "4"))
        self.button3.setText(_translate("MainWindow", "3"))
        self.button2.setText(_translate("MainWindow", "2"))
        self.button0.setText(_translate("MainWindow", "0"))
        self.clearButton.setText(_translate("MainWindow", "C"))
        self.minusButton.setText(_translate("MainWindow", "-"))
        self.mulButton.setText(_translate("MainWindow", "*"))
        self.divButton.setText(_translate("MainWindow", "/"))
        self.button1.setText(_translate("MainWindow", "1"))
        self.buttonDot.setText(_translate("MainWindow", "."))
        self.plusButton.setText(_translate("MainWindow", "+"))
        self.powerButton.setText(_translate("MainWindow", "^"))
        self.sinButton.setText(_translate("MainWindow", "sin(x)"))
        self.cosButton.setText(_translate("MainWindow", "cos(x)"))
        self.expButton.setText(_translate("MainWindow", "exp(x)"))
        self.rightBracket.setText(_translate("MainWindow", "("))
        self.leftBracket.setText(_translate("MainWindow", ")"))
        self.sqrtButton.setText(_translate("MainWindow", "sqrt"))
        self.X0L.setText(_translate("MainWindow", "X0:"))
        self.X1L.setText(_translate("MainWindow", "X1:"))
        self.EnterIntervalL.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">Enter Interval:</span></p></body></html>"))
        self.NLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-weight:600;\">Number Of Iterations: </span></p></body></html>"))
        self.ErrorLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-weight:600;\">Error Tolerance:</span></p></body></html>"))
        self.CalculateButton.setText(_translate("MainWindow", "Calculate"))

    def set_dark_mode(self, MainWindow):
        """
        Apply a dark mode style to the UI elements.
        """
        # Define color variables
        background_color = "#333"  # Dark gray
        text_color = "#fff"  # White
        button_color = "#444"  # Dark gray
        button_hover_color = "#009189"  # Green
        border_color = "#555"  # Dark gray
        selection_color = "#006c64"  # Green
        corner_radius = "8px"  # Rounded corner radius

        # Define the dark mode palette
        dark_palette = QtGui.QPalette()
        dark_palette.setColor(QtGui.QPalette.Window, QtGui.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.WindowText, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.Base, QtGui.QColor(25, 25, 25))
        dark_palette.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.ToolTipBase, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.ToolTipText, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.Text, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.Button, QtGui.QColor(53, 53, 53))
        dark_palette.setColor(QtGui.QPalette.ButtonText, QtCore.Qt.white)
        dark_palette.setColor(QtGui.QPalette.BrightText, QtCore.Qt.red)
        dark_palette.setColor(QtGui.QPalette.Link, QtGui.QColor(42, 130, 218))
        dark_palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor(42, 130, 218))
        dark_palette.setColor(QtGui.QPalette.HighlightedText, QtCore.Qt.black)

        # Set the dark mode palette to the application
        app = QtWidgets.QApplication.instance()
        app.setPalette(dark_palette)

        # Apply custom stylesheets to the MainWindow
        MainWindow.setStyleSheet(f"""
            QWidget {{
                background-color: {background_color};
                color: {text_color};
            }}
            QPushButton {{
                background-color: {button_color};
                border: none;
                padding: 4px 8px;
                border-radius: {corner_radius};
                color: {text_color};
            }}
            QPushButton:hover {{
                background-color: {button_hover_color};
            }}
            QLineEdit, QTextEdit {{
                background-color: {button_color};
                border: 1px solid {border_color};
                padding: 5px;
                border-radius: {corner_radius};
                selection-background-color: {selection_color};
            }}
            QComboBox {{
                background-color: {button_color};
                border: 1px solid {border_color};
                padding: 5px;
                border-radius: {corner_radius};
            }}
            QComboBox::drop-down {{
                subcontrol-origin: padding;
                subcontrol-position: top right;
                width: 20px;
                border-left-width: 1px;
                border-left-color: {border_color};
                border-top-right-radius: {corner_radius};
                border-bottom-right-radius: {corner_radius};
                background-color: {button_color};
            }}
            QComboBox::down-arrow {{
                image: url(down_arrow_icon.png);  /* Replace with your arrow icon path */
                width: 10px;
                height: 10px;
            }}
            QComboBox QAbstractItemView {{
                background-color: {button_color};
                selection-background-color: {button_hover_color};
                color: {text_color};
            }}
        """)

    #-------------------------------------------------------------------

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())