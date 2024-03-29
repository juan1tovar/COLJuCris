import matplotlib
# matplotlib.use("TkAgg")  # Linux
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D  # Linux
import numpy as np
from copy import copy

matplotlib.use("TkAgg")
matplotlib.rc('axes', edgecolor='blue', facecolor='314c73', grid=True)
# matplotlib.rc('axes', labelcolor='#D4DCEE')  # compila bien pero no tiene efecto
matplotlib.rc('grid', linestyle='--', alpha=0.65)
matplotlib.rc('xtick', color='#D4DCEE')
matplotlib.rc('ytick', color='#D4DCEE')
matplotlib.rc('figure', facecolor='#303340')
matplotlib.rc('legend', facecolor='#426ca8')
plt.rcParams['figure.constrained_layout.use'] = True


class grafica ():

    def __init__(self, is_3D):
        self.is_3D = is_3D
        self.fig = plt.figure()
        if self.is_3D:
            self.axe = self.fig.add_subplot(projection='3d')
        else:
            self.axe = self.fig.add_subplot()
        # self.fig.legend()

    def show(self):
        # plt.legend()
        plt.show()

    # def update(self):
    #     plt.legend()
    #     self.fig.show()

    def add_axes(self, axe):
        self.fig.add_axes(axe)

    def clf(self):
        self.fig.clf()

    def cla(self):
        self.axe.cla()

    def getAxeLines(self):
        return self.axe.lines

    def add3d(self, nparray, label):
        self.axelines = self.axe.plot(nparray[:, 0], nparray[:, 1], nparray[:, 2], label=label)
        print('grafica 3d añadida')
        # print('nparray3d')
        # print(nparray)

    def add2D(self, nparray, label):
        self.axelines = self.axe.plot(nparray[:, 0], nparray[:, 1], '+-', label=label)
        print('grafica 2D añadida')
        # print('nparray2d')
        # print(nparray)

    def addPMM(self, pmm, label, bycolumns=True):
        if not (bycolumns):
            pmm = np.transpose(pmm)
        nparray = copy(pmm)
        nparray[:, 0] = pmm[:, 1]
        nparray[:, 1] = pmm[:, 2]
        nparray[:, 2] = pmm[:, 0]
        self.add3d(nparray[:, 0:3], label)

    def addPM(self, pmm, label, bycolumns=True):
        if not (bycolumns):
            pmm = np.transpose(pmm)
        nparray = copy(pmm)
        nparray[:, 0] = (pmm[:, 1]**2+pmm[:, 2]**2)**0.5
        nparray[:, 1] = pmm[:, 0]
        self.add2D(nparray[:, 0:2], label)

    def addVer(self, pmm, label, bycolumns=True):
        if self.is_3D:
            self.addPMM(pmm, label, bycolumns)
        else:
            self.addPM(pmm, label, bycolumns)

    def addSecc(self, sec, label):
        nparray = sec.cord_conc(0)
        self.axelines = self.axe.plot(nparray[:, 0], nparray[:, 1], label=label)        
        nparray = sec.refXY
        self.axelines = self.axe.plot(nparray[:, 0], nparray[:, 1], 'o', label=label)
        print('Sección 2D añadida')

    def addSolicitaciones(self, solicitaciones, label, bycolumns=True):      
        if not (bycolumns):
            solicitaciones = np.transpose(solicitaciones)
        if self.is_3D:
            self.axelines = self.axe.plot(solicitaciones[:, 0], solicitaciones[:, 1], solicitaciones[:, 2], 'o', label=label)
            print('solicitaciones 3D añadidas')
        else:
            nparray = copy(solicitaciones)
            nparray[:, 0] = (solicitaciones[:, 1]**2+solicitaciones[:, 2]**2)**0.5
            nparray[:, 1] = solicitaciones[:, 0]
            self.axelines = self.axe.plot(nparray[:, 0], nparray[:, 1], 'o', label=label)
            print('solicitaciones 2D añadidas')
