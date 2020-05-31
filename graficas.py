import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from copy import copy

matplotlib.rc('axes', edgecolor='blue', facecolor='314c73', grid=True)
# matplotlib.rc('axes', labelcolor='#D4DCEE')  # compila bien pero no tiene efecto
matplotlib.rc('grid', linestyle='--', alpha=0.65)
matplotlib.rc('xtick', color='#D4DCEE')
matplotlib.rc('ytick', color='#D4DCEE')
matplotlib.rc('figure', facecolor='#303340')
matplotlib.rc('legend', facecolor='#426ca8')


def show():
    plt.legend()
    plt.show()


def add3d(nparray):
    plt.figure().add_subplot(projection='3d')
    plt.plot(nparray[:, 0], nparray[:, 1], nparray[:, 2])
    print('nparray3d')
    print(nparray)


def add2d(nparray, label):
    plt.plot(nparray[:, 0], nparray[:, 1], '+-', label=label)
    print('nparray2d')
    print(nparray)


def addPMM(pmm, bycolumns=True):
    if not (bycolumns):
        pmm = np.transpose(pmm)
    nparray = copy(pmm)
    nparray[:, 0] = pmm[:, 1]
    nparray[:, 1] = pmm[:, 2]
    nparray[:, 2] = pmm[:, 0]
    add3d(nparray[:, 0:3])


def addPM(pmm, bycolumns=True):
    if not (bycolumns):
        pmm = np.transpose(pmm)
    nparray = copy(pmm)
    nparray[:, 0] = (pmm[:, 1]**2+pmm[:, 2]**2)**0.5
    nparray[:, 1] = pmm[:, 0]
    add2d(nparray[:, 0:2])
