import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from copy import copy

matplotlib.rc('axes', edgecolor='blue', facecolor='#303340', grid=True)
# matplotlib.rc('axes', labelcolor='#D4DCEE')  # compila bien pero no tiene efecto
matplotlib.rc('grid', linestyle='--', alpha=0.65)
matplotlib.rc('xtick', color='#D4DCEE')
matplotlib.rc('ytick', color='#D4DCEE')
matplotlib.rc('figure', facecolor='#303340')


def plot3d(nparray):
    plt.figure().add_subplot(projection='3d')
    plt.plot(nparray[:, 0], nparray[:, 1], nparray[:, 2])
    print('nparray3d')
    print(nparray)
    plt.show()


def plot2d(nparray):
    plt.plot(nparray[:, 0], nparray[:, 1], '+-')
    print('nparray2d')
    print(nparray)
    plt.show()


def plotPMM(pmm, bycolumns=True):
    if not (bycolumns):
        pmm = np.transpose(pmm)
    nparray = copy(pmm)
    nparray[:, 0] = pmm[:, 1]
    nparray[:, 1] = pmm[:, 2]
    nparray[:, 2] = pmm[:, 0]
    plot3d(nparray[:, 0:3])


def plotPM(pmm, bycolumns=True):
    if not (bycolumns):
        pmm = np.transpose(pmm)
    nparray = copy(pmm)
    nparray[:, 0] = (pmm[:, 1]**2+pmm[:, 2]**2)**0.5
    nparray[:, 1] = pmm[:, 0]
    plot2d(nparray[:, 0:2])
