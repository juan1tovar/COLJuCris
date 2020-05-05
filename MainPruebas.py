import sys
import debugexc
import numpy as np
from Seccion import seccion
# from Seccion import calc_ang
import curvas
import graficas as graf
# import curvas
sys.excepthook = debugexc.info
Fc = 28
X = 0.4
Y = 0.6
angulo = 99
c = 0.2
sXxY = seccion(Fc, X, Y)
np.set_printoptions(precision=3, suppress=True)
print("dimensión x=%.3f m dimensión y=%.3f"
      % (sXxY.x, sXxY.y))
print("f'c=%.1f MPa, E=%.2f MPa, defc=%.3f"
      % (sXxY.fc, sXxY.Ec_MPa(), sXxY.defConc))
# print("coordenas de las equinas a 0°")
# print(sXxY.cord_conc(0))class
# print("coordenas de las equinas a angulo")
# print(sXxY.cord_conc(angulo))
sXxY.set_ref_sim(5, 5, '#7', '#7', opcion='esquinero')
# print('Refuerzo XY, Varilla, Área')
# print(sXxY.refXY, sXxY.varillas)
print("coordenas de las varillas rotadas angulo")
print(sXxY.cord_ref(angulo))
# print(var[0:2] for var in sXxY.Vari)
# print(np.array(sXxY.iteracion(angulo, 0.52)))
# sXxY.areas_varillas()
# sXxY.fy = 420*np.ones(len(sXxY.AreasRef))
# angulo, c = sXxY.buscar_punto(1775, 115.8, -467.2)
# print(angulo, c)
res = sXxY.resultante(angulo, c)
print(res)
# fires = sXxY.fi_result(angulo, c)
# print(fires)
# print(calc_ang(fires[1], fires[2]))
# print(calc_d(sXxY, angulo))
# pmm = np.array(curvas.vertical(sXxY, -279))
# print('pmm')
# print(pmm)
# graf.plotPMM(pmm)
# graf.plotPMM(np.transpose(pmm), bycolumns=False)
# graf.plotPM(pmm)
parlem = np.array(curvas.horizontal(sXxY, 1775))
print(parlem)
graf.plot2d(parlem)
