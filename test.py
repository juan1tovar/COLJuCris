import sys
import debugexc
import numpy as np
import Seccion as s
# from Seccion import calc_ang
import curvas
import graficas
# import interface
import public_data as pd


sys.excepthook = debugexc.info
vcc = pd.column_design_factors

Fc = 28
X = .6
Y = .3
angulo = 45
# c = 0.2
Pu = 100
recu = 0.04
sXxY = s.seccion(Fc, X, Y)
np.set_printoptions(precision=3, suppress=True)
print("dimensión x=%.3f m dimensión y=%.3f"
      % (sXxY.x, sXxY.y))
print("f'c=%.1f MPa, E=%.2f MPa, defc=%.3f"
      % (sXxY.fc, sXxY.Ec_MPa(), vcc['defConc']))
# print("coordenas de las equinas a 0°")
# print(sXxY.cord_conc(0))class
# print("coordenas de las equinas a angulo")class
# print(sXxY.cord_conc(angulo))
# breakpoint()
var1 = s.seccion_varilla('#5', 420, 200000)
var2 = s.seccion_varilla('#5', 420, 200000)
varf = s.seccion_varilla('#3', 420, 200000)
sXxY.set_ref_sim(4, 4, var1, var2, varf,
                 recu, opcion='esquinero')
# print('Refuerzo XY, Varilla, Área')
# print(sXxY.refXY, sXxY.varillas)
print("coordenas de las varillas rotadas angulo=", angulo)
print(sXxY.cord_ref(angulo))
# print(var[0:2] for var in sXxY.Vari)
# print(np.array(sXxY.iteracion(angulo, 0.52)))
# sXxY.areas_varillas()
# sXxY.fy = 420*np.ones(len(sXxY.AreasRef))
# angulo, c = sXxY.buscar_punto(1775, 115.8, -467.2)
# print(angulo, c)
# res = sXxY.resultante(angulo, c)
# print(res)
# fires = sXxY.fi_result(angulo, c)
# print(fires)
# print(calc_ang(fires[1], fires[2]))
# print(calc_d(sXxY, angulo))
graf3D = graficas.grafica(is_3D=True)
# ventana = interface.tkwindow(graf)
pmm = np.array(curvas.vertical(sXxY, angulo))
graf3D.addPMM(pmm, str(angulo)+'°')
graf3D.show()
grafver = graficas.grafica(is_3D=False)
grafver.addPM(pmm, str(angulo)+'°')
grafver.show()
grafhor = graficas.grafica(is_3D=False)
parlem = np.array(curvas.horizontal(sXxY, Pu, metodo='fatcircle'))
grafhor.add2D(parlem, "parlem")
solici = np.array(curvas.horizontal(sXxY, Pu, metodo='ang_sol'))
grafhor.add2D(solici, "ang sol")
column = np.array(curvas.horizontal(sXxY, Pu, metodo='ang_col'))
grafhor.add2D(column, "ang_col")
grafhor.show()
# ventana.mainloop()

print('ok')
