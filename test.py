# !D:\Coding\Proyectos\COLJuCris\venv\Scripts\python.exe
import sys
import debugexc
import numpy as np
import Seccion as s
# from Seccion import calc_ang
import curvas
import graficas
# import gui
import public_data as pd
# import design_col as dc


sys.excepthook = debugexc.info
vcc = pd.column_design_factors

Fc = 21
X = .4
Y = .25
angulo = -10.18
# c = 0.2
Pu = -53
PxX = -3.04
PxY = -6.27
recu = 0.04
sec = s.seccion(Fc, X, Y)
np.set_printoptions(precision=3, suppress=True,)
print("dimensión x=%.3f m dimensión y=%.3f"
      % (sec.x, sec.y))
print("f'c=%.1f MPa, E=%.2f MPa, defc=%.3f"
      % (sec.fc, sec.Ec_MPa(), vcc['defConc']))
print("coordenas de las equinas a 0°")
print(sec.cord_conc(0))
print("coordenas de las equinas a angulo")
print(sec.cord_conc(angulo))
# breakpoint()
var1 = s.seccion_varilla('#6', 420, 200000)
var2 = s.seccion_varilla('#6', 420, 200000)
varf = s.seccion_varilla('#3', 420, 200000)
sec.set_ref_sim(2, 3, var1, var2, varf,
                recu, opcion='esquinero')
# print('Refuerzo XY, Varilla, Área')
# print(sec.refXY, sec.varillas)
# print("coordenas de las varillas rotadas angulo=", angulo)
# print(sec.cord_ref(angulo))
# print(var[0:2] for var in sXxY.Vari)
# print(np.array(sXxY.iteracion(angulo, 0.52)))
# sec.areas_varillas()
# sec.fy = 420*np.ones(len(sXxY.AreasRef))
# angulo, c = s.buscar_punto(sec, Pu, PxX, PxY)
# print(angulo, c)
# breakpoint()
res = s.result_conc(sec, -4.448965545, 0.306521003568861)
print('res', res)
# fires = sec.fi_result(angulo, c)
# print(fires)
# print(calc_ang(fires[1], fires[2]))
# print(calc_d(sec, angulo))
graf3D = graficas.grafica(is_3D=True)
# ventana = interface.tkwindow(graf)
pmm = np.array(curvas.vertical(sec, angulo))
graf3D.addPMM(pmm, str(angulo)+'°')
graf3D.show()
# grafver = graficas.grafica(is_3D=False)
# grafver.addPM(pmm, str(angulo)+'°')
# grafver.show()
# grafhor = graficas.grafica(is_3D=False)
# parlem = np.array(curvas.horizontal(sec, Pu, metodo='fatcircle'))
# grafhor.add2D(parlem, "parlem")
# solici = np.array(curvas.horizontal(sec, Pu, metodo='ang_sol'))
# grafhor.add2D(solici, "ang sol")
# column = np.array(curvas.horizontal(sec, Pu, metodo='ang_eje'))
# grafhor.add2D(column, "ang_col")
# grafhor.show()

# print(*grafhor.axe.lines)
# grafhor.axe.lines[0].remove()
# print(*grafhor.axe.lines)
# grafhor.show()

# breakpoint()
# ventana.mainloop()

# sol = dc.solicitation(100, 10, 20, 5, 90, 40)
# print("solicitaciones: ", sol.str())
# a, c = s.buscar_punto(sec, sol.P, sol.Mx, sol.My)
# print(f"punto encontrado: angulo={a:.2f}° c={c:.3f}m")
# res = s.fi_result(sec, a, c)
# print(f"res: φPn={res[0]:.1f}, φMnx={res[1]:.1f}, φMny={res[2]:.1f}")
# # print("P=%.1f, Mx=%.1f, My=%.1f (porcentaje operator)"
# #       % (res[0], res[1], res[2]))
# # print("P={0:.1f}, Mx={1:.1f}, My={2:.1f} (format() method)".format(*s.fi_result(sec, a, c)))
# ind = dc.indice_flco(sec, sol)
# print(f"índice={ind:.3f}")
# # print("índice={:.1f} (format() method)".format(ind))


print('ok')
