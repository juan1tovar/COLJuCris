import Seccion as sec
from math import sin, cos, radians, log
import numpy as np


# def calc_es_with_fi(s, fi):
#     return min(s.fy)/s.Es+(fi-s.ficomp)*(s.deftrac-min(s.fy)/s.Es)/(s.fitracc-s.ficomp)


def calc_d(s, angulo):
    xyConc = s.cord_conc(angulo)
    xyRef = s.cord_ref(angulo)
    return max(xyConc[:, 1])-min(xyRef[:, 1])


def buscar_punto_with_fi(s, fi, asol, angulo=None, c=None):
    es = max(s.fy)/s.Es+(fi-s.ficomp)*(s.deftrac-max(s.fy)/s.Es)/(s.fitracc-s.ficomp)
    x = cos(radians(asol))
    y = sin(radians(asol))
    pos = ((y/x < -1) or (y/x > 1)) and (y < 0)

    if (angulo is None) and (c is None):
        angulo = (-asol+90) % 360  # Primera aproximación del ángulo
        c = calc_d(s, angulo)/(1+es/s.defConc)  # d*defConc/(defConc+defsteel)
    elif angulo is None:
        angulo = (-asol+90) % 360  # Primera aproximación del ángulo
    else:
        c = calc_d(s, angulo)/(1+es/s.defConc)

    res = s.fi_result(angulo, c)
    errA = abs(asol-sec.calc_ang(res[1], res[2], positive=pos))
    if errA < 0.001:
        return [angulo, c]

    deltaC = 1
    deltaA = 1
    ies = 0
    while (deltaC > 0.0001) or (deltaA > 0.001):
        a1 = angulo
        c1 = c
        a2 = s.buscar_ang(x, y, c1, ac1=a1)
        c2 = calc_d(s, a2)/(1+es/s.defConc)
        a3 = s.buscar_ang(x, y, c2, ac1=a2)
        c3 = calc_d(s, a3)/(1+es/s.defConc)
        deltaA = abs(a3-a2)
        deltaC = abs(c3-c2)
        # print(a1, a2, a3)
        # print(c1, c2, c3)
        # print(deltaA, deltaC)
        if (abs(a3-a1) < 0.01) and (abs(c3-c1) < 0.0001):
            angulo = (a3+a2)/2
            c = (c3+c2)/2
        else:
            angulo = a3
            c = c3
        ies = ies+1
    # breakpoint()
    print('ies=', ies)
    return [angulo, c]


def vertical(s, asol, sv=3):
    # sv = 3  # factor para suavisar la curva
    c, fipmm = [], []
    x, y = cos(radians(asol)), sin(radians(asol))
    c.append(1000000)
    rpmax = (s.ficomp*s.Pnmax(), 0, 0, s.ficomp)  # ; print('rpmax')
    amax, cmax = s.buscar_punto(rpmax[0], x, y)  # ; print('amax, cmax')
    acom, ccom = buscar_punto_with_fi(s, s.ficomp, asol, angulo=amax, c=0.7*cmax)
    atra, ctra = buscar_punto_with_fi(s, s.fitracc, asol, angulo=amax, c=0.9*ccom)
    a0, c0 = s.buscar_punto(0, x, y)  # ; print('a0 c0')
    c.extend(np.arange(cmax, ccom, -(cmax-ccom)/sv))
    if c[-1]-ccom < 0.0005:
        c.pop(-1)
    c.extend(np.arange(c[-1]-(c[-1]-ccom)/sv, ccom, -(c[-1]-ccom)/sv))
    if c[-1]-ccom < 0.0005:
        c.pop(-1)
    c.extend(np.arange(ccom, ctra, -(ccom-ctra)/2/sv))
    if c[-1]-ctra < 0.0005:
        c.pop(-1)
    c.extend(np.arange(ctra, c0, -(ctra-c0)/sv))
    if c[-1]-c0 < 0.0005:
        c.pop(-1)
    c.extend(np.arange(c0, 0, -(c0-0)/sv))
    if c[-1] < 0.005:
        c.pop(-1)
    c.append(0)
    fipmm.append(rpmax)
    # se podría reducir la cantidad de iteraciones de buscar_ang si
    # cada iteración inicia desde el último ángulo encontrado
    print('c')
    print(c)
    for i in range(1, sv*6):
        fipmm.append(s.fi_result(s.buscar_ang(x, y, c[i]), c[i]))
    fipmm.append(s.fi_result(0, 0))
    return fipmm


# cálculo de n en la ecuación x^n+y^n=1 donde y^n=k*x^n
#  n = log(1/(1+k))/log(x)
#  k^c = 1+k
#  c = -log(x)/log(y/x)
def n_fatcircle(x, y):
    if x >= 1.0 or y >= 1.0 or x+y < 1.0:
        raise 'n_fatcircle: X o Y están fuera del rango'
    if x == y:
        return log(0.5)/log(x)
    c = -log(x)/log(y/x)
    k = 1
    # def f(k): return k**c-k-1
    # def df(k): return c*k**(c-1)-1
    # def g(k): return k-f(k)/df(k)
    # while abs(f(k)) > 0.001:
    #     k = g(k)
    while abs(k**c-k-1) > 0.001:
        k = k-(k**c-k-1)/(c*k**(c-1)-1)
    return log(1/(1+k))/log(x)


def fiMn_cardinal(s, Pu, simetria='doble'):
    fiMn = [0, 0, 0, 0]
    fiMn[0] = s.fi_result(*s.buscar_punto(Pu, 1, 0))[1]
    fiMn[1] = s.fi_result(s.buscar_punto(Pu, 0, 1))[2]
    if simetria == 'doble':
        fiMn[2] = -fiMn[0]
        fiMn[3] = -fiMn[1]
        return fiMn
    elif simetria == 'ejeX':
        fiMn[2] = s.fi_result(s.buscar_punto(Pu, -1, 0))[1]
        fiMn[3] = -fiMn[1]
        return fiMn
    elif simetria == 'ejeY':
        fiMn[2] = -fiMn[0]
        fiMn[3] = s.fi_result(s.buscar_punto(Pu, 0, -1))[2]
        return fiMn
    elif simetria == 'ninguna':
        fiMn[2] = s.fi_result(s.buscar_punto(Pu, -1, 0))[1]
        fiMn[3] = s.fi_result(s.buscar_punto(Pu, 0, -1))[2]
        return fiMn
    else:
        raise 'fiMn_cardinal: simetría no definida'


def n_fatcircle_45(s, Pu, Mx, My):
    a, c = s.buscar_punto(Pu, Mx, My)
    res = s.resultante(a, c)
    return log(0.5)/log(res[1]/Mx)


# div   : cantidad de divisiones en un cuadrante de curva
# M[0]  : fi*Mn cuando el angulo de reacción es 0
# M[1]  : fi*Mn cuando el angulo de reacción es 90
# M[2]  : fi*Mn cuando el angulo de reacción es 180
# M[3]  : fi*Mn cuando el angulo de reacción es 270
def fatcircle(s, div, angi, angf, Pu, M):
    angi = angi % 360  # angulo de reacción donde inicia la curva
    angf = angf % 360  # angulo de reacción donde finaliza la curva
    if angf < angi:
        raise 'fatcircle: el ángulo final es menor que el inicial'
    a = angi
    da = 90/div
    fiMn = []
    while a < 90 and a < angf:
        n = n_fatcircle_45(s, Pu, M[0], M[1])
        x = cos(a)
        y = (1-x**n)**(1/n)
        fiMn.append([x*M[0], y*M[1]])
        a = a+da
    while a < 180 and a < angf:
        n = n_fatcircle_45(s, Pu, M[2], M[1])
        x = cos(a)
        y = (1-x**n)**(1/n)
        fiMn.append([x*M[2], y*M[1]])
        a = a+da
    while a < 270 and a < angf:
        n = n_fatcircle_45(s, Pu, M[2], M[3])
        x = cos(a)
        y = (1-x**n)**(1/n)
        fiMn.append([x*M[2], y*M[3]])
        a = a+da
    while a < 360 and a < angf:
        n = n_fatcircle_45(s, Pu, M[0], M[3])
        x = cos(a)
        y = (1-x**n)**(1/n)
        fiMn.append([x*M[0], y*M[3]])
        a = a+da


def horizontal(s, Pu, metodo='fatcircle', div=10, angi=0, angf=360):
    MxMy = {
        'fatcircle': fatcircle,
        # 'exacto': h_exacta
    }.get(metodo)
    if MxMy is None:
        raise 'Método de curva horizontal no definida'
    return MxMy(s, div, angi, angf, Pu, fiMn_cardinal(s, Pu))
