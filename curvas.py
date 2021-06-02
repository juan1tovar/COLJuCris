import Seccion as s
from math import sin, cos, tan, radians, log, pi
import numpy as np
import public_data as pd


vcc = pd.column_design_factors
il = pd.iteration_limit


# calcula la deformación unitara en la varilla a tracción
# para que de un determinado valor de fi.
def calc_es_with_fi(fi, fy, Es):
    return (fy/Es+(fi-vcc['ficomp'])
            * (vcc['deftrac']-fy/Es)/(vcc['fitracc']-vcc['ficomp']))


def calc_d(seccion, angulo):
    xyConc = seccion.cord_conc(angulo)
    xyRef = seccion.cord_ref(angulo)
    return max(xyConc[:, 1])-min(xyRef[:, 1])


def buscar_punto_with_fi(seccion, fi, asol, angulo=None, c=None):
    # fy = 420
    # Es = 200000
    es = calc_es_with_fi(fi, 420, 200000)
    x = cos(radians(asol))
    y = sin(radians(asol))
    pos = ((y/x < -1) or (y/x > 1)) and (y < 0)

    if (angulo is None) and (c is None):
        asol = s.calc_ang(x, y, pos)
        c = calc_d(seccion, asol)/(1+es/vcc['defConc'])  # d*defConc/(defConc+defsteel)
        res = s.resultante(seccion, asol, c)
        es = calc_es_with_fi(fi, res[4].fy, res[4].Es)
        angulo = s.calc_ang(res[1], res[2], pos)  # Primera aproximación del ángulo
        c = calc_d(seccion, angulo)/(1+es/vcc['defConc'])  # d*defConc/(defConc+defsteel)
    elif angulo is None:
        asol = s.calc_ang(x, y, pos)
        res = s.resultante(seccion, asol, c)
        es = calc_es_with_fi(fi, res[4].fy, res[4].Es)
        angulo = s.calc_ang(res[1], res[2], pos)  # Primera aproximación del ángulo
    else:
        c = calc_d(seccion, angulo)/(1+es/vcc['defConc'])

    res = s.fi_result(seccion, angulo, c)
    errA = abs(asol-s.calc_ang(res[1], res[2], positive=pos))
    if errA < il['errA']:
        return [angulo, c]

    es = calc_es_with_fi(fi, res[4].fy, res[4].Es)
    deltaC = 1
    deltaA = 1
    ies = 0
    while (deltaC > il['errC']) or (deltaA > il['errA']):
        a1 = angulo
        c1 = c
        a2 = s.buscar_ang(seccion, x, y, c1, a1)
        c2 = calc_d(seccion, a2)/(1+es/vcc['defConc'])
        a3 = s.buscar_ang(seccion, x, y, c2, a2)
        c3 = calc_d(seccion, a3)/(1+es/vcc['defConc'])
        deltaA = abs(a3-a2)
        deltaC = abs(c3-c2)
        # print('a1=%.4f a2=%.4f a3=%.4f' % (a1, a2, a3))
        # print('c1=%.4f c2=%.4f c3=%.4f' % (c1, c2, c3))
        # print('deltaA=%.4f deltaC=%.4f' % (deltaA, deltaC))
        if (abs(a3-a1) < il['errA']) and (abs(c3-c1) < il['errC']):
            angulo = (a3+a2)/2
            c = (c3+c2)/2
        else:
            angulo = a3
            c = c3
        ies = ies+1
        # breakpoint()
        print('ies=', ies)
    return [angulo, c]


def vertical(seccion, asol, sv=3):
    # sv = 3  # factor para suavisar la curva
    c, fipmm = [], []
    x, y = cos(radians(asol)), sin(radians(asol))
    c.append(1000000)
    rpmax = (vcc['ficomp']*s.Pnmax(seccion), 0, 0, vcc['ficomp'])  # ; print('rpmax')
    amax, cmax = s.buscar_punto(seccion, rpmax[0], x, y)  # ; print('amax, cmax')
    acom, ccom = buscar_punto_with_fi(seccion, vcc['ficomp'], asol, angulo=amax, c=0.7*cmax)
    atra, ctra = buscar_punto_with_fi(seccion, vcc['fitracc'], asol, angulo=amax, c=0.9*ccom)
    a0, c0 = s.buscar_punto(seccion, 0, x, y)  # ; print('a0 c0')
    c.extend(np.arange(cmax, ccom, -(cmax-ccom)/sv))
    if c[-1]-ccom < 0.0005:
        c.pop(-1)
    c.extend(np.arange(c[-1]-(c[-1]-ccom)/sv, ccom, -(c[-1]-ccom)/sv))
    if c[-1]-ccom < 0.0005:
        c.pop(-1)
    c.extend(np.arange(ccom, ctra, -(ccom-ctra)/2/sv))
    if c[-1]-ctra < 0.0005:
        c.pop(-1)
    if c0 > ctra:
        c0 = 2*ctra/3
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
    a = amax
    for i in range(1, sv*6):
        a = s.buscar_ang(seccion, x, y, c[i], a)
        fipmm.append(s.fi_result(seccion, a, c[i])[:4])
    fipmm.append(s.fi_result(seccion, 0, 0)[:4])
    return fipmm


# cálculo de n en la ecuación x^n+y^n=1 donde y^n=k*x^n
#  n = log(1/(1+k))/log(x)
#  k^c = 1+k
#  c = -log(x)/log(y/x)
def n_fatcircle(x, y):
    if x >= 1.0 or y >= 1.0 or x+y < 1.0:
        raise NameError('n_fatcircle: X o Y están fuera del rango')
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


def fiMn_cardinal(seccion, Pu, simetria='doble'):
    fiMn = [0, 0, 0, 0]
    fiMn[0] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, 1, 0))[1]
    fiMn[1] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, 0, 1))[2]
    if simetria == 'doble':
        fiMn[2] = -fiMn[0]
        fiMn[3] = -fiMn[1]
        return fiMn
    elif simetria == 'ejeX':
        fiMn[2] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, -1, 0))[1]
        fiMn[3] = -fiMn[1]
        return fiMn
    elif simetria == 'ejeY':
        fiMn[2] = -fiMn[0]
        fiMn[3] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, 0, -1))[2]
        return fiMn
    elif simetria == 'ninguna':
        fiMn[2] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, -1, 0))[1]
        fiMn[3] = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, 0, -1))[2]
        return fiMn
    else:
        raise NameError('fiMn_cardinal: simetría no definida')


def n_fatcircle_45(seccion, Pu, Mx, My):
    res = s.fi_result(seccion, *s.buscar_punto(seccion, Pu, Mx, My))
    # print('res=', res)
    return log(0.5)/log(res[1]/Mx)


# div   : cantidad de divisiones en un cuadrante de curva
# M[0]  : fi*Mn cuando el angulo de reacción es 0
# M[1]  : fi*Mn cuando el angulo de reacción es 90
# M[2]  : fi*Mn cuando el angulo de reacción es 180
# M[3]  : fi*Mn cuando el angulo de reacción es 270
def fatcircle(seccion, div, angi, angf, Pu, M=[]):
    if M == []:
        M = fiMn_cardinal(seccion, Pu)
    # print(M)
    angi = angi % 360  # angulo de reacción donde inicia la curva
    angf = angf % 360  # angulo de reacción donde finaliza la curva
    if angf-angi <= 0:
        angf = angf + 360
    a = radians(angi)
    da = radians(90/div)
    angf = radians(angf)
    fiMn = []

# aún se puede mejorar la distribución de los puntos para que se concentre
# en los extremos y en el medio.

# NO FUNCIONA CORRECTAMENTE SI EL ÁNGULO INICIAL ES MENOR A 0 se requeriría un switch

    n = n_fatcircle_45(seccion, Pu, M[0], M[1])
    while a % (2*pi) <= pi/2 and a <= angf:
        x = 1/(1+abs(tan(a)))**(1/n)
        y = (1-x**n)**(1/n)
        # print('fatcircle a=%.2f x=%.3f y=%.4f n=%.2f' % (degrees(a), x, y, n))
        fiMn.append([x*M[0], y*M[1]])
        a = a+da

    n = n_fatcircle_45(seccion, Pu, M[2], M[1])
    while a % (2*pi) <= pi and a <= angf:
        x = 1/(1+abs(tan(a)))**(1/n)
        y = (1-x**n)**(1/n)
        # print('fatcircle a=%.2f x=%.3f y=%.4f n=%.2f' % (degrees(a), x, y, n))
        fiMn.append([x*M[2], y*M[1]])
        a = a+da

    n = n_fatcircle_45(seccion, Pu, M[2], M[3])
    while a % (2*pi) <= 3*pi/2 and a <= angf:
        x = 1/(1+abs(tan(a)))**(1/n)
        y = (1-x**n)**(1/n)
        # print('fatcircle a=%.2f x=%.3f y=%.4f n=%.2f' % (degrees(a), x, y, n))
        fiMn.append([x*M[2], y*M[3]])
        a = a+da

    n = n_fatcircle_45(seccion, Pu, M[0], M[3])
    while a % (2*pi) <= 2*pi and a <= angf:
        x = 1/(1+abs(tan(a)))**(1/n)
        y = (1-x**n)**(1/n)
        # print('fatcircle a=%.2f x=%.3f y=%.4f n=%.2f' % (degrees(a), x, y, n))
        fiMn.append([x*M[0], y*M[3]])
        a = a+da
    return fiMn


def h_ang_sol(seccion, div, angi, angf, Pu):
    angi = angi % 360  # angulo de reacción donde inicia la curva
    angf = angf % 360  # angulo de reacción donde finaliza la curva
    if angf-angi <= 0:
        angf = angf + 360
    ar = radians(angi)
    da = radians(90/div)
    angf = radians(angf)
    fiMn = []
    a, c = s.buscar_punto(seccion, Pu, cos(ar), sin(ar))
    while ar <= angf:
        # print('h_ang_sol a=%.2f ar=%.2f c=%.2f' % (a, degrees(ar), c))
        # breakpoint()
        fiMn.append(s.fi_result(seccion, a, c)[1:3])
        ar = ar+da
        a, c = s.buscar_punto(seccion, Pu, cos(ar), sin(ar), c=c)
    return fiMn


def h_ang_ejeN(seccion, div, angi, angf, Pu):
    angi = angi % 360  # angulo de reacción donde inicia la curva
    angf = angf % 360  # angulo de reacción donde finaliza la curva
    if angf-angi <= 0:
        angf = angf + 360
    a = angi
    da = 90/div
    angf = angf
    fiMn = []
    c = s.buscar_c(seccion, Pu, a)
    while a <= angf:
        fiMn.append(s.fi_result(seccion, a, c)[1:3])
        a = a+da
        c = s.buscar_c(seccion, Pu, a, c)
    return fiMn


def horizontal(seccion, Pu, metodo='fatcircle', div=10, angi=0, angf=360):
    MxMy = {
        'fatcircle': fatcircle,
        'ang_sol': h_ang_sol,
        'ang_eje': h_ang_ejeN
    }.get(metodo)
    if MxMy is None:
        raise NameError('Método de curva horizontal no definida')
    return MxMy(seccion, div, angi, angf, Pu)
