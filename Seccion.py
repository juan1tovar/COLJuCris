import numpy as np
# import math
from math import radians, sin, cos, copysign, atan, tan, degrees
from operator import itemgetter
import public_data as pd
np.set_printoptions(precision=3, suppress=True)


# def np_err_handler(type, flag):
#     print("Floating point error (%s), with flag %s" % (type, flag))
#     print('form more info search for numpy.seterrcall')
#     breakpoint()


# np.seterrcall(np_err_handler)
# np.seterr(all='call')


def rotar(coordenadas, grados, invertir=False):
    # la operación de módulo me aproxima los números negativos
    # por eso me tocó meter todo en esa fórmula
    a = radians(-grados)
    t = np.array([[cos(a), -sin(a)],
                  [sin(a), cos(a)]])
    if invertir:
        # breakpoint()
        return t.dot(coordenadas)
    else:
        # breakpoint()
        return coordenadas.dot(t)


def calc_ang(x, y, positive=False):
    ang = None
    if x == 0:
        if y > 0:
            ang = 90
        else:
            ang = -90
    else:
        ang = degrees(atan(y/x))

    if x < 0:
        ang = ang + 180

    if positive and (ang < 90):
        ang = ang + 360
    return ang


class seccion_varilla:

    def __init__(self, etiqueta, fy, Es):
        self.etiq = etiqueta
        self.diam = pd.steel_bars[etiqueta][0]/1000
        self.area = pd.steel_bars[etiqueta][1]/1000/1000
        self.fy = fy  # MPa
        self.Es = Es  # MPa


class refuerzo:

    def __init__(self):
        self.refXY = []
        self.varillas = []
        self.refActual = ''  # Es un flag para saber si ya se asignó un refuerzo
        # self.rec = recubrimiento

    # tablaRefuerzo = public_data.varillas

    def cord_ref(self, angulo):
        if self.refXY == []:
            raise NameError('No se han ingresado varillas de refuerzo')
        return rotar(self.refXY, angulo)

    def As(self):
        if self.refXY == []:
            raise NameError('No se han ingresado varillas de refuerzo')

        As = 0
        for var in self.varillas:
            As = As + var.area
        return As

    # def fy_max(self):
    #     raise NameError("Sin implementar")
    #     self.fy = 420*np.ones(len(self.AreasRef))

    def set_ref_sim(self, NVarX, NVarY, v_esquina, v_alterna, v_fleje, recub,
                    opcion='intercalado'):
        func = getattr(self, opcion, 'OpcInvalida')

        if func == 'OpcInvalida':
            raise NameError('No existe esa opción de refuerzo simétrico')

        func(NVarX, NVarY, v_esquina, v_alterna, v_fleje, recub)
        self.refActual = opcion

    def intercalado(self, NVarX, NVarY, v_esquina, v_alterna, v_fleje, recub):
        DE = v_esquina.diam
        DA = v_alterna.diam
        DF = v_fleje.diam
        Xmax = self.x/2-recub-DF
        Ymax = self.y/2-recub-DF
        sepx = (Xmax*2-DE)/(NVarX-1)
        if sepx < 0:
            raise NameError("Intercalado: sepx<0")
        sepy = (Ymax*2-DE)/(NVarY-1)
        if sepy < 0:
            raise NameError("Intercalado: sepy<0")
        varx = [list(range(int(NVarX/2)+NVarX % 2)),
                list(range(int(NVarX/2)+NVarX % 2)),
                list(range(int(NVarX/2)+NVarX % 2)),
                list(range(int(NVarX/2)+NVarX % 2))]

        vary = [list(range(int(NVarY/2)+NVarY % 2)),
                list(range(int(NVarY/2)+NVarY % 2)),
                list(range(int(NVarY/2)+NVarY % 2)),
                list(range(int(NVarY/2)+NVarY % 2))]

        for n in range(0, int(NVarX/2)+NVarX % 2, 2):
            x = Xmax-DE/2-sepx*n
            y = Ymax-DE/2
            varx[0][n] = [x, y, v_esquina]
            varx[1][n] = [-x, y, v_esquina]
            varx[2][n] = [-x, -y, v_esquina]
            varx[3][n] = [x, -y, v_esquina]
            if n+1 == int(NVarX/2)+NVarX % 2:
                break
            x = Xmax-DE/2-sepx*(n+1)
            y = Ymax-DA/2
            varx[0][n+1] = [x, y, v_alterna]
            varx[1][n+1] = [-x, y, v_alterna]
            varx[2][n+1] = [-x, -y, v_alterna]
            varx[3][n+1] = [x, -y, v_alterna]

        for n in range(0, int(NVarY/2)+NVarY % 2, 2):
            x = Xmax-DE/2
            y = Ymax-DE/2-sepy*n
            vary[0][n] = [x, y, v_esquina]
            vary[1][n] = [-x, y, v_esquina]
            vary[2][n] = [-x, -y, v_esquina]
            vary[3][n] = [x, -y, v_esquina]
            if n+1 == int(NVarY/2)+NVarY % 2:
                break
            x = Xmax-DA/2
            y = Ymax-DE/2-sepy*(n+1)
            vary[0][n+1] = [x, y, v_alterna]
            vary[1][n+1] = [-x, y, v_alterna]
            vary[2][n+1] = [-x, -y, v_alterna]
            vary[3][n+1] = [x, -y, v_alterna]

        ReSi = []
        ReSi.extend(vary[3][:])
        vary[0].reverse()
        ReSi.extend(vary[0][NVarY % 2:])
        ReSi.extend(varx[0][1:])
        varx[1].reverse()
        ReSi.extend(varx[1][NVarX % 2:])
        ReSi.extend(vary[1][1:])
        vary[2].reverse()
        ReSi.extend(vary[2][NVarY % 2:])
        ReSi.extend(varx[2][1:])
        varx[3].reverse()
        ReSi.extend(varx[3][NVarX % 2:-1])
        self.refXY = np.array([var[0:2] for var in ReSi])
        self.varillas = [var[2] for var in ReSi]

    def esquinero(self, NVarX, NVarY, v_esquina, v_alterna, v_fleje, recub):
        DE = v_esquina.diam
        DA = v_alterna.diam
        DF = v_fleje.diam
        Xmax = self.x/2-recub-DF
        Ymax = self.y/2-recub-DF
        sepx = (Xmax*2-DE)/(NVarX-1)
        if sepx < 0:
            raise NameError("Esquinero: sepx<0")
        sepy = (Ymax*2-DE)/(NVarY-1)
        if sepy < 0:
            raise NameError("Esquinero: sepy<0")
        ReSi = []
        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([x, -y, v_esquina])
        x = Xmax-DA/2

        for i in range(1, NVarY-1):
            y = Ymax-DE/2-sepy*i
            ReSi.append([x, -y, v_alterna])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([x, y, v_esquina])
        y = Ymax-DA/2

        for i in range(1, NVarX-1):
            x = Xmax-DE/2-sepx*i
            ReSi.append([x, y, v_alterna])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([-x, y, v_esquina])
        x = Xmax-DA/2

        for i in range(1, NVarY-1):
            y = Ymax-DE/2-sepy*i
            ReSi.append([-x, y, v_alterna])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([-x, -y, v_esquina])
        y = Ymax-DA/2

        for i in range(1, NVarX-1):
            x = Xmax-DE/2-sepx*i
            ReSi.append([-x, -y, v_alterna])

        self.refXY = np.array([var[0:2] for var in ReSi])
        self.varillas = [var[2] for var in ReSi]


class concreto:

    def __init__(self, fc):
        self.fc = fc  # MPa

    def Ec_MPa(self):
        return 4700*self.fc**0.5


class seccion(concreto, refuerzo):

    def __init__(self, fc, LadoX, LadoY):
        concreto.__init__(self, fc)
        refuerzo.__init__(self)
        self.x = LadoX  # metros
        self.y = LadoY  # metros

    def psi(self):
        return atan(self.y/self.x)

    def Ag(self):
        return self.x*self.y

    def esquinas(self):
        return np.array([[0+self.x/2, 0-self.y/2],
                         [0+self.x/2, 0+self.y/2],
                         [0-self.x/2, 0+self.y/2],
                         [0-self.x/2, 0-self.y/2]])

    def cord_conc(self, angulo):
        return rotar(self.esquinas(), angulo)


def Pnmax(seccion, factor=0.75):
    Ag = seccion.Ag()
    AsFy = 0
    for var in seccion.varillas:
        AsFy = AsFy+var.area*var.fy*1000
    return factor*(0.85*(Ag-seccion.As())*seccion.fc*1000+AsFy)
    # 0.75 es un factor que depende del código o del material


def result_conc(seccion, angulo, ejeN):
    cor = sorted(seccion.cord_conc(angulo), key=itemgetter(0))
    cor = sorted(cor, key=itemgetter(1))
    ymax = cor[3][1]
    yEje = ymax-ejeN
    if ejeN == 0:
        return [0, 0, 0, yEje]
    betac = ejeN*max(min(0.85, 0.85-0.05*(seccion.fc-28)/7), 0.65)

    # htri : altura (height) del triángulo
    # btri : base del triángulo
    # hpar : altura del paralelogramo
    # bpar : base del paralelogramo
    # htra : altura del trapecio
    htri, hpar, htra, btri, bpar, btra, ta = 0, 0, 0, 0, 0, 0, 0

    ht = ymax-cor[2][1]
    if abs(angulo) % 90 == 0 or ht == 0:
        ht = 1e-10
        # htri, btri, htra, btra = 0, 0, 0, 0, 0
        if (abs(angulo) % 180 < 1) or (abs(angulo) % 180 > 179):
            ta = 0.00001
            hpar = min(betac, seccion.y)
            bpar = seccion.x
        else:
            ta = 1000000
            hpar = min(betac, seccion.x)
            bpar = seccion.y
    else:
        htri = min(ht, betac)
        hpar = min(cor[2][1]-cor[1][1], betac-htri)
        htra = min(ht, betac-hpar-htri)
        angbase = abs(angulo) % 180.0  # Este angulo se usa para calcular el bloque de Witney
        if angbase > 90:
            angbase = 180-angbase
        if radians(angbase) < seccion.psi():
            bpar = seccion.x/cos(radians(angbase))
            ta = tan(radians(angbase))
        else:
            bpar = seccion.y/cos(radians(90-angbase))
            ta = tan(radians(90-angbase))
        btri = bpar/ht*htri
        btra = bpar/ht*(ht-htra)

    # print('ta=',ta)
    # if ta == 0:
    #     ta = 0.00001

    A = [btri*htri/2,
         bpar*hpar,
         (bpar+btra)/2*htra]
    x3 = cor[3][0]
    x2 = cor[2][0]
    x1 = cor[1][0]
    X = [x3+htri/ht*(x2-x3+(x1-x3)*ht/(ymax-cor[1][1]))/3,
         copysign(bpar/2, x1)+x2+(cor[0][0]-x2)*hpar/2/(cor[2][1]-cor[0][1]),
         x1-copysign(1, x1)*(btra**2/2+htra*(btra/ta+bpar/2*ta)
                             + htra**2/3*(1/ta**2-ta**2))/((bpar+btra)/2)]
    Y = [ymax-2*htri/3,
         cor[2][1]-hpar/2,
         cor[1][1]-htra/3*(bpar+2*btra)/(bpar+btra)]
    Ag = sum(A)
    # breakpoint()
    A = np.array(A)
    Fuerza = Ag*0.85*seccion.fc*1000
    xyro = np.array([[sum(A*np.array(X))/Ag],
                     [sum(A*np.array(Y))/Ag]])
    xy = rotar(xyro, angulo, invertir=True)
    # print('Refuerzo: Fuerza=%.1f FxX=%.1f FxY=%.1f',
    #       (Fuerza, Fuerza*xy[0, 0], Fuerza*xy[1, 0]))
    # breakpoint()
    return [Fuerza, Fuerza*xy[0, 0], Fuerza*xy[1, 0], yEje]


def result_ref(seccion, angulo, ejeN, yEje, defConc=0.003, alfa=1.0):
    # defConc   Ɛ deformación unitaria de falla del concreto
    # alfa      α factor de amplificación por sobreresistencia del acero
    cor = seccion.cord_ref(angulo)
    defmin = 1
    varilla = None
    fuerza = []
    if ejeN == 0:
        ejeN = 0.0000001
    for xy, var in zip(cor, seccion.varillas):
        du = (xy[1]-yEje)*defConc/ejeN
        if du < defmin:
            defmin = du
            varilla = var
        if du > 0:
            fuerza.append(min(var.fy*alfa, var.Es*du)*var.area*1000)
        else:
            fuerza.append(max(-var.fy*alfa, var.Es*du)*var.area*1000)
    fuerza = np.array(fuerza)
    Fuerza = np.sum(fuerza)
    FxX = fuerza.dot(seccion.refXY[:, 0])
    FxY = fuerza.dot(seccion.refXY[:, 1])
    # print('Refuerzo: Fuerza=%.1f FxX=%.1f FxY=%.1f', (Fuerza, FxX, FxY))
    return [Fuerza, FxX, FxY, defmin, varilla]


vcc = pd.column_design_factors


def resultante(seccion, angulo, ejeN, defConc=vcc['defConc'], alfa=vcc['alfa']):
    # defConc   Ɛ deformación unitaria de falla del concreto
    # alfa      α factor de amplificación por sobreresistencia del acero
    if seccion.refXY == [] or seccion.refXY == []:
        raise NameError('No se han ingresado varillas de refuerzo')

    if ejeN < 0:
        raise NameError("c<0, profundidad del eje neutro no válida")
    # xyConc = sec.cord_conc(angulo)
    # if ejeN >= max(xyConc[:, 1])*2*sec.defConc/(sec.defConc-max(sec.fy)/sec.Es):
    #     return [sec.Pnmax(), 0, 0, max(sec.fy)/sec.Es]
    # if abs(angulo) < 0.0001:
        # breakpoint()

    rConc = result_conc(seccion, angulo, ejeN)
    rRefu = result_ref(seccion, angulo, ejeN, rConc[3], defConc, alfa)
    F = rConc[0]+rRefu[0]
    FxX = rConc[1]+rRefu[1]
    FxY = rConc[2]+rRefu[2]
    # breakpoint()
    return [F, FxX, FxY, rRefu[3], rRefu[4]]


# Este método se debe ampliar para poder calcular fi en función de Pu
# que disminuye a partir de pu=0 hasta Pu=0.10f'c
def fi_result(seccion, angulo, EjeN, defConc=vcc['defConc'], deftrac=vcc['deftrac'],
              ficomp=vcc['ficomp'], fitracc=vcc['fitracc'], alfa=vcc['alfa']):
    # defConc = 0.003 # Ɛ  deformación unitaria de falla del concreto
    # deftrac = 0.005 # Ɛ  deformación unitaria para falla por tracción
    # ficomp = 0.65   # ɸc para sección controlada por compresión
    # fitracc = 0.90  # ɸt para sección controlada por tracción
    # alfa = 1.0      # α  factor de amplificación por sobreresistencia del acero
    res = resultante(seccion, angulo, EjeN, defConc, alfa)
    # res[3]    # deformación de la varilla extrema
    # res[4]    # varilla extrema
    fi = ((fitracc-ficomp)/(deftrac-res[4].fy/res[4].Es)
          * (-res[3]-res[4].fy/res[4].Es)+ficomp)
    fi = min(fitracc, max(ficomp, fi))
    # print('fi resultante: angulo=%.4f c=%.4f, fi=%.2f' % (angulo, EjeN, fi))
    return res[0]*fi, res[1]*fi, res[2]*fi, fi, res[4]


il = pd.iteration_limit


def buscar_punto(sec, P, PxX, PxY, angulo=None, c=None):
    # sec.areas_varillas()
    # sec.fy_varillas()

    pos = PxY < 0

    if (angulo is None) and (c is None):
        asol = calc_ang(PxX, PxY, pos)
        c = buscar_c(sec, P, asol)
        res = fi_result(sec, asol, c)
        angulo = calc_ang(res[1], res[2], pos)  # primera aproximación del ángulo
        c = buscar_c(sec, P, angulo)
        # print('414 angulo=%.4f c=%.4f' % (angulo, c)')
    elif angulo is None:
        asol = calc_ang(PxX, PxY, pos)
        res = fi_result(sec, asol, c)
        angulo = calc_ang(res[1], res[2], pos)  # primera aproximación del ángulo
        c = buscar_c(sec, P, angulo)
        # print('420 angulo=%.4f c=%.4f' % (angulo, c)')
    else:
        c = buscar_c(sec, P, angulo)
        # print('423 angulo=%.4f c=%.4f' % (angulo, c)')

    res = fi_result(sec, angulo, c)
    errP = abs(P-res[0])
    errA = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
    if (errP < il['errP']*sec.Ag()*sec.fc) and (errA < il['errA']):
        # print('429 angulo=%.4f c=%.4f' % (angulo, c)')
        return [angulo, c]

    deltaC = 1
    deltaA = 1
    # print('434 angulo=%.4f c=%.4f' % (angulo, c)')
    ip = 0
    ipc = 0
    while (deltaC > il['errC']) or (deltaA > il['errA']):
        # print('437 angulo=%.4f c=%.4f' % (angulo, c)')
        # breakpoint()
        a1 = angulo
        c1 = c
        a2 = buscar_ang(sec, PxX, PxY, c1, a1)
        c2 = buscar_c(sec, P, a2, c1)
        a3 = buscar_ang(sec, PxX, PxY, c2, a2)
        c3 = buscar_c(sec, P, a3, c2)
        deltaA = abs(a3-a2)
        deltaC = abs(c3-c2)
        while ((abs(a3-a1) < il['errA']) and (abs(c3-c1) < il['errC'])
               and ((deltaC > il['errC']) or (deltaA > il['errA']))):
            # res = fi_result(sec, a1, c1)
            # errP11 = abs(P-res[0])
            # errA11 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # res = fi_result(sec, a2, c2)
            # errP22 = abs(P-res[0])
            # errA22 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # res = fi_result(sec, a3, c3)
            # errP33 = abs(P-res[0])
            # errA33 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # print('ipc=%d\n' % ipc
            #       + 'a1=%.4f c1=%.4f errP=%.4f errA=%.4f\n' % (a1, c1, errP11, errA11)
            #       + 'a2=%.4f c2=%.4f errP=%.4f errA=%.4f\n' % (a2, c2, errP22, errA22)
            #       + 'a3=%.4f c3=%.4f errP=%.4f errA=%.4f\n' % (a3, c3, errP33, errA33))
            a1 = (a3+a2)/2
            c1 = buscar_c(sec, P, a1, c3)
            a2 = buscar_ang(sec, PxX, PxY, c1, a1)
            c2 = buscar_c(sec, P, a2, c1)
            a3 = buscar_ang(sec, PxX, PxY, c2, a2)
            c3 = buscar_c(sec, P, a3, c2)
            ipc = ipc+1
            # res = fi_result(sec, a1, c1)
            # errP11 = abs(P-res[0])
            # errA11 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # res = fi_result(sec, a2, c2)
            # errP22 = abs(P-res[0])
            # errA22 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # res = fi_result(sec, a3, c3)
            # errP33 = abs(P-res[0])
            # errA33 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            # print('ipc=%d\n' % ipc
            #       + 'a1=%.4f c1=%.4f errP=%.4f errA=%.4f\n' % (a1, c1, errP11, errA11)
            #       + 'a2=%.4f c2=%.4f errP=%.4f errA=%.4f\n' % (a2, c2, errP22, errA22)
            #       + 'a3=%.4f c3=%.4f errP=%.4f errA=%.4f\n' % (a3, c3, errP33, errA33))
            while (((a2-a1)*(a1-a3) > 0) or ((c2-c1)*(c1-c3) > 0)
                   and ((deltaC > il['errC']) or (deltaA > il['errA']))):
                a1 = (a1+a2)/2
                c1 = buscar_c(sec, P, a1, c2)
                a2 = buscar_ang(sec, PxX, PxY, c1, a1)
                c2 = buscar_c(sec, P, a2, c1)
                a3 = buscar_ang(sec, PxX, PxY, c2, a2)
                c3 = buscar_c(sec, P, a3, c2)
                deltaA = abs(a3-a2)
                deltaC = abs(c3-c2)
                ipc = ipc+1
                # res = fi_result(sec, a1, c1)
                # errP11 = abs(P-res[0])
                # errA11 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
                # res = fi_result(sec, a2, c2)
                # errP22 = abs(P-res[0])
                # errA22 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
                # res = fi_result(sec, a3, c3)
                # errP33 = abs(P-res[0])
                # errA33 = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
                # print('--ipc=%d\n' % ipc
                #       + 'a1=%.4f c1=%.4f errP=%.4f errA=%.4f\n' % (a1, c1, errP11, errA11)
                #       + 'a2=%.4f c2=%.4f errP=%.4f errA=%.4f\n' % (a2, c2, errP22, errA22)
                #       + 'a3=%.4f c3=%.4f errP=%.4f errA=%.4f\n' % (a3, c3, errP33, errA33))
                if ipc > il['ip']:
                    res = fi_result(sec, a1, c1)
                    errP = abs(P-res[0])
                    errA = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
                    if (errP < il['errP']*sec.Ag()*sec.fc) and (errA < il['errA']):
                        # print('429 angulo=%.4f c=%.4f' % (angulo, c)')
                        return [a1, c1]
                    raise NameError('objetivo: P=%.1f angulo=%.4f\n' % (P, asol)
                                    + 'resultado: P=%.1f angulo=%.4f\n'
                                    % (res[0], calc_ang(res[1], res[2], pos))
                                    + 'Buscar ángulo y eje neutro no converge\n'
                                    + 'iteraciones ip=%d\n' % ip
                                    + 'iteraciones ipc=%d' % ipc)
            deltaA = abs(a3-a2)
            deltaC = abs(c3-c2)
        if (abs(a3-a1) < il['errA']) and (abs(c3-c1) < il['errC']):
            angulo = (a3+a2)/2
            c = buscar_c(sec, P, angulo, c3)
        else:
            angulo = a3
            c = c3
        # print('485 angulo=%.4f  c=%.4f\n' % (angulo, c)
        #       + '        a1=%.4f c1=%.4f\n' % (a1, c1)
        #       + '        a2=%.4f c2=%.4f\n' % (a2, c2)
        #       + '        a3=%.4f c3=%.4f\n' % (a3, c3))
        ip = ip+1
        # print('ip=', ip)
        if ip > il['ip']:
            res = fi_result(sec, angulo, c)
            errP = abs(P-res[0])
            errA = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
            if (errP < il['errP']*sec.Ag()*sec.fc) and (errA < il['errA']):
                # print('429 angulo=%.4f c=%.4f' % (angulo, c)')
                return [angulo, c]
            print('ipc=', ipc)
            raise NameError('objetivo: P=%.1f angulo=%.4f\n' % (P, asol)
                            + 'resultado: P=%.1f angulo=%.4f\n'
                            % (res[0], calc_ang(res[1], res[2], pos))
                            + 'Buscar ángulo y eje neutro iteraciones máximas ip>%d' % ip)
    # print('ip=', ip)
    # print('ipc=', ipc)
    # res = fi_result(sec, angulo, c)
    # errP = abs(P-res[0])
    # errA = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
    # print('objetivo: P=%.1f angulo=%.4f\n' % (P, asol)
    #       + 'resultado: P=%.1f angulo=%.4f\n'
    #       % (res[0], calc_ang(res[1], res[2], pos))
    #       + 'ai=%.4f ci=%.4f errP=%.4f errA=%.4f\n' % (angulo, c, errP, errA))
    return [angulo, c]


def buscar_c(sec, P, angulo, c=None):

    c1 = None
    if c is None:
        c1 = 0.4*min(sec.x, sec.y)
        c2 = c1+0.1*min(sec.x, sec.y)
    else:
        c1 = max(0, c-0.05*min(sec.x, sec.y))
        c2 = c
    p1 = fi_result(sec, angulo, c1)[0]
    p2 = fi_result(sec, angulo, c2)[0]
    # breakpoint()
    # print('477 P=%.1f p1=%.1f p2=%.1f' % (P, p1, p2))
    # print('478 ci=%.4f c1=%.4f c2=%.4f' % ((c2-c1)/(p2-p1)*(P-p1)+c1, c1, c2))
    ic = 0
    while abs(P-p2) > (il['errP']*sec.Ag()*sec.fc):
        ci = (c2-c1)/(p2-p1)*(P-p1)+c1
        if ci < 0:
            ci = 0
        # print('484 P=%.1f p1=%.1f p2=%.1f' % (P, p1, p2))
        # print('485 ci=%.4f c1=%.4f c2=%.4f' % (ci, c1, c2))
        # breakpoint()
        c1 = c2
        p1 = p2
        c2 = ci
        p2 = fi_result(sec, angulo, c2)[0]
        ic = ic+1
        if ic > il['ic']:
            raise NameError('P=%.1f p1=%.1f p2=%.1f' % (P, p1, p2)
                            + 'ci=%.4f c1=%.4f c2=%.4f' % ((c2-c1)/(p2-p1)*(P-p1)+c1, c1, c2)
                            + 'Buscar eje neutro no converge ic>%d' % ic)
    # print('c2', angulo, c2)
    # print('ic=', ic)
    return c2


def buscar_ang(sec, PxX, PxY, c, a=None):
    if c == 0:
        raise NameError("c=0, ángulo inderterminado")
    pos = PxY < 0
    asol = calc_ang(PxX, PxY, pos)
    # breakpoint()
    if a is None:
        res = resultante(sec, asol, c)
        a = calc_ang(res[1], res[2], False)  # primera aproximación
    else:
        a = 45
    fasol = 90-asol
    acs = fasol-fasol % 90  # la pendiente es negativa entonces acs es menor aci
    aci = acs+90
    if a < acs or aci < a:
        a = acs+a % 90
    res = resultante(sec, a, c)
    ar = calc_ang(res[1], res[2], pos)
    if abs(asol-ar) < il['errA']/10:
        return a
    if abs(asol-ar) > 90:
        raise NameError("caudránte equivocado")
    ari = ar-ar % 90
    ars = ari+90
    if asol < ari or ars < asol:
        raise NameError("caudránte equivocado")
    if asol > ar:
        aci = a
        ari = ar
        ars = ar-ar % 90+90
    else:
        acs = a
        ars = ar
        ari = ar-ar % 90
    a = a-copysign(0.1, asol-ar)
    res = resultante(sec, a, c)
    ar = calc_ang(res[1], res[2], pos)
    ia = 0
    while abs(asol-ar) > il['errA']/10:
        # print('640 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
        # print('641 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
        # if ars-ari == 0:
        #     breakpoint()
        # print('pendiente=', (acs-aci)/(ars-ari))
        if (ars+ari)/2-asol > 0:  # el siguiente punto se calcula con el más cercano a la solición
            ar2 = ari
            a2 = aci
        else:
            ar2 = ars
            a2 = acs
        if asol > ar:
            aci = a
            ari = ar
        else:
            acs = a
            ars = ar
        if asol < ari or ars < asol:
            raise NameError("Valor fuera de los límites")
        a = a+(asol-ar)*(a-a2)/(ar-ar2)
        if a < acs or aci < a:
            a = (aci+acs)/2
        res = resultante(sec, a, c)
        ar = calc_ang(res[1], res[2], pos)

        ia = ia+1
        if ia > il['ia']:
            raise NameError('660 asol=%.4f ar=%.4f ari=%.4f ars=%.4f\n'
                            % (asol, ar, ari, ars)
                            + '662 asol=%.4f a=%.4f aci=%.4f acs=%.4f\n'
                            % (asol, a, aci, acs)
                            + 'Buscar ángulo no converge ia=%d' % (ia))
    # print('ia=%d asol=%.4f ar=%.4f a=%.4f c=%.4f' % (ia, asol, ar, a, c))
    return a
