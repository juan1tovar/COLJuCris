import numpy as np
# import math
from math import radians, sin, cos, copysign, atan, tan, degrees
from operator import itemgetter
np.set_printoptions(precision=3, suppress=True)


def np_err_handler(type, flag):
    print("Floating point error (%s), with flag %s" % (type, flag))
    print('form more info search for numpy.seterrcall')
    breakpoint()


np.seterrcall(np_err_handler)
np.seterr(all='call')


def rotar(coordenadas, grados, invertir=False):
    # la operación de módulo me aproxima los números negativos
    # por eso me tocó meter todo en esa fórmula
    a = radians(-grados)
    t = np.array([[cos(a), -sin(a)],
                  [sin(a), cos(a)]])
    if invertir:
        return t.dot(coordenadas)
    else:
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


class seccion:

    def __init__(self, fc, LadoX, LadoY):
        self.x = LadoX  # metros
        self.y = LadoY  # metroskb
        self.psi = atan(self.y/self.x)
        self.fc = fc    # MPa
        self.fy = []    # MPa
        self.Es = 200000    # MPa
        self.defConc = 0.003    # Ɛ deformación unitaria
        self.deftrac = 0.005    # Ɛ deformación unitaria falla por tracción
        self.ficomp = 0.65       # ɸc para sección controlada por compresión
        self.fitracc = 0.90     # ɸt para sección controlada por tracción
        self.refActual = ''
        self.AreasRef = []
        self.refXY = []
        self.varillas = []
        self.esquinas = np.array([[0+self.x/2, 0-self.y/2],
                                  [0+self.x/2, 0+self.y/2],
                                  [0-self.x/2, 0+self.y/2],
                                  [0-self.x/2, 0-self.y/2]])

# -------------------------- Métodos de concreto -------------------------
    def Ec_MPa(self):
        return 4700*self.fc**0.5

    def Ag(self):
        return self.x*self.y

    def cord_conc(self, angulo):
        return rotar(self.esquinas, angulo)


# -------------------------- Métodos de refuerzo -------------------------
    tablaRefuerzo = {'Tag': ['D[mm]', 'A[mm²]'],
                     '#3': [9.5, 71],
                     '#4': [12.7, 129],
                     '#5': [15.9, 199],
                     '#6': [19.1, 284],
                     '#7': [22.2, 387],
                     '#8': [25.4, 510],
                     '#9': [28.7, 645],
                     '#10': [32.3, 819]}

    def intercalado(self, NVarX, NVarY, DEsqui, DAlter, DFlej='#3', recub=0.04):
        DE = self.tablaRefuerzo[DEsqui][0]/1000
        DA = self.tablaRefuerzo[DAlter][0]/1000
        DF = self.tablaRefuerzo[DFlej][0]/1000
        Xmax = self.x/2-recub-DF
        Ymax = self.y/2-recub-DF
        sepx = (Xmax*2-DE)/(NVarX-1)
        if sepx < 0:
            raise "Intercalado: sepx<0"
        sepy = (Ymax*2-DE)/(NVarY-1)
        if sepy < 0:
            raise "Intercalado: sepy<0"
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
            varx[0][n] = [x, y, DEsqui]
            varx[1][n] = [-x, y, DEsqui]
            varx[2][n] = [-x, -y, DEsqui]
            varx[3][n] = [x, -y, DEsqui]
            if n+1 == int(NVarX/2)+NVarX % 2:
                break
            x = Xmax-DE/2-sepx*(n+1)
            y = Ymax-DA/2
            varx[0][n+1] = [x, y, DAlter]
            varx[1][n+1] = [-x, y, DAlter]
            varx[2][n+1] = [-x, -y, DAlter]
            varx[3][n+1] = [x, -y, DAlter]

        for n in range(0, int(NVarY/2)+NVarY % 2, 2):
            x = Xmax-DE/2
            y = Ymax-DE/2-sepy*n
            vary[0][n] = [x, y, DEsqui]
            vary[1][n] = [-x, y, DEsqui]
            vary[2][n] = [-x, -y, DEsqui]
            vary[3][n] = [x, -y, DEsqui]
            if n+1 == int(NVarY/2)+NVarY % 2:
                break
            x = Xmax-DA/2
            y = Ymax-DE/2-sepy*(n+1)
            vary[0][n+1] = [x, y, DAlter]
            vary[1][n+1] = [-x, y, DAlter]
            vary[2][n+1] = [-x, -y, DAlter]
            vary[3][n+1] = [x, -y, DAlter]

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
        self.areas_varillas()
        self.fy_varillas()

    def esquinero(self, NVarX, NVarY, DEsqui, DAlter, DFlej='#3', recub=0.04):
        DE = self.tablaRefuerzo[DEsqui][0]/1000
        DA = self.tablaRefuerzo[DAlter][0]/1000
        DF = self.tablaRefuerzo[DFlej][0]/1000
        Xmax = self.x/2-recub-DF
        Ymax = self.y/2-recub-DF
        sepx = (Xmax*2-DE)/(NVarX-1)
        if sepx < 0:
            raise "Esquinero: sepx<0"
        sepy = (Ymax*2-DE)/(NVarY-1)
        if sepy < 0:
            raise "Esquinero: sepy<0"
        ReSi = []
        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([x, -y, DEsqui])
        x = Xmax-DA/2

        for i in range(1, NVarY-1):
            y = Ymax-DE/2-sepy*i
            ReSi.append([x, -y, DAlter])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([x, y, DEsqui])
        y = Ymax-DA/2

        for i in range(1, NVarX-1):
            x = Xmax-DE/2-sepx*i
            ReSi.append([x, y, DAlter])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([-x, y, DEsqui])
        x = Xmax-DA/2

        for i in range(1, NVarY-1):
            y = Ymax-DE/2-sepy*i
            ReSi.append([-x, y, DAlter])

        x = Xmax-DE/2
        y = Ymax-DE/2
        ReSi.append([-x, -y, DEsqui])
        y = Ymax-DA/2

        for i in range(1, NVarX-1):
            x = Xmax-DE/2-sepx*i
            ReSi.append([-x, -y, DAlter])

        self.refXY = np.array([var[0:2] for var in ReSi])
        self.varillas = [var[2] for var in ReSi]
        self.areas_varillas()
        self.fy_varillas()

    def cord_ref(self, angulo):

        if self.refXY == []:
            raise NameError('No se han ingresado varillas de refuerzo')

        return rotar(self.refXY, angulo)

    def areas_varillas(self):

        if self.refXY == []:
            raise NameError('No se han ingresado varillas de refuerzo')

        self.AreasRef = []
        for var in self.varillas:
            self.AreasRef.append(self.tablaRefuerzo[var][1])
        # return areas

    def fy_varillas(self):
        self.fy = 420*np.ones(len(self.AreasRef))

# ---------------------------- Métodos de la sección ----------------------------
    def set_ref_sim(self, NVarX, NVarY, DEsqui, DAlter,
                    opcion='intercalado', DFlej='#3', recub=0.04):
        func = getattr(self, opcion, 'OpcInvalida')

        if func == 'OpcInvalida':
            raise NameError('No existe esa opción de refuerzo simétrico')

        func(NVarX, NVarY, DEsqui, DAlter, DFlej, recub)
        self.refActual = opcion

    def Pnmax(self):
        As = sum(self.AreasRef)
        Ag = self.Ag()
        AsFy = 0
        for i, a in enumerate(self.AreasRef):
            AsFy = AsFy+self.fy[i]*a/1000
        return 0.75*(0.85*(Ag-As/10**6)*self.fc*1000+AsFy)
        # 0.75 es un factor que depende del código o del material

    def result_conc(self, angulo, ejeN):
        cor = sorted(self.cord_conc(angulo), key=itemgetter(0))
        cor = sorted(cor, key=itemgetter(1))
        ymax = cor[3][1]
        yEje = ymax-ejeN
        betac = ejeN*max((min((0.85, 0.85-0.05*(self.fc-28)/7)), 0.65))

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
                hpar = min((betac, self.y))
                bpar = self.x
            else:
                ta = 1000000
                hpar = min((betac, self.x))
                bpar = self.y
        else:
            htri = min((ht, betac))
            hpar = min((cor[2][1]-cor[1][1], betac-htri))
            htra = min((ht, betac-hpar-htri))
            angbase = abs(angulo) % 180.0  # Este angulo se usa para calcular el bloque de Witney
            if angbase > 90:
                angbase = 180-angbase
            if radians(angbase) < self.psi:
                bpar = self.x/cos(radians(angbase))
                ta = tan(radians(angbase))
            else:
                bpar = self.y/cos(radians(90-angbase))
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
             x1-copysign(1, x1)*(btra**2/2+htra*(btra/ta+bpar/2*ta) +
                                 htra**2/3*(1/ta**2-ta**2))/((bpar+btra)/2)]
        Y = [ymax-2*htri/3,
             cor[2][1]-hpar/2,
             cor[1][1]-htra/3*(bpar+2*btra)/(bpar+btra)]
        Ag = sum(A)
        # breakpoint()
        A = np.array(A)
        Fuerza = Ag*0.85*self.fc*1000
        xyro = np.array([[sum(A*np.array(X))/Ag],
                        [sum(A*np.array(Y))/Ag]])
        xy = rotar(xyro, angulo, invertir=True)
        # breakpoint()
        return [Fuerza, Fuerza*xy[0, 0], Fuerza*xy[1, 0], yEje]

    def result_ref(self, angulo, ejeN, yEje):
        cor = self.cord_ref(angulo)
        defmin = 1
        fuerza = []
        for i, v in enumerate(cor):
            du = (v[1]-yEje)*self.defConc/ejeN
            if du < defmin:
                defmin = du
            if du > 0:
                fuerza.append(min(self.fy[i], self.Es*du)*self.AreasRef[i]/1000)
            else:
                fuerza.append(max(-self.fy[i], self.Es*du)*self.AreasRef[i]/1000)
        fuerza = np.array(fuerza)
        Fuerza = np.sum(fuerza)
        FxX = fuerza.dot(self.refXY[:, 0])
        FxY = fuerza.dot(self.refXY[:, 1])
        # print('Refuerzo', Fuerza, FxX, FxY)
        return [Fuerza, FxX, FxY, defmin]

    def resultante(self, angulo, ejeN):
        if self.refXY == []:
            raise NameError('No se han ingresado varillas de refuerzo')
        if self.fy == []:
            raise NameError('No se ha asignado fy de cada varilla')
        if self.AreasRef == []:
            raise NameError('No se ha determinado el area de varillas de refuerzo')

        if ejeN <= 0:
            AsFy = 0
            for i, a in enumerate(self.AreasRef):
                AsFy = AsFy+self.fy[i]*a/1000
            return [-AsFy, 0, 0, -self.deftrac]
        # xyConc = self.cord_conc(angulo)
        # if ejeN >= max(xyConc[:, 1])*2*self.defConc/(self.defConc-max(self.fy)/self.Es):
        #     return [self.Pnmax(), 0, 0, max(self.fy)/self.Es]
        # if abs(angulo) < 0.0001:
            # breakpoint()

        rConc = self.result_conc(angulo, ejeN)
        rRefu = self.result_ref(angulo, ejeN, rConc[3])
        F = rConc[0]+rRefu[0]
        FxX = rConc[1]+rRefu[1]
        FxY = rConc[2]+rRefu[2]
        # breakpoint()
        return [F, FxX, FxY, rRefu[3]]

# Este método se debe ampliar para poder calcular con un fi general y
#  fi en función de Pu que disminuye a partir de pu=0 hasta Pu=0.10f'c
    def fi_result(self, angulo, EjeN):
        res = self.resultante(angulo, EjeN)
        fi = ((self.fitracc-self.ficomp)/(self.deftrac-max(self.fy)/self.Es)
              * (-res[3]-max(self.fy)/self.Es)+self.ficomp)
        fi = min((self.fitracc, max((self.ficomp, fi))))
        # print('fi resultante: angulo, c, fi', angulo, EjeN, fi)
        return res[0]*fi, res[1]*fi, res[2]*fi, fi

    def buscar_punto(self, P, PxX, PxY, angulo=None, c=None):
        # self.areas_varillas()
        # self.fy_varillas()

        pos = PxY < 0

        if (angulo is None) and (c is None):
            asol = calc_ang(PxX, PxY, pos)
            c = self.buscar_c(P, asol)
            res = self.resultante(asol, c)
            angulo = calc_ang(res[1], res[2], pos)  # primera aproximación del ángulo
            c = self.buscar_c(P, angulo)
            # print('1', angulo, c)
        elif angulo is None:
            asol = calc_ang(PxX, PxY, pos)
            res = self.resultante(asol, c)
            angulo = calc_ang(res[1], res[2], pos)  # primera aproximación del ángulo
            c = self.buscar_c(P, angulo)
            # print('2', angulo, c)
        else:
            c = self.buscar_c(P, angulo)
            # print('3', angulo, c)

        res = self.fi_result(angulo, c)
        errP = abs(P-res[0])
        errA = abs(calc_ang(PxX, PxY, pos)-calc_ang(res[1], res[2], pos))
        if (errP < 0.1) and (errA < 0.01):
            # print('4', angulo, c)
            return [angulo, c]

        deltaC = 1
        deltaA = 1
        # print('5', angulo, c)
        ip = 0
        while (deltaC > 0.0001) or (deltaA > 0.01):
            # print('6', angulo, c)
            # breakpoint()
            a1 = angulo
            c1 = c
            a2 = self.buscar_ang(PxX, PxY, c1, a1)
            c2 = self.buscar_c(P, a2, c1)
            a3 = self.buscar_ang(PxX, PxY, c2, a2)
            c3 = self.buscar_c(P, a3, c2)
            deltaA = abs(a3-a2)
            deltaC = abs(c3-c2)
            if (abs(a3-a1) < 0.01) and (abs(c3-c1) < 0.0001):
                angulo = (a3+a2)/2
                c = (c3+c2)/2
            else:
                angulo = a3
                c = c3
            # print('7', angulo, c)
            ip = ip+1
            print('ip=', ip)
            if ip > 10:
                raise 'Buscar ángulo y eje neutro no converge'
        return [angulo, c]

    def buscar_c(self, P, angulo, c=None):

        c1 = None
        if c is None:
            c1 = 0.4*min((self.x, self.y))
            c2 = c1+0.1*min((self.x, self.y))
        else:
            c1 = c-0.05*min((self.x, self.y))
            c2 = c
        p1 = self.fi_result(angulo, c1)[0]
        p2 = self.fi_result(angulo, c2)[0]
        # breakpoint()
        # print('380 P=%.1f p1=%.1f p2=%.1f' % (P, p1, p2))
        # print('381 ci=%.4f c1=%.4f c2=%.4f' % ((c2-c1)/(p2-p1)*(P-p1)+c1, c1, c2))
        ic = 0
        while abs(P-p2) > (0.0001*self.Ag()*self.fc):
            ci = (c2-c1)/(p2-p1)*(P-p1)+c1
            # print('422 P p1 p2', P, p1, p2)
            # print('423 ci c1 c2', ci, c1, c2)
            # breakpoint()
            c1 = c2
            p1 = p2
            c2 = ci
            p2 = self.fi_result(angulo, c2)[0]
            ic = ic+1
            if ic > 50:
                raise 'Buscar eje neutro no converge'
        # print('c2', angulo, c2)
        print('ic=', ic)
        return c2

    def buscar_ang(self, PxX, PxY, c, a=None):
        pos = PxY < 0
        asol = calc_ang(PxX, PxY, pos)
        # breakpoint()
        if a is None:
            res = self.resultante(asol, c)
            a = calc_ang(res[1], res[2], False)  # primera aproximación
        res = self.resultante(a, c)
        ar = calc_ang(res[1], res[2], pos)
        if abs(asol-ar) < 0.001:
            return a
        if a < -180:
            a = max(-180, a)
        if a > 450:
            a = min(450, a)
        aci = a
        ari = ar
        acs = a-copysign(0.1, asol-ar)
        res = self.resultante(acs, c)
        ars = calc_ang(res[1], res[2], pos)  # primera aproximación
        if acs < -180:
            acs = max(-180, acs)
        if acs > 450:
            acs = min(450, acs)
        if ars < ari:
            a = acs
            ar = ars
            acs = aci
            ars = ari
            aci = a
            ari = ar
        a = aci+(asol-ari)*(acs-aci)/(ars-ari)
        res = self.resultante(a, c)
        ar = calc_ang(res[1], res[2], pos)
        if a < -180:
            a = max(-180, a)
        if a > 450:
            a = min(450, a)
        ia = 1
        while (asol < ari or ars < asol) and abs(asol-ar) > 0.001:
            # print('508 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
            # print('509 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
            if ar < asol:
                if ari < asol:
                    if ari < ar:
                        aci = a
                        ari = ar
                    else:
                        print('516 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
                        print('517 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
                        raise 'Buscar ángulo no converge'
                else:
                    aci = a
                    ari = ar
            else:
                if asol < ars:
                    if ar < ars:
                        acs = a
                        ars = ar
                    else:
                        print('528 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
                        print('529 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
                        raise 'Buscar ángulo no converge'
                else:
                    acs = a
                    ars = ar
            if ars < ari:
                a = acs
                ar = ars
                acs = aci
                ars = ari
                aci = a
                ari = ar
            a = aci+(asol-ari)*(acs-aci)/(ars-ari)
            # if a > 360 or a < 0:
            #     a = ajustar_angulo(a, not pos)
            res = self.resultante(a, c)
            ar = calc_ang(res[1], res[2], pos)
            ia = ia+1
            if ia > 100:
                raise 'Buscar ángulo no converge'
        while abs(asol-ar) > 0.001:
            # print('550 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
            # print('551 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
            # if ars-ari == 0:
            #     breakpoint()
            # print('pendiente=', (acs-aci)/(ars-ari))
            if ari < ar and ar < asol:
                aci = a
                ari = ar
            elif asol < ar and ar < ars:
                acs = a
                ars = ar
            else:
                print('556 asol=%.4f ar=%.4f ari=%.4f ars=%.4f' % (asol, ar, ari, ars))
                print('563 asol=%.4f a=%.4f aci=%.4f acs=%.4f' % (asol, a, aci, acs))
                raise 'Buscar ángulo no converge'
            a = aci+(asol-ari)*(acs-aci)/(ars-ari)
            if a < -180:
                a = max(-180, a)
            if a > 450:
                a = min(450, a)
            res = self.resultante(a, c)
            ar = calc_ang(res[1], res[2], pos)
            ia = ia+1
            if ia > 100:
                raise 'Buscar ángulo no converge ia>50'
        print('ia=%d asol=%.4f ar=%.4f a=%.4f c=%.4f' % (ia, asol, ar, a, c))
        return a
