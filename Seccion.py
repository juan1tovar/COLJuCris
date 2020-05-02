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
    a = radians(-(grados % 360))
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

    if positive and (ang < 0):
        ang = ang + 360
    return ang


def ajustar_angulo(ang, positive=True):
    ang = copysign(abs(ang) % 360, ang)
    if (positive and (ang < 0)) or (not positive and (ang < 90)):
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
        sepy = (Ymax*2-DE)/(NVarY-1)
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
        sepy = (Ymax*2-DE)/(NVarY-1)
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
        cor = sorted(self.cord_conc(angulo), key=itemgetter(1))
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
            ht = 0.00000000001
            # htri, btri, htra, btra = 0, 0, 0, 0, 0
            if (angulo == 0) or (angulo == 180):
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
             x1-copysign(1, x1)*(htra/ta+btra/2+((ta-1/ta)*btra*htra/2
                                 + htra**2/3*(ta**2-1/ta**2))/(bpar+btra)/2)]
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

        rConc = self.result_conc(angulo, ejeN)
        rRefu = self.result_ref(angulo, ejeN, rConc[3])
        F = rConc[0]+rRefu[0]
        FxX = rConc[1]+rRefu[1]
        FxY = rConc[2]+rRefu[2]
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

        pos = ((PxY/PxX < -1) or (PxY/PxX > 1)) and (PxY < 0)

        if (angulo is None) and (c is None):
            asol = calc_ang(PxX, PxY, pos)
            c = self.buscar_c(P, asol)
            res = self.resultante(asol, c)
            angulo = calc_ang(res[1], res[2], positive=pos)  # primera aproximación del ángulo
            c = self.buscar_c(P, angulo)
            # print('1', angulo, c)
        elif angulo is None:
            asol = calc_ang(PxX, PxY, positive=pos)
            res = self.resultante(asol, c)
            angulo = calc_ang(res[1], res[2], positive=pos)  # primera aproximación del ángulo
            # print('2', angulo, c)
        else:
            c = self.buscar_c(P, angulo)
            # print('3', angulo, c)

        res = self.fi_result(angulo, c)
        errP = abs(P-res[0])
        errA = abs(calc_ang(PxX, PxY, positive=pos)
                   - calc_ang(res[1], res[2], positive=pos))
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
            a2 = self.buscar_ang(PxX, PxY, c1, ac1=a1)
            c2 = self.buscar_c(P, a2, c=c1)
            a3 = self.buscar_ang(PxX, PxY, c2, ac1=a2)
            c3 = self.buscar_c(P, a3, c=c2)
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
            if ip > 50:
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
        # print('380 P p1 p2', P, p1, p2)
        # print('381 ci c1 c2', (c2-c1)/(p2-p1)*(P-p1)+c1, c1, c2)
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
        pos = ((PxY/PxX < -1) or (PxY/PxX > 1)) and (PxY < 0)
        asol = calc_ang(PxX, PxY, pos)

        aci = asol
        res = self.resultante(aci, c)
        ari = calc_ang(res[1], res[2], not pos)  # primera aproximación
        if a is None:
            a = ari
            acs = ari
        else:
            acs = a
        res = self.resultante(acs, c)
        ars = calc_ang(res[1], res[2], pos)
        ar = ars
        if ars < ari:
            acs = aci
            ars = ari
            aci = a
            ari = ar
            a = aci-copysign(0.1, asol-ari)
        else:
            a = acs-copysign(0.1, asol-ars)
        res = self.resultante(a, c)
        ar = calc_ang(res[1], res[2], pos)
        ia = 0
        while (asol < aci or acs < asol) and abs(asol-ar) > 0.001:
            if a < asol:
                if aci < asol:
                    if aci < a:
                        aci = a
                        ari = ar
                    else:
                        raise 'Buscar ángulo no converge'
                else:
                    aci = a
                    ari = ar
            else:
                if asol < acs:
                    if a < acs:
                        acs = a
                        ars = ar
                    else:
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
            if a > 360 or a < 0:
                a = ajustar_angulo(a, not pos)
            res = self.resultante(a, c)
            ar = calc_ang(res[1], res[2], pos)
            ia = ia+1
            if ia > 50:
                raise 'Buscar ángulo no converge'
        while abs(asol-ar1) > 0.001:
            print('419 asol=%.4f ar1=%.4f ar2=%.4f' % (asol, ar1, ar2))
            print('420 asol=%.4f ac1=%.4f ac2=%.4f' % (asol, ac1, ac2))
            # if ar2-ar1 == 0:
            #     breakpoint()
            # print('pendiente=', (ac2-ac1)/(ar2-ar1))
            acn = ac2
            arn = ar2
            ac2 = ac1+(asol-ar1)*(ac2-ac1)/(ar2-ar1)
            if ac2 > 360 or ac2 < 0:
                ac2 = ajustar_angulo(ac2, not pos)
            res = self.resultante(ac2, c)
            ar2 = calc_ang(res[1], res[2], pos)
            # div = 1.0
            # while abs(asol-ar1) < abs(asol-ar2):
            #     print('412 asol=%.4f ar1=%.4f ar2=%.4f div=%d' % (asol, ar1, ar2, div))
            #     print('413 asol=%.4f ac1=%.4f ac2=%.4f div=%d' % (asol, ac1, ac2, div))
            #     # breakpoint()
            #     div = div*2
            #     ac2 = ac1+(asol-ar1)*(ac2-ac1)/(ar2-ar1)/(div+.1)
            #     if abs(ac2) > 360:
            #         ac2 = ajustar_angulo(ac2, not pos)
            #     res = self.resultante(ac2, c)
            #     ar2 = calc_ang(res[1], res[2], pos)
            #     print('div =', div)
            ac1 = acn
            ar1 = arn
            # print('erra3', ac1, ar1)
            ia = ia+1
            if ia > 50:
                raise 'Buscar ángulo no converge'
        print('ia=%d c=%.4f' % (ia, c))
        # print('a2', ac1, ar1)
        return ac1