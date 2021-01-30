import Seccion as s
import public_data as pd


vcc = pd.column_design_factors


class solicitation():

    def __init__(self, P, Vx, Vy, T, Mx, My):
        self.P = P    # Fz Axial
        self.Vx = Vx  # Fx Cortante, Shear
        self.Vy = Vy  # Fy Cortante, Shear
        self.T = T    # Mz-z Momento torsor, torsion
        self.Mx = Mx  # Fz*x=My-y Momento flector, flexion
        self.My = My  # Fz*y=Mx-x Momento flector, flexion

    def pmm(self):
        return self.P, self.Mx, self.My

    def ang_sol(self):
        return s.calc_ang(self.Mx, self.My)

    def str(self):
        return "P=%.1f, Vx=%.1f, Vy=%.1f, T=%.1f, Mx=%.1f, My=%.1f"\
            % (self.P, self.Vx, self.Vy, self.T, self.Mx, self.My)


# Ínide de sobreessfuerzo a flexocompresión
def indice_flco(sec, sol, fic=vcc['ficomp'], fit=vcc['fitracc'],
                ec=vcc['defConc'], et=vcc['deftrac'], alfa=vcc['alfa']):
    a, c = s.buscar_punto(sec, sol.P, sol.Mx, sol.My)
    res = s.fi_result(sec, a, c, defConc=ec, deftrac=et,
                      ficomp=fic, fitracc=fit, alfa=alfa)
    return (sol.Mx**2 + sol.My**2)**0.5 / (res[1]**2 + res[2]**2)**0.5
