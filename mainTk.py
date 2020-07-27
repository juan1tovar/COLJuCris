import sys
import debugexc
import numpy as np
from Seccion import seccion
# from Seccion import calc_ang
import curvas
import graficas
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler

sys.excepthook = debugexc.info
np.set_printoptions(precision=3, suppress=True)


class chart_frame (tk.Frame, graficas.grafica):

    def __init__(self, Frame, is_3D):
        tk.Frame.__init__(self, Frame)
        graficas.grafica.__init__(self, is_3D)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.mpl_connect("key_press_event", self.on_key_press)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

    def on_key_press(self, event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, self.canvas, self.toolbar)

    def draw(self):
        self.canvas.draw()
        self.toolbar.update()


root = tk.Tk()
root.title("ColJuCris")
root.geometry("800x800")
f_datos = tk.Frame(root)
f_graficas = tk.Frame(root)
f_3D = chart_frame(f_graficas, is_3D=True)
f_hor = chart_frame(f_graficas, is_3D=False)
f_ver = chart_frame(f_graficas, is_3D=False)
# f_3D = tk.Frame(f_graficas)
# f_hor = tk.Frame(f_graficas)
# f_ver = tk.Frame(f_graficas)
# graf = graficas.grafica()
# canvas_3D = FigureCanvasTkAgg(graf.fig3D, f_3D)
# canvas_3D.draw()
# canvas_3D.mpl_connect("key_press_event", on_key_press_3D)
# toolbar_3D = NavigationToolbar2Tk(canvas_3D, f_3D)
# toolbar_3D.update()
# canvas_hor = FigureCanvasTkAgg(graf.fighor, f_hor)
# canvas_hor.draw()
# canvas_hor.mpl_connect("key_press_event", on_key_press_hor)
# toolbar_hor = NavigationToolbar2Tk(canvas_hor, f_hor)
# toolbar_hor.update()
# canvas_ver = FigureCanvasTkAgg(graf.figver, f_ver)
# canvas_ver.draw()
# canvas_ver.mpl_connect("key_press_event", on_key_press_ver)
# toolbar_ver = NavigationToolbar2Tk(canvas_ver, f_ver)
# toolbar_ver.update()


# FUNCTIONS---------------------------------------------------------------------------
def get_col(Fc, X, Y, vX, vY, d1, d2, df, opc, rec):
    sec = seccion(Fc, X, Y)
    sec.set_ref_sim(vX, vY, d1, d2, opc, df, rec)
    return sec


def bf_graficar():
    sec = get_col(
        cv_fc.get(),
        cv_ladoX.get(),
        cv_ladoY.get(),
        cv_varX.get(),
        cv_varY.get(),
        cv_dia1.get(),
        cv_dia2.get(),
        cv_diaf.get(),
        cv_conf_ref.get(),
        cv_recu.get()
    )
    if cv_bool3D.get() or cv_boolver.get():
        pmm = np.array(curvas.vertical(sec, cv_angulo.get()))
        if cv_bool3D.get():
            f_3D.addPMM(pmm)
            f_3D.draw()
        if cv_boolver.get():
            f_ver.addPM(pmm, cv_angulo.get())
            f_ver.draw()
    if cv_boolhor.get():
        if cv_bool_fatcircle.get():
            parlem = np.array(curvas.horizontal(sec, cv_Pu.get(), metodo='fatcircle'))
            f_hor.add2D(parlem, "parlem")
        if cv_bool_ang_sol.get():
            solici = np.array(curvas.horizontal(sec, cv_Pu.get(), metodo='ang_sol'))
            f_hor.add2D(solici, "ang sol")
        if cv_bool_column.get():
            column = np.array(curvas.horizontal(sec, cv_Pu.get(), metodo='ang_col'))
            f_hor.add2D(column, "ang_col")
        f_hor.draw()
    print("axes añadidos")
    # graf.update()


def bf_clear():
    # graf.clear()
    f_3D.cla()
    f_hor.cla()
    f_ver.cla()
    f_3D.draw()
    f_hor.draw()
    f_ver.draw()


def bf_salir():
    quit()


# VARIABLES DE CONTROL (cv)------------------------------------------------------------
cv_fc = tk.DoubleVar(value=28)
cv_ladoX = tk.DoubleVar(value=0.3)
cv_ladoY = tk.DoubleVar(value=0.3)
cv_angulo = tk.DoubleVar(value=0)
cv_c = tk.DoubleVar(value=.01)
cv_Pu = tk.DoubleVar(value=100)
cv_varX = tk.IntVar(value=4)
cv_varY = tk.IntVar(value=4)
cv_dia1 = tk.StringVar(value='#5')
cv_dia2 = tk.StringVar(value='#5')
cv_diaf = tk.StringVar(value='#3')
cv_conf_ref = tk.StringVar(value='esquinero')
cv_recu = tk.DoubleVar(value=.04)

cv_bool3D = tk.BooleanVar(value=True)
cv_boolhor = tk.BooleanVar(value=True)
cv_boolver = tk.BooleanVar(value=True)
cv_bool_fatcircle = tk.BooleanVar(value=True)
cv_bool_ang_sol = tk.BooleanVar(value=True)
cv_bool_column = tk.BooleanVar(value=True)
# LABELS------------------------------------------------------------------------------
et_fc = tk.Label(f_datos, text="f'c [MPa]=")
et_ladoX = tk.Label(f_datos, text="Lado X [m]=")
et_ladoY = tk.Label(f_datos, text="Lado Y [m]=")
et_angulo = tk.Label(f_datos, text="Ángulo [°]=")
et_c = tk.Label(f_datos, text="c [m]=")
et_Pu = tk.Label(f_datos, text="Pu [KN]=")
et_varX = tk.Label(f_datos, text="Varillas en X=")
et_varY = tk.Label(f_datos, text="Varillas en Y=")
et_dia1 = tk.Label(f_datos, text="Ø_1 [#/8'']=")
et_dia2 = tk.Label(f_datos, text="Ø_2 [#/8'']=")
et_diaf = tk.Label(f_datos, text="Ø_fleje [#/8'']=")
et_conf_ref = tk.Label(f_datos, text="Op. Refuerzo=")
et_recu = tk.Label(f_datos, text="Recubrimiento [cm]=")
# TEXT BOXES--------------------------------------------------------------------------
tb_fc = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_fc)
tb_ladoX = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_ladoX)
tb_ladoY = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_ladoY)
tb_angulo = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_angulo)
tb_c = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_c)
tb_Pu = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_Pu)
tb_varX = tk.Entry(f_datos, width=2, borderwidth=3, textvariable=cv_varX)
tb_varY = tk.Entry(f_datos, width=2, borderwidth=3, textvariable=cv_varY)
tb_dia1 = tk.Entry(f_datos, width=3, borderwidth=3, textvariable=cv_dia1)
tb_dia2 = tk.Entry(f_datos, width=3, borderwidth=3, textvariable=cv_dia2)
tb_diaf = tk.Entry(f_datos, width=3, borderwidth=3, textvariable=cv_diaf)
tb_conf_ref = tk.Entry(f_datos, width=10, borderwidth=3, textvariable=cv_conf_ref)
tb_recu = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_recu)
# BUTTONS-----------------------------------------------------------------------------
btt_graficar = tk.Button(root, text="Graficar", command=bf_graficar)
btt_clear = tk.Button(root, text="Limpiar", command=bf_clear)
btt_salir = tk.Button(root, text="salir", command=bf_salir)
# UBICACIÓN DE LOS WIDGETS------------------------------------------------------------
et_fc.grid(row=0, column=0, sticky='E')
tb_fc.grid(row=0, column=1, sticky='W')
et_ladoX.grid(row=1, column=0, sticky='E')
tb_ladoX.grid(row=1, column=1, sticky='W')
et_ladoY.grid(row=2, column=0, sticky='E')
tb_ladoY.grid(row=2, column=1, sticky='W')
et_angulo.grid(row=3, column=0, sticky='E')
tb_angulo.grid(row=3, column=1, sticky='W')
et_c.grid(row=0, column=2, sticky='E')
tb_c.grid(row=0, column=3, sticky='W')
et_Pu.grid(row=1, column=2, sticky='E')
tb_Pu.grid(row=1, column=3, sticky='W')
et_varX.grid(row=2, column=2, sticky='E')
tb_varX.grid(row=2, column=3, sticky='W')
et_varY.grid(row=3, column=2, sticky='E')
tb_varY.grid(row=3, column=3, sticky='W')
et_dia1.grid(row=0, column=4, sticky='E')
tb_dia1.grid(row=0, column=5, sticky='W')
et_dia2.grid(row=1, column=4, sticky='E')
tb_dia2.grid(row=1, column=5, sticky='W')
et_diaf.grid(row=2, column=4, sticky='E')
tb_diaf.grid(row=2, column=5, sticky='W')
et_recu.grid(row=3, column=4, sticky='E')
tb_recu.grid(row=3, column=5, sticky='W')
et_conf_ref.grid(row=4, column=4, sticky='E')
tb_conf_ref.grid(row=4, column=5, sticky='W')
f_datos.grid(columnspan=3)
btt_graficar.grid(row=1, column=0)
btt_clear.grid(row=1, column=1)
btt_salir.grid(row=1, column=2)
f_graficas.grid(columnspan=3, sticky='NSEW')
f_hor.grid(rowspan=2, sticky='NSEW')
f_3D.grid(row=0, column=1, sticky='NSEW')
f_ver.grid(row=1, column=1, sticky='NSEW')

root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)
root.columnconfigure(2, weight=1)
root.rowconfigure(2, weight=1)
f_graficas.columnconfigure(0, weight=1)
f_graficas.columnconfigure(1, weight=1)
f_graficas.rowconfigure(0, weight=1)
f_graficas.rowconfigure(1, weight=1)

root.mainloop()
