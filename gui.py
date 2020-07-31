import numpy as np
import Seccion as s
# from Seccion import calc_ang
import curvas
import graficas
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler

np.set_printoptions(precision=3, suppress=True)


class mainWin():

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("ColJuCris")
        self.root.resizable(0, 0)
        tk.Label(self.root, text="Bienvenido\n" +
                 "Abra la aplicación que desee usar").pack()
        tk.Button(self.root, text='Columna',
                  command=self.opencol).pack()

    def opencol(self):
        columnWin(self.root)


class columnWin():

    win_id = 0

    def __init__(self, parent):
        self.parent = parent
        columnWin.win_id += 1
        self.col_id = columnWin.win_id

        self.win = tk.Toplevel(parent)
        self.win.title(f"Columna {self.col_id}")
        self.win.geometry("750x550")
        self.datos = datosFrame(self)
        self.datos.frame.grid(row=0, column=0)
        self.set_col()

        chart_btts = chartButtonsFrame(self)
        chart_btts.frame.grid(row=1, column=0)

        self.f_refchart = chartFrame(self.win, is_3D=False)
        self.f_refchart.grid(row=0, column=1, rowspan=2)

        self.win.columnconfigure(1, weight=1)
        self.win.rowconfigure(0, weight=1)
        self.win.rowconfigure(1, weight=1)

        # self.win.resizable(0, 0)
        self.win.transient(master=parent)
        self.parent.wait_window(self.win)

    def set_col(self):
        d = self.datos
        sec = s.seccion(d.cv_fc.get(), d.cv_ladoX.get(), d.cv_ladoY.get())
        var1 = s.seccion_varilla(d.cv_dia1.get(), 420, 200000)
        var2 = s.seccion_varilla(d.cv_dia2.get(), 420, 200000)
        varf = s.seccion_varilla(d.cv_diaf.get(), 420, 200000)
        sec.set_ref_sim(d.cv_varX.get(), d.cv_varY.get(), var1, var2, varf,
                        d.cv_recu.get(), d.cv_conf_ref.get())
        self.sec = sec


class datosFrame():

    def __init__(self, colwin):
        self.frame = tk.Frame(colwin.win)
        self.colwin = colwin

        btt_defcol = tk.Button(self.frame, text='Actualizar\ncolumna',
                               command=self.update_col)

        self.cv_fc = tk.DoubleVar(value=28)
        self.cv_ladoX = tk.DoubleVar(value=0.6)
        self.cv_ladoY = tk.DoubleVar(value=0.3)
        self.cv_varX = tk.IntVar(value=4)
        self.cv_varY = tk.IntVar(value=4)
        self.cv_dia1 = tk.StringVar(value='#5')
        self.cv_dia2 = tk.StringVar(value='#5')
        self.cv_diaf = tk.StringVar(value='#3')
        self.cv_conf_ref = tk.StringVar(value='esquinero')
        self.cv_recu = tk.DoubleVar(value=.04)

        lb_fc = tk.Label(self.frame, text="f'c [MPa]=")
        lb_ladoX = tk.Label(self.frame, text="Lado X [m]=")
        lb_ladoY = tk.Label(self.frame, text="Lado Y [m]=")
        lb_varX = tk.Label(self.frame, text="Varillas en X=")
        lb_varY = tk.Label(self.frame, text="Varillas en Y=")
        lb_dia1 = tk.Label(self.frame, text="Ø_1 [#/8'']=")
        lb_dia2 = tk.Label(self.frame, text="Ø_2 [#/8'']=")
        lb_diaf = tk.Label(self.frame, text="Ø_fleje [#/8'']=")
        lb_conf_ref = tk.Label(self.frame, text="Op. Refuerzo=")
        lb_recu = tk.Label(self.frame, text="Recubrimiento [cm]=")

        tb_fc = tk.Entry(self.frame, width=5, borderwidth=3, textvariable=self.cv_fc)
        tb_ladoX = tk.Entry(self.frame, width=5, borderwidth=3, textvariable=self.cv_ladoX)
        tb_ladoY = tk.Entry(self.frame, width=5, borderwidth=3, textvariable=self.cv_ladoY)
        tb_varX = tk.Entry(self.frame, width=2, borderwidth=3, textvariable=self.cv_varX)
        tb_varY = tk.Entry(self.frame, width=2, borderwidth=3, textvariable=self.cv_varY)
        tb_dia1 = tk.Entry(self.frame, width=3, borderwidth=3, textvariable=self.cv_dia1)
        tb_dia2 = tk.Entry(self.frame, width=3, borderwidth=3, textvariable=self.cv_dia2)
        tb_diaf = tk.Entry(self.frame, width=3, borderwidth=3, textvariable=self.cv_diaf)
        tb_conf_ref = tk.Entry(self.frame, width=10, borderwidth=3,
                               textvariable=self.cv_conf_ref)
        tb_recu = tk.Entry(self.frame, width=5, borderwidth=3, textvariable=self.cv_recu)

        lb_fc.grid(row=0, column=0, sticky='E')
        tb_fc.grid(row=0, column=1, sticky='W')
        lb_ladoX.grid(row=1, column=0, sticky='E')
        tb_ladoX.grid(row=1, column=1, sticky='W')
        lb_ladoY.grid(row=2, column=0, sticky='E')
        tb_ladoY.grid(row=2, column=1, sticky='W')
        lb_varX.grid(row=3, column=0, sticky='E')
        tb_varX.grid(row=3, column=1, sticky='W')
        lb_varY.grid(row=4, column=0, sticky='E')
        tb_varY.grid(row=4, column=1, sticky='W')
        lb_dia1.grid(row=5, column=0, sticky='E')
        tb_dia1.grid(row=5, column=1, sticky='W')
        lb_dia2.grid(row=6, column=0, sticky='E')
        tb_dia2.grid(row=6, column=1, sticky='W')
        lb_diaf.grid(row=7, column=0, sticky='E')
        tb_diaf.grid(row=7, column=1, sticky='W')
        lb_recu.grid(row=8, column=0, sticky='E')
        tb_recu.grid(row=8, column=1, sticky='W')
        lb_conf_ref.grid(row=9, column=0, columnspan=2)
        tb_conf_ref.grid(row=10, column=0, columnspan=2, sticky='EW')
        btt_defcol.grid(row=11, column=0, columnspan=2)

    def update_col(self):
        self.colwin.set_col()
        self.colwin.f_refchart.draw()


class chartButtonsFrame():

    def __init__(self, colwin):
        self.frame = tk.Frame(colwin.win)
        self.colwin = colwin
        self.sec = colwin.sec
        self.col_id = colwin.col_id

        self.bool_3D = False
        self.bool_ver = False
        self.bool_hor = False
        # self.w_3D = None
        # self.w_ver = None
        # self.w_hor = None

        f_datos = tk.Frame(self.frame)
        f_datos.grid(row=0, column=0, sticky='EW')

        self.cv_asol = tk.DoubleVar(value=45)
        self.cv_Pu = tk.DoubleVar(value=100)
        # cv_c = tk.DoubleVar(value=.01)

        tk.Label(f_datos, text="Ángulo [°]=")\
            .grid(row=0, column=0, sticky='E')
        tk.Label(f_datos, text="Pu [KN]=")\
            .grid(row=0, column=2, sticky='E')
        # lb_c = tk.Label(f_datos, text="c [m]=")
        # lb_c.grid(row=0, column=2, sticky='E')

        tk.Entry(f_datos, width=3, borderwidth=3, textvariable=self.cv_asol)\
            .grid(row=0, column=1, sticky='W')
        tk.Entry(f_datos, width=6, borderwidth=3, textvariable=self.cv_Pu)\
            .grid(row=0, column=3, sticky='W')
        # tb_c = tk.Entry(f_datos, width=5, borderwidth=3, textvariable=cv_c)
        # tb_c.grid(row=0, column=3, sticky='W')

        f_datos.columnconfigure(0, weight=1)
        f_datos.columnconfigure(1, weight=1)
        f_datos.columnconfigure(2, weight=1)
        f_datos.columnconfigure(3, weight=1)

        f_buttons = tk.Frame(self.frame)
        f_buttons.grid(row=1, column=0)

        tk.Label(f_buttons, text="Dibujar").grid(row=0, column=0)
        tk.Label(f_buttons, text="Limpiar").grid(row=0, column=1)
        tk.Label(f_buttons, text="Cerrar").grid(row=0, column=2)

        tk.Button(f_buttons, text='3D', command=self.draw_3D)\
            .grid(row=1, column=0)
        tk.Button(f_buttons, text='3D', command=self.clear_3D)\
            .grid(row=1, column=1)
        tk.Button(f_buttons, text='X', command=self.close_3D)\
            .grid(row=1, column=2)
        tk.Button(f_buttons, text='Vertical', command=self.draw_ver)\
            .grid(row=2, column=0)
        tk.Button(f_buttons, text='Vertical', command=self.clear_ver)\
            .grid(row=2, column=1)
        tk.Button(f_buttons, text='X', command=self.close_ver)\
            .grid(row=2, column=2)
        tk.Button(f_buttons, text='FatCircle', command=self.draw_fatcircle)\
            .grid(row=3, column=0)
        tk.Button(f_buttons, text='Por ángulo\nsolicitacion', command=self.draw_ang_sol)\
            .grid(row=4, column=0)
        tk.Button(f_buttons, text='Por ángulo\neje neutro', command=self.draw_ang_eje)\
            .grid(row=5, column=0)
        tk.Button(f_buttons, text='Horizontales', command=self.clear_hor)\
            .grid(row=3, column=1, rowspan=3, sticky='NS')
        tk.Button(f_buttons, text='X', command=self.close_hor)\
            .grid(row=3, column=2, rowspan=3, sticky='NS')
        tk.Button(f_buttons, text='Todo', command=self.draw_all)\
            .grid(row=6, column=0)
        tk.Button(f_buttons, text='Todo', command=self.clear_all)\
            .grid(row=6, column=1)
        tk.Button(f_buttons, text='X', command=self.close_all)\
            .grid(row=6, column=2)

    def draw_3D(self):
        if not self.bool_3D:
            self.w_3D = chartToplevel(self.colwin.win, is_3D=True)
            self.w_3D.geometry("400x400")
            self.w_3D.title(f"Columna {self.col_id}, P(My,Mx) [KN, m]")
            self.bool_3D = True
        pmm = np.array(curvas.vertical(self.sec, self.cv_asol.get()))
        self.w_3D.addPMM(pmm, '%.1f°' % self.cv_asol.get())
        self.w_3D.draw()

    def draw_ver(self):
        if not self.bool_ver:
            self.w_ver = chartToplevel(self.colwin.win, is_3D=False)
            self.w_ver.geometry("400x400")
            self.w_ver.title(f"Columna {self.col_id}, P(M) [KN, m]")
            self.bool_ver = True
        pmm = np.array(curvas.vertical(self.sec, self.cv_asol.get()))
        self.w_ver.addPM(pmm, '%.1f°' % self.cv_asol.get())
        self.w_ver.draw()

    def draw_fatcircle(self):
        if not self.bool_hor:
            self.w_hor = chartToplevel(self.colwin.win, is_3D=False)
            self.w_hor.geometry("400x400")
            self.w_hor.title(f"Columna {self.col_id}, My(Mx) [KN, m]")
            self.bool_hor = True
        parlem = np.array(curvas.horizontal(self.sec, self.cv_Pu.get(), metodo='fatcircle'))
        self.w_hor.add2D(parlem, "parlem")
        self.w_hor.draw()

    def draw_ang_sol(self):
        if not self.bool_hor:
            self.w_hor = chartToplevel(self.colwin.win, is_3D=False)
            self.w_hor.geometry("400x400")
            self.w_hor.title(f"Columna {self.col_id}, My(Mx) [KN, m]")
            self.bool_hor = True
        solici = np.array(curvas.horizontal(self.sec, self.cv_Pu.get(), metodo='ang_sol'))
        self.w_hor.add2D(solici, "ang sol")
        self.w_hor.draw()

    def draw_ang_eje(self):
        if not self.bool_hor:
            self.w_hor = chartToplevel(self.colwin.win, is_3D=False)
            self.w_hor.geometry("400x400")
            self.w_hor.title(f"Columna {self.col_id}, My(Mx) [KN, m]")
            self.bool_hor = True
        column = np.array(curvas.horizontal(self.sec, self.cv_Pu.get(), metodo='ang_eje'))
        self.w_hor.add2D(column, "ang eje")
        self.w_hor.draw()

    def draw_all(self):
        self.draw_3D()
        self.draw_ver()
        self.draw_fatcircle()
        self.draw_ang_sol()
        self.draw_ang_eje()

    def clear_3D(self):
        if self.bool_3D:
            self.w_3D.cla()
            self.w_3D.draw()

    def clear_ver(self):
        if self.bool_ver:
            self.w_ver.cla()
            self.w_ver.draw()

    def clear_hor(self):
        if self.bool_hor:
            self.w_hor.cla()
            self.w_hor.draw()

    def clear_all(self):
        self.clear_3D()
        self.clear_ver()
        self.clear_hor()

    def close_3D(self):
        if self.bool_3D:
            self.w_3D.destroy()
            self.bool_3D = False

    def close_ver(self):
        if self.bool_ver:
            self.w_ver.destroy()
            self.bool_ver = False

    def close_hor(self):
        if self.bool_hor:
            self.w_hor.destroy()
            self.bool_hor = False

    def close_all(self):
        self.close_3D()
        self.close_ver()
        self.close_hor()


class chartCanvas(graficas.grafica):

    def __init__(self, parent, is_3D):
        graficas.grafica.__init__(self, is_3D)
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.draw()
        self.canvas.mpl_connect("key_press_event", self.on_key_press)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

    def on_key_press(self, event):
        # print(f"you pressed {event.key}")
        key_press_handler(event, self.canvas, self.toolbar)

    def draw(self):
        self.axe.legend()
        self.canvas.draw()
        print("canvas dibujado")


class chartToplevel(tk.Toplevel, chartCanvas):

    def __init__(self, parent, is_3D):
        tk.Toplevel.__init__(self, parent)
        chartCanvas.__init__(self, self, is_3D)
        # parent.wait_window(self)


class chartFrame(tk.Frame, chartCanvas):

    def __init__(self, parent, is_3D):
        tk.Frame.__init__(self, parent)
        chartCanvas.__init__(self, self, is_3D)


# def bf_salir():
#     quit()

# btt_salir = tk.Button(root, text="salir", command=bf_salir)
# btt_salir.grid(row=0, column=3)

# f_graficas = tk.Frame(root)
#
# f_3D = chart_frame(f_graficas, is_3D=True)
# f_ver = chart_frame(f_graficas, is_3D=False)
# f_hor = chart_frame(f_graficas, is_3D=False)
#
# f_graficas.grid(row=3, column=0, sticky='NSEW')
