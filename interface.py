import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
# import graficas


def tkwindow(graf):
    fig = graf.fig

    ventana = tk.Tk()
    ventana.title("ColJuCris")
    ventana.geometry("800x800")

    frame1 = tk.Frame(ventana)

    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)

    canvas = FigureCanvasTkAgg(fig, frame1)
    canvas.draw()
    canvas.mpl_connect("key_press_event", on_key_press)
    toolbar = NavigationToolbar2Tk(canvas, frame1)
    toolbar.update()

    tx_print = tk.StringVar()

    etiqueta = tk.Label(ventana, text="ingres su texto", bg="lightblue", padx=50)
    et_print = tk.Label(ventana, textvariable=tx_print)

    tb_print = tk.Entry(ventana, textvariable=tx_print, width=10, borderwidth=3)

    def bf_print():
        et_print["text"] = tb_print.get()

    def bf_clear():
        tb_print.delete(0, tk.END)

    btt_print = tk.Button(ventana, text="print", command=bf_print)
    btt_clear = tk.Button(ventana, text="clear", command=bf_clear)
    btt_quit = tk.Button(ventana, text="salir", command=quit)

    etiqueta.pack(side=tk.TOP)
    etiqueta.grid(row=0, column=0, columnspan=3)
    et_print.grid(row=3, column=1, columnspan=3)
    tb_print.grid(row=1, column=1, padx=2, pady=2)
    btt_print.grid(row=2, column=0, padx=3, pady=3)
    btt_clear.grid(row=2, column=2, padx=3, pady=3)
    frame1.grid(row=4, column=0, columnspan=3)
    # canvas.get_tk_widget().grid(row=4, column=0, columnspan=3)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    btt_quit.grid(row=5, column=0, columnspan=3)

    ventana.mainloop()
    return ventana
