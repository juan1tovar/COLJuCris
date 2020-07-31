import sys
import tkinter as tk


def info(type, value, tb):
    if hasattr(sys, 'ps1') or not sys.stderr.isatty():
        # we are in interactive mode or we don't have a tty-like
        # device, so we call the default hook
        sys.__excepthook__(type, value, tb)
    else:
        import traceback
        import pdb
        # we are NOT in interactive mode, print the exception...
        print('antes de print traceback')
        traceback.print_exception(type, value, tb)
        # print
        print('antes de pdb')
        print('type:', type)
        print('value:', value)
        print('tb:', tb)
        # ...then start the debugger in post-mortem mode.
        # pdb.pm() # deprecated
        pdb.post_mortem(tb)  # more "modern"


sys.excepthook = info


class TkErrorCatcher:

    '''
    In some cases tkinter will only print the traceback.
    Enables the program to catch tkinter errors normally

    To use
    import tkinter
    tkinter.CallWrapper = TkErrorCatcher
    '''

    def __init__(self, func, subst, widget):
        self.func = func
        self.subst = subst
        self.widget = widget
        # self.win = None

    def __call__(self, *args):
        try:
            if self.subst:
                args = self.subst(*args)
            return self.func(*args)
        except SystemExit as msg:
            print('raise SystemExit exception')
            raise SystemExit(msg)
        except Exception:
            # print('raise exception')
            # print('Exception:', Exception)
            self.type, self.value, self.tb = sys.exc_info()

            self.ErrorBox()
            # raise

    def ErrorBox(self):
        root = self.widget
        while root.master is not None:
            root = root.master

        self.errorwin = tk.Toplevel(root)

        self.errorwin.title('Exception handler')
        tk.Label(self.errorwin, text=self.value).grid(columnspan=2)
        btt_debug = tk.Button(self.errorwin, text='Depurar', command=self.raiseException)
        btt_debug.grid(row=1, column=0)
        btt_skip = tk.Button(self.errorwin, text='Continuar', command=self.continuar)
        btt_skip.grid(row=1, column=1)

        self.errorwin.grab_set()
        self.errorwin.resizable(0, 0)
        root.wait_window(self.errorwin)

    def raiseException(self):
        self.errorwin.destroy()
        sys.excepthook(self.type, self.value, self.tb)

    def continuar(self):
        self.errorwin.destroy()


tk.CallWrapper = TkErrorCatcher
