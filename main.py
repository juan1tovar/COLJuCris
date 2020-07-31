import sys
import debugexc
import gui


sys.excepthook = debugexc.info


def main():
    app = gui.mainWin()
    app.root.mainloop()
    return(0)


if __name__ == '__main__':
    main()
