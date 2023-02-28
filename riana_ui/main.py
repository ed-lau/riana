# Path: riana/riana_ui.py
# -*- coding: utf-8 -*-
# RIANA GUI

import sys
from typing import NamedTuple
import tkinter as tk
from tkinter import ttk
import sv_ttk
import queue
import threading
import logging
import tqdm


from riana import __version__
from riana import riana_integrate
from riana_ui.riana_ui_integrate import Frame1
from riana_ui.riana_ui_model import Frame2
from riana_ui.riana_ui_plot import Frame3


class TextHandler():
    pass


class Menubar(tk.Menu):
    """ Menu bar """
    def __init__(self, parent):
        tk.Menu.__init__(self, parent)

        file_menu = tk.Menu(self, tearoff=False)
        self.add_cascade(label="File",underline=0, menu=file_menu)
        file_menu.add_command(label="Exit RIANA", underline=1, command=self.quit)

        about_menu = tk.Menu(self, tearoff=False)
        self.add_cascade(label="About", underline=0, menu=about_menu)
        about_menu.add_command(label="About RIANA", underline=1, command=self.open_about)

    def open_about(self):
        top = tk.Toplevel(self)
        top.geometry("250x250")
        top.title("About RIANA")
        ttk.Label(top, text=f'RIANA {__version__}', font=('Mistral 18 bold')).place(x=50, y=50)

    def quit(self):
        sys.exit(0)

class Application(tk.Tk):
    """ Main application. """

    def __init__(self):
        """ Initialize the application. """
        super().__init__()

        self.title(f'RIANA {__version__}')
        self.geometry('1200x700')
        self.minsize('700', '700')

        # Create the main frame and notebook layout
        self.notebook = ttk.Notebook(self)

        self.Frame1 = Frame1(self.notebook)
        self.Frame2 = Frame2(self.notebook)
        self.Frame3 = Frame3(self.notebook)

        self.notebook.add(self.Frame1, text='Integrate')
        self.notebook.add(self.Frame2, text='Model')
        self.notebook.add(self.Frame3, text='Plot')

        self.notebook.pack(fill='both', expand=True)

        # ttk theme options
        style = ttk.Style()
        style.theme_use('aqua')
        # # style.configure("TButton", padding=6, relief="flat",
        # #                 background="#ccc")
        # # style.map("TButton", background=[('active', 'white')])
        # # style.configure("TFrame",
        # #                 background="#121212",
        # #                 foreground="#121212",
        # #                 borderwidth=0,
        # #                 relief=FLAT)
        #
        # # style.configure("TLabel", background="#ccc")
        # #
        # style.configure("TNotebook", background="#111111")
        # # style.configure("TNotebook.Tab", background="#ccc")
        # style.map("TNotebook.Tab", background=[('active', 'white')])
        # # style.configure("TEntry", background="#ccc")
        # sv_ttk.use_dark_theme()



def main():


    app = Application()

    # Menu
    menu = tk.Menu(app)
    app.config(menu=Menubar(app))

    app.mainloop()

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


if __name__ == "__main__":
    main()