# Path: riana/riana_ui_plot.py
# -*- coding: utf-8 -*-
# RIANA GUI

import tkinter as tk
from tkinter import ttk

class Frame3(ttk.Frame):
    def __init__(self, container):
        super().__init__(container)

        self.labelC = ttk.Label(self, text="This is on Frame Three")
        self.labelC.grid(column=1, row=1)
