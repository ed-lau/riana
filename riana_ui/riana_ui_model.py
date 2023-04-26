# Path: riana/riana_ui_model.py
# -*- coding: utf-8 -*-
# RIANA GUI

import tkinter as tk
from tkinter import ttk, filedialog, FLAT, BOTH, LEFT, TOP, END, BOTTOM, X


class Frame2(ttk.Frame):
    def __init__(self, container):
        super().__init__(container)

        self.labelB = ttk.Label(self, text="Step 3: RIANA Model")
        self.labelB.pack()

        self.riana_path = None

        self.create_tab2_widgets()

    def create_tab2_widgets(self):
        # Create left and right frames
        tab2_left_frame = ttk.LabelFrame(self,
                                    width=200,
                                    height=400,
                                    relief='flat',
                                    # borderwidth=1,
                                    text="Input",
                                    )
        tab2_left_frame.pack(side='left', fill='both', expand=True)

        tab2_right_frame = ttk.LabelFrame(self,
                                     width=650,
                                     height=400,
                                     relief='flat',
                                     # borderwidth=1,
                                     text="Output",
                                     )
        tab2_right_frame.pack(side="right", fill='both', expand=True)

        # ---- Section separator ----
        ttk.Label(tab2_left_frame, text="Select Input File(s)").pack()
        ttk.Separator(tab2_left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w', )

        # ---- Select Percolator ID file ----
        self.select_riana_integration_files = ttk.Button(tab2_left_frame,
                                            text="Select RIANA Integration file(s)",
                                            command=self.riana_files_dialog,
                                            width=20,
                                            state="normal"
                                            )
        self.select_riana_integration_files.pack(side="top")

        # ---- Select output folder ----
        self.select_model_output = ttk.Button(tab2_left_frame,
                                        text="Select Output Folder",
                                        command=self.model_output_folder_dialog,
                                        width=20,
                                        state="normal"
                                        )
        self.select_model_output.pack(side="top")

        # ---- Section separator ----
        ttk.Label(tab2_left_frame, text="Integrate Options").pack()
        ttk.Separator(tab2_left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Display output ----
        self.tab2_output = tk.Text(tab2_right_frame,
                              width=400,
                              height=5,
                              )
        self.tab2_output.pack()
        self.tab2_output.insert('end', "Output will be displayed here")



        ttk.Separator(tab2_right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w', )
        # self.select_riana_integration_files.bind("<Button-1>", self.print_selected_files)

    # ---- Section separator ----
    def print_selected_files(self):
        if self.riana_path is not None:
            for each_path in self.riana_path:
                tk.LabelFrame(self, text=each_path).pack()

    # RIANA integration files dialog
    def riana_files_dialog(self):

        self.riana_path = filedialog.askopenfilenames(initialdir="./tests/local/fit/",
                                                      title="Select one or more RIANA integration files",)
        # filetypes = #(("jpeg files","*.jpg"), ("all files","*.*")) )
        self.tab2_output.insert('end', self.riana_path)
        print(self.riana_path)

        # enable button
        # if self.id_path is not None and self.mzml_path is not None:
        #     self.run_button["state"] = "normal"

    # output folder dialog
    def model_output_folder_dialog(self):
        """ Returns a selected directory name. """
        self.model_output_file = filedialog.askdirectory(initialdir='./out/ui_test/',
                                                   title='Select output folder',
                                                   )
        self.tab2_output.insert(END, self.model_output_file)
        print(self.model_output_file)