# Path: riana/riana_ui_model.py
# -*- coding: utf-8 -*-
# RIANA GUI

import tkinter as tk
from tkinter import ttk, filedialog, FLAT, BOTH, LEFT, TOP, END, BOTTOM, X, StringVar, IntVar, DoubleVar
from riana import riana_fit, models
from typing import NamedTuple
import sys
import rx
import os
import ast
import pandas as pd
from pandastable import Table, TableModel, config
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class ModelingVars(NamedTuple):
    """ Variables for modeling. """
    riana_path: list[str]
    model: str
    label: str
    aa: str
    kp: bool
    kr: float
    rp: float
    q_value: float
    depth: int
    ria: float
    out: str
    thread: int
    plotcurves: bool
    gui: bool

class Frame2(ttk.Frame):
    def __init__(self, container):
        super().__init__(container)


        self.labelB = ttk.Label(self, text="Step 3: RIANA Model")
        self.labelB.pack()

        self.riana_path = None
        self.model_output_path = None

        self.model = StringVar()
        self.model_options = ["simple", "guan", "fornasiero"]
        self.model.set("simple")

        self.label_type = StringVar()
        self.label_type_options = ['aa', 'hw', 'o18']
        self.label_type.set("aa")

        self.kp = DoubleVar()
        self.kp.set(0.5)

        self.kr = DoubleVar()
        self.kr.set(0.05)

        self.rp = DoubleVar()
        self.rp.set(10)

        self.final_ria = tk.DoubleVar()
        self.final_ria.set(0.1)

        self.q_value = tk.DoubleVar()
        self.q_value.set(0.01)

        self.aa = StringVar()
        self.aa.set("K")

        self.create_tab2_widgets()

        # sys.stdout = TextRedirector(self.master.master.console, "stdout")
        # sys.stderr = TextRedirector(self.master.master.console, "stderr")


    def handle_model_selection(self, value):
        # disable kp, kr, rp if simple model is selected
        self.kp_entry.config(state='normal')
        self.kr_entry.config(state='normal')
        self.rp_entry.config(state='normal')

        if self.model.get() == "simple":
            self.kp_entry.config(state='disabled')
            self.kr_entry.config(state='disabled')
            self.rp_entry.config(state='disabled')

        elif self.model.get() == 'guan':
            self.kp_entry.config(state='normal')
            self.kr_entry.config(state='disabled')
            self.rp_entry.config(state='disabled')

    def handle_kr_selection(self, value):

        try:
            float(self.kr.get())
        except ValueError:
            self.kr.set(0.05)
            return

        if self.kr.get() > 1:
            self.kr.set(1)
        elif self.kr.get() < 0:
            self.kr.set(0)

        print(self.kr.get())
        self.kr_label_text.set(f'kr: {self.kr.get()}')

    def create_tab2_widgets(self):
        # Create left and right frames
        tab2_left_frame = ttk.LabelFrame(self,
                                    width=200,
                                    height=400,
                                    relief='flat',
                                    # borderwidth=1,
                                    text="Input",
                                    )
        tab2_left_frame.pack(side='left', fill='both', expand=False)

        self.tab2_right_frame = ttk.LabelFrame(self,
                                     width=650,
                                     height=400,
                                     relief='flat',
                                     # borderwidth=1,
                                     text="Output",
                                     )
        self.tab2_right_frame.pack(side="right", fill='both', expand=True)

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
        ttk.Label(tab2_left_frame, text="Model Options").pack()
        ttk.Separator(tab2_left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Select model to use ----
        self.model_label = ttk.Label(tab2_left_frame,
                                       text=f'Model to use: {self.model.get()}',
                                       )
        self.model_label.pack()
        self.model_entry = ttk.OptionMenu(tab2_left_frame,
                                           self.model,
                                           self.model.get(),
                                          *self.model_options,
                                           command=self.handle_model_selection

                                        )
        self.model_entry.pack(side="top")

        # ---- Select label type ----
        self.label_type_label = ttk.Label(tab2_left_frame,
                                        text=f'Label type: {self.label_type.get()}',
                                        )
        self.label_type_label.pack()
        self.label_type_entry = ttk.OptionMenu(tab2_left_frame,
                                             self.label_type,
                                                self.label_type.get(),
                                               *self.label_type_options,
                                               # command=self.handle_label_type_selection
                                                )
        self.label_type_entry.pack(side="top")

        # ---- Select aa ----
        self.aa_label = ttk.Label(tab2_left_frame,
                                        text=f'aa: {self.aa.get()}',
                                        )
        self.aa_label.pack()
        self.aa_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.aa,
                                                width=10,
                                                )
        self.aa_entry.pack(side="top")

        # ---- Select kp ----
        self.kp_label = ttk.Label(tab2_left_frame,
                                        text=f'kp: {self.kp.get()}',
                                        )
        self.kp_label.pack()
        self.kp_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.kp,
                                                width=10,
                                                )
        self.kp_entry.pack(side="top")

        # ---- Select kr ----
        self.kr_label_text = StringVar()
        self.kr_label_text.set(f'kr: {self.kr.get()}')
        self.kr_label = ttk.Label(tab2_left_frame,
                                        text=self.kr_label_text.get(),
                                        )
        self.kr_label.pack()
        self.kr_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.kr,
                                                width=10,
                                  # command=self.handle_kr_selection,
                                                )
        self.kr_entry.bind('<KeyRelease>', self.handle_kr_selection)
        self.kr_entry.pack(side="top")

        # ---- Select rp ----
        self.rp_label = ttk.Label(tab2_left_frame,
                                        text=f'rp: {self.rp.get()}',
                                        )
        self.rp_label.pack()
        self.rp_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.rp,
                                                width=10,
                                                )
        self.rp_entry.pack(side="top")

        # ---- Select final ria ----
        self.final_ria_label = ttk.Label(tab2_left_frame,
                                        text=f'final ria: {self.final_ria.get()}',
                                        )
        self.final_ria_label.pack()
        self.final_ria_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.final_ria,
                                                width=10,
                                                )
        self.final_ria_entry.pack(side="top")

        # ---- Select q-value ----
        self.q_value_label = ttk.Label(tab2_left_frame,
                                        text=f'q-value: {self.q_value.get()}',
                                        )
        self.q_value_label.pack()
        self.q_value_entry = ttk.Entry(tab2_left_frame,
                                                textvariable=self.q_value,
                                                width=10,
                                                )
        self.q_value_entry.pack(side="top")



        # link to RIANA model
        self.model_run_button = ttk.Button(tab2_left_frame,
                                     text="Run RIANA Model",
                                     command=self.handle_model_button,
                                     width=20,
                                     state="disabled"
                                     )
        self.model_run_button.pack(side=BOTTOM)





        ttk.Label(self.tab2_right_frame, text="Input").pack()
        ttk.Separator(self.tab2_right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w', )
        # self.select_riana_integration_files.bind("<Button-1>", self.print_selected_files)
        self.input_frame = ttk.Frame(self.tab2_right_frame,
                                        width=600,
                                        height=300,
                                        )
        self.input_frame.pack()

        # ---- Section separator ----
        ttk.Label(self.tab2_right_frame, text="Results").pack()
        ttk.Separator(self.tab2_right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')



        # ---- Using pandastable ----

        self.model_result_view = ttk.Frame(self.tab2_right_frame,
                                     width=600,
                                     height=300,
                                     )
        self.model_result_view.pack(fill=BOTH, expand=1)
        # self.load_result_table()
        self.model_result_view.pack_forget()

        # ---- Inspect view from selected row on the data ----
        self.model_inspect_view = ttk.Frame(self.tab2_right_frame,
                                     width=600,
                                     height=200,
                                     )
        # self.inspect_view.pack(fill=BOTH, expand=1)

        self.model_inspect_text = ttk.Label(self.model_inspect_view, text='')
        self.model_inspect_text.pack()

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
        self.master.master.console.insert('end', self.riana_path)
        print(self.riana_path)

        # enable button
        if self.riana_path is not None and self.model_output_path is not None:
             self.model_run_button["state"] = "normal"

        # For each file, create a new label in the right frame to display the name, as well as a button to remove it
        for each_path in self.riana_path:
            ttk.Label(self.input_frame, text=each_path).pack(side="top", anchor="w")


    # output folder dialog
    def model_output_folder_dialog(self):
        """ Returns a selected directory name. """
        self.model_output_path = filedialog.askdirectory(initialdir='./out/ui_test/',
                                                   title='Select output folder',
                                                   )
        self.master.master.console.insert(END, self.model_output_path)
        print(self.model_output_path)

        # enable button
        if self.riana_path is not None and self.model_output_path is not None:
             self.model_run_button["state"] = "normal"

    # dummy progress bar
    def progress(self):
        self.prog_bar = ttk.Progressbar(
            self.master, orient="horizontal",
            length=400, mode="indeterminate"
        )
        self.prog_bar.pack(side=BOTTOM, fill=X)

    def handle_model_button(self):
        """
        Handle the model button click, run model fit
        :return:
        """
        rx.empty().subscribe(
            on_completed=self.on_click,
            scheduler=self.master.master.pool_scheduler
        )

    def on_click(self) -> None:
        """
        Handle the model button click, run model fit
        :return:
        """
        self.progress()
        self.prog_bar.start()

        # Disable all buttons
        self.model_run_button["text"] = "Running..."
        self.model_run_button["state"] = "disabled"
        self.master.tab(0, state="disabled")

        self.master.master.status_label["text"] = "Status: Running..."

        modeling_vars = ModelingVars(
            riana_path=self.riana_path,
            model=self.model.get(),
            label=self.label_type.get(),
            aa='KR',
            kp=0.5,
            kr=0.05,
            rp=10,
            q_value=0.01,
            depth=9,
            ria=0.08,
            out=self.model_output_path,
            thread=1,
            plotcurves=False,
            gui=True,
        )

        riana_fit.fit_all(modeling_vars)

        # Callback for task complete
        self.master.master.after(5, self.on_task_complete)

    def on_task_complete(self):
        """
        After the model run is complete, stop the progress bar and update the status
        :return:
        """
        self.prog_bar.stop()
        self.prog_bar.destroy()
        # self.tab2_output.insert(END, msg)
        # with open(os.path.join(self.out, 'tmt.log'), 'r') as f:
        #     self.output.insert(INSERT, f.read())

        # self.queue.task_done()

        # Clear the result view when the button is pressed
        self.model_result_view.pack_forget()
        self.model_inspect_view.pack_forget()


        # Disable all buttones
        self.model_run_button["text"] = "Model"
        self.model_run_button["state"] = "normal"
        self.master.tab(0, state="normal")

        self.master.master.status_label["text"] = "Status: Finished"

        # load the model result table
        self.load_model_result_table()
        self.model_result_view.pack(fill=BOTH, expand=1)

        return None  # runs on main thread


    def load_model_result_table(self):
        """
        Function to load the model result table
        :return:
        """
        try:  # try to read in data from RIANA fit
            fitted_file_path = os.path.join(self.model_output_path, 'riana_fit_peptides.txt')
            # print(output_file_path)
            self.master.master.console.insert(END, f"Reading in data from {fitted_file_path}\n")
            df = pd.read_csv(fitted_file_path, sep='\t', index_col='concat').reset_index()
            # df = df.rename(columns={'Unnamed: 0': 'sequence'})

        except:  # if not, read in sample data
            df = TableModel.getSampleData()

        self.model_result_table = Table(self.model_result_view,
                                  dataframe=df,
                                  showtoolbar=False,
                                  showstatusbar=False,
                                  )

        # self.result_table.setTheme('default')
        self.model_result_table.editable = False

        # self.result_table.cols.colheader.bgcolor = '#ECECEC'

        self.model_result_table.show()
        options = {'font': 'SF Pro Text',
                   'fontsize': 13,
                   'rowheight': 20,
                   'grid_color': '#FFFFFF',
                   'rowselectedcolor': '#B7D6F0',
                   'colselectedcolor': '#B7D6F0',
                   'boxoutlinecolor': '#FFFFFF',
                   'cellbackgr': '#FFFFFF',
                   }
        config.apply_options(options, self.model_result_table)

        # Set header attributes after show()
        self.model_result_table.colheader.bgcolor = '#ECECEC'
        self.model_result_table.colheader.textcolor = 'black'
        self.model_result_table.rowheader.bgcolor = '#ECECEC'
        self.model_result_table.rowheader.textcolor = 'black'

        def handle_model_table_left_click(event):

            rowclicked_single = self.model_result_table.get_row_clicked(event)
            # print(rowclicked_single)

            # Clear the canvas
            try:
                self.canvas.get_tk_widget().pack_forget()
                for item in self.canvas.get_tk_widget().find_all():
                    self.canvas.get_tk_widget().delete(item)
                self.toolbar.pack_forget()
            except AttributeError:
                pass

            self.model_result_table.setSelectedRow(rowclicked_single)
            self.model_result_table.redraw()

            # Refresh the inspect view every time there is a new selection
            self.model_inspect_view.pack_forget()
            self.model_inspect_view.pack(fill=BOTH, expand=1)
            self.trace_plot = None

            # Get the selected row data
            selected_data = self.model_result_table.getSelectedRowData()
            # print(selected_data['pep_id'])
            print(selected_data)

            # Display the text on the interface
            model_out_text= f'Row: {self.model_result_table.getSelectedRow()}, ' \
                      f'Peptide: {selected_data["concat"].values[0]}, ' \
                      #f'Protein: {selected_data["protein id"].values[0]} ' \
                      # f'Retention time: {intensities_df_subset["rt"].values[0]}, ' \
                      #  f'Mass: {intensities_df_subset["m0"].values[0]}'

            self.model_inspect_text.config(text=model_out_text)
            self.model_inspect_text.pack()

            # Turn the literal string into a list of floats
            tp = ast.literal_eval(selected_data['t'].values[0])
            tp = np.array([float(i) for i in tp])

            fs = ast.literal_eval(selected_data['fs'].values[0])
            fs = np.array([float(i) for i in fs])
            print(tp)
            print(fs)

            # Plot the model using the same function as the CLI
            self.fig = models.plot_model(protein=selected_data["protein id"].values[0],
                                    peptide=selected_data["concat"].values[0],
                                    k_deg=selected_data["k_deg"].values[0],
                                    r_squared=selected_data["R_squared"].values[0],
                                    sd=selected_data["sd"].values[0],
                                    t_series=tp,
                                    fs_series=fs,
                                    start_time=0,
                                    end_time=np.max(tp),
                                    model_to_use=models.one_exponent,
                                    model_pars={'k_p': 1,
                                                'k_r': 1,
                                                'r_p': 1})


            # # Display the plots on the interface
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.model_inspect_view)
            self.canvas.draw()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.model_inspect_view)
            self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            self.toolbar.update()
            self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.model_result_table.bind("<Button-1>", handle_model_table_left_click)
