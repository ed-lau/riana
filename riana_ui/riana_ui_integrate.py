# Path: riana/riana_ui.py
# -*- coding: utf-8 -*-
# RIANA GUI

from typing import NamedTuple
import threading
import sys
import tkinter as tk
import tqdm
from tkinter import ttk, filedialog, FLAT, BOTH, LEFT, TOP, END, BOTTOM
import queue
from riana import riana_integrate
import pandas as pd
from pandastable import Table, TableModel, config

class IntegrationVars(NamedTuple):
    """ Variables for integration. """
    id_path: str
    mzml_path: str
    out: str
    unique: bool
    q_value: float
    r_time: float
    mass_tol: int
    mass_defect: str
    thread: int
    write_intensities: bool
    sample: str
    iso: str


class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):

        self.widget.configure(state="normal")
        #self.widget.delete("1.0", "end-1l")  # delete previous text
        self.widget.delete("1.0", "end")  # delete previous text
        self.widget.insert("end", str, (self.tag,))
        self.widget.configure(state="disabled")
        # Autoscroll to the bottom
        self.widget.yview("end")

class Frame1(ttk.Frame):
    """ RIANA integrate """
    def __init__(self, container):
        """ Initialize the frame. """

        super().__init__(container)

        # Frame label
        self.labelA = ttk.Label(self, text="Step 1: RIANA Integrate")
        self.labelA.pack()


        # Mimicking argparse options (move this to frame1 later)
        self.id_path = None
        self.mzml_path = None
        self.out = 'out'

        self.unique = tk.BooleanVar()

        self.q_value = tk.DoubleVar()
        self.q_value.set(0.01)

        self.sample = tk.StringVar()
        self.sample.set('time1')

        self.r_time = tk.DoubleVar()
        self.r_time.set(0.5)

        self.mass_tol = tk.IntVar()
        self.mass_tol.set(25)

        self.mass_defect = tk.StringVar()
        self.mass_defect.set('D')
        self.mass_defect_options: list[str] = ['D', 'C13', 'SILAC']

        self.thread = 4
        self.write_intensities = False

        self.iso = tk.StringVar()
        self.iso.set('0,6')


        self.create_tab1_widgets()

    def update_q_value(self,
                       value: float):
        self.q_value.set(value)
        tqdm.tqdm.write(str(self.q_value.get()))
        self.q_value_label.config(text=f'Peptide FDR cutoff: {round(self.q_value.get(), 3)}')

    def update_r_time(self,
                       value: float):
        self.r_time.set(value)
        tqdm.tqdm.write(str(self.r_time.get()))
        self.r_time_label.config(text=f'Retention time range: {round(self.r_time.get(), 1)}')

    def update_sample(self, value):
        # self.sample.set(self.type_sample.get().replace("\n", ""))
        # self.type_sample.delete("1.0", END)
        # self.type_sample.set(self.sample.get().rstrip())
        # self.type_sample.delete("end-1c", END)
        # self.type_sample.xview_moveto(0.0)
        print(self.sample.get().rstrip())

    def update_isotopomers(self, value):
        print(self.iso.get())

    def create_tab1_widgets(self):

        # Create left and right frames
        left_frame = ttk.LabelFrame(self,
                                    width=200,
                                    height=400,
                                    relief=FLAT,
                                    #borderwidth=1,
                                    text="Input",
                                    )
        left_frame.pack(side=LEFT, fill=BOTH, expand=True)

        right_frame = ttk.LabelFrame(self,
                                     width=650,
                                     height=400,
                                     relief=FLAT,
                                     #borderwidth=1,
                                     text="Output",
                                     )
        right_frame.pack(side="right", fill=BOTH, expand=True)

        # ---- Section separator ----
        ttk.Label(left_frame, text="Select Input Files").pack()
        ttk.Separator(left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w',)

        # ---- Select Percolator ID file ----
        self.select_percolator = ttk.Button(left_frame,
                                            text="Select Percolator ID file",
                                            command=self.percolator_file_dialog,
                                            width=20,
                                            state="normal"
                                            )
        self.select_percolator.pack(side="top")

        # ---- Select path to mzML files ----
        self.select_mzml = ttk.Button(left_frame,
                                      text="Select Path to mzML files",
                                      command=self.browse_mzml,
                                      width=20,
                                      state="normal"
                                      )
        self.select_mzml.pack(side="top")


        # ---- Section separator ----
        ttk.Label(left_frame, text="Integrate Options").pack()
        ttk.Separator(left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Choose sample name ----
        self.sample_label = ttk.Label(left_frame,
                                      text='Sample name:',
                                      )
        self.sample_label.pack()
        self.type_sample = ttk.Entry(left_frame,
                                        textvariable=self.sample,
                                        width=20,
                                        )

        self.type_sample.pack(side="top")
        self.type_sample.bind('<KeyRelease>', self.update_sample)

        # ---- Select isotopomers ----
        self.isotopomer_label = ttk.Label(left_frame,
                                            text='Isotopomers:',
                                            )
        self.isotopomer_label.pack()
        self.type_isotopomer = ttk.Entry(left_frame,
                                            textvariable=self.iso,
                                            width=20,
                                            )
        self.type_isotopomer.pack(side="top")
        self.type_isotopomer.bind('<KeyRelease>', self.update_isotopomers)

        # ---- Select q value threshold ----
        self.q_value_label = ttk.Label(left_frame,
                                       text=f'Peptide FDR cutoff: {round(self.q_value.get(), 3)}',
                                       )
        self.q_value_label.pack()
        self.select_q_value = ttk.Scale(left_frame,
                                        from_=0.0,
                                        to=0.1,
                                        orient=tk.HORIZONTAL,
                                        variable=self.q_value,
                                        command=self.update_q_value,
                                        length=200,
                                        )
        self.select_q_value.pack(side="top")

        # ---- Select retention time threshold
        self.r_time_label = ttk.Label(left_frame,
                                       text=f'Retention time range: {round(self.r_time.get(), 3)}',
                                       )
        self.r_time_label.pack()
        self.select_r_time = ttk.Scale(left_frame,
                                        from_=0.1,
                                        to=10.0,
                                        orient=tk.HORIZONTAL,
                                        variable=self.r_time,
                                        command=self.update_r_time,
                                        length=200,
                                        )
        self.select_r_time.pack(side="top")

        # ---- Select unique peptide
        self.unique_label = ttk.Label(left_frame,
                                       text=f'Integrate only unique peptides:',
                                       )
        self.unique_label.pack()
        self.select_unique = ttk.Checkbutton(left_frame,
                                             variable=self.unique,
                                             onvalue=True,
                                             offvalue=False,
                                             command=lambda: print(self.unique.get()),

                                        )
        self.select_unique.pack(side="top")


        # ---- Section separator ----
        ttk.Label(left_frame, text="Advanced Options").pack()
        ttk.Separator(left_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Select mass defect ----
        self.mass_defect_label = ttk.Label(left_frame,
                                        text=f'Mass defect:',
                                        )
        self.mass_defect_label.pack()
        self.select_mass_defect = ttk.OptionMenu(left_frame,
                                                    self.mass_defect,
                                                    self.mass_defect_options[0],
                                                    *self.mass_defect_options,
                                                 )
        self.select_mass_defect.pack(side="top")

        # link to RIANA integrate
        self.run_button = ttk.Button(left_frame,
                                     text="Run RIANA Integrate",
                                     command=self.run_riana_integrate,
                                     width=20,
                                     state="disabled"
                                     )
        self.run_button.pack(side=BOTTOM)

        # ---- Section separator ----
        ttk.Label(right_frame, text="Console").pack()
        ttk.Separator(right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')


        # ---- Display output ----
        self.output = tk.Text(right_frame,
                              width=400,
                              height=5,
                              )
        self.output.pack()
        self.output.insert(END, "Output will be displayed here")

        self.status_label = ttk.Label(right_frame,
                                      text="Status: Idle",
                                      )
        self.status_label.pack()


        # Redirect console output to GUI
        # sys.stdout = TextRedirector(self.output, "stdout")
        # sys.stderr = TextRedirector(self.output, "stderr")

        # Quit button
        # self.quit = ttk.Button(left_frame, text="Quit RIANA",
        #                       command=self.master.destroy)
        # self.quit.pack(side="bottom")

        # ---- Section separator ----
        ttk.Label(right_frame, text="Results").pack()
        ttk.Separator(right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Using pandastable ----

        self.result_view = ttk.Frame(right_frame,
                                     width=400,
                                     height=300,
                                     )
        self.result_view.pack(fill=BOTH, expand=1)
        df = TableModel.getSampleData()
        self.result_table = Table(self.result_view,
                                  dataframe=df,
                                  showtoolbar=False,
                                  showstatusbar=False,
                                  )

        # self.result_table.setTheme('default')
        self.result_table.editable = False

        # self.result_table.cols.colheader.bgcolor = '#ECECEC'

        self.result_table.show()
        options = {'font': 'SF Pro Text',
                   'fontsize': 13,
                   'rowheight': 20,
                   'grid_color': '#FFFFFF',
                   'rowselectedcolor': '#B7D6F0',
                   'colselectedcolor': '#B7D6F0',
                   'boxoutlinecolor': '#FFFFFF',
                   'cellbackgr': '#FFFFFF',
                   }
        config.apply_options(options, self.result_table)

        # Set header attributes after show()
        self.result_table.colheader.bgcolor = '#ECECEC'
        self.result_table.colheader.textcolor = 'black'
        self.result_table.rowheader.bgcolor = '#ECECEC'
        self.result_table.rowheader.textcolor = 'black'


        # ---- Using tree view
        self.result_tree = ttk.Treeview(right_frame,
                                        columns = df.columns,
                                        show='headings',
                                        height=300,
                                        )

        df_col = df.columns

        # all the column name are generated dynamically.
        self.result_tree["columns"] = df_col
        counter = len(df)

        # generating for loop to create columns and give heading to them through df_col var.
        for x in range(len(df_col)):
            self.result_tree.column(x, width=100)
            self.result_tree.heading(x, text=df_col[x])
            # generating for loop to print values of dataframe in treeview column.
        for i in range(counter):
            self.result_tree.insert('', 'end', values=df.iloc[i, :].tolist())


        self.tree_scroll = ttk.Scrollbar(right_frame,

                                            orient="vertical",
                                            command=self.result_tree.yview)

        self.result_tree.configure(yscrollcommand=self.tree_scroll.set)
        self.tree_scroll.pack(side="right", fill="y")
        self.result_tree.pack()



    # percolator file dialog
    def percolator_file_dialog(self):

        self.id_path = filedialog.askopenfilename(initialdir="/",
                                             title="Select A File")
        # filetypes = #(("jpeg files","*.jpg"), ("all files","*.*")) )
        self.output.insert(END, self.id_path)
        print(self.id_path)

        # enable button
        if self.id_path is not None and self.mzml_path is not None:
            self.run_button["state"] = "normal"

    # mzml folder dialog
    def browse_mzml(self):
        self.mzml_path = filedialog.askdirectory()
        # self.label = ttk.Label(self)
        # self.label.pack()
        # self.label.configure(text = self.mzml)
        self.output.insert(END, self.mzml_path)

        # enable button
        if self.id_path is not None and self.mzml_path is not None:
            self.run_button["state"] = "normal"

    # dummy progress bar
    def progress(self):
        self.prog_bar = ttk.Progressbar(
            self.master, orient="horizontal",
            length=400, mode="indeterminate"
        )
        self.prog_bar.pack(side=BOTTOM)

    # start pytmt
    def run_riana_integrate(self):
        self.progress()
        self.prog_bar.start()
        self.queue = queue.Queue()

        integration_vars = IntegrationVars(
            id_path=self.id_path,
            mzml_path=self.mzml_path,
            out=self.out,
            unique=self.unique.get(),
            q_value=self.q_value.get(),
            r_time=self.r_time.get(),
            mass_tol=self.mass_tol.get(),
            mass_defect=self.mass_defect.get(),
            thread=self.thread,
            write_intensities=self.write_intensities,
            sample=self.sample.get(),
            iso=self.iso.get(),
        )

        IntegrationThreadedTask(self.queue, integration_vars).start()

        self.master.after(100, self.process_queue)
        # disable button
        self.run_button["text"] = "Running..."
        self.run_button["state"] = "disabled"
        self.status_label["text"] = "Status: Running..."

    def process_queue(self):
        try:
            msg = self.queue.get_nowait()
            # Show result of the task if needed
            self.prog_bar.stop()
            self.output.insert(END, msg)
            # with open(os.path.join(self.out, 'tmt.log'), 'r') as f:
            #     self.output.insert(INSERT, f.read())

            self.queue.task_done()
            self.run_button["text"] = "Integrate"
            self.run_button["state"] = "normal"
            # self.status_label["text"] = "Status: Finished"

        except queue.Empty:
            self.master.after(100, self.process_queue)
            # re-enable button


class IntegrationThreadedTask(threading.Thread):
    def __init__(self, queue, args):
        super().__init__()
        self.queue = queue
        self.args = args

    def run(self):
        riana_integrate.integrate_all(self.args)  # Simulate long running process
        self.queue.put("Task finished")

        #
        # read_res = pd.read_csv(os.path.join(self.out, self.sample.get()) + '_riana.txt')
        # print(read_res.head())

