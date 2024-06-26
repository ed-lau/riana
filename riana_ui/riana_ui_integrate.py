# Path: riana/riana_ui.py
# -*- coding: utf-8 -*-
# RIANA GUI

from typing import NamedTuple
import os
import sys
import ast
import tkinter as tk
import tqdm
from tkinter import ttk, filedialog, FLAT, BOTH, LEFT, TOP, END, BOTTOM, X
from riana import riana_integrate
import pandas as pd
import numpy as np
from pandastable import Table, TableModel, config
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from riana import constants
from scipy.signal import find_peaks
import rx

from console import TextRedirector

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
    iso: list[int]
    smoothing: int
    gui: bool



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
        self.output_file = 'out'

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
        self.mass_defect.set('C13')
        self.mass_defect_options: list[str] = ['C13', 'D', 'SILAC']
        # self.mass_defect_display: list[str] = [f'C13: {constants.C13_MASSDIFF}',
        #                                        f'Deuterium: {constants.D_MASSDIFF}',
        #                                        f'SILAC K/R: {constants.SILAC_MASSDIFF}']


        self.iso = tk.StringVar()
        self.iso.set('0,1,2,3,4,5')

        self.smoothing = tk.IntVar()
        self.smoothing.set(0)
        self.smoothing_options: list[int] = [0, 3, 5, 7, 9]
        # self.smoothing_display: list[str] = ['None', '3', '5', '7', '9']

        self.create_tab1_widgets()

        # For plotting intensities
        # self.intensities_df_subset = None
        self.canvas = None
        self.toolbar = None
        self.trace_plot = None
        self.trace_plot_2 = None
        self.bar_plot = None
        self.fig = None

        # sys.stdout = TextRedirector(self.master.master.console, "stdout")
        # sys.stderr = TextRedirector(self.master.master.console, "stderr")

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
        # Check the value is a list of numbers separated by spaces
        # try:
        #     self.iso.set([int(_) for _ in self.iso.get().split(',')])
        # except ValueError:
        #     self.iso.set([0])


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
        left_frame.pack(side=LEFT, fill=BOTH, expand=False)

        right_frame = ttk.LabelFrame(self,
                                     width=650,
                                     height=400,
                                     relief=FLAT,
                                     #borderwidth=1,
                                     text="Output",
                                     )
        right_frame.pack(side="right", fill=BOTH, expand=True)

        # ---- Section separator ----
        ttk.Label(left_frame, text="Select I/O Paths").pack()
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
                                      command=self.mzml_dialog,
                                      width=20,
                                      state="normal"
                                      )
        self.select_mzml.pack(side="top")

        # ---- Select output folder ----
        self.select_output = ttk.Button(left_frame,
                                        text="Select Output Folder",
                                        command=self.output_folder_dialog,
                                        width=20,
                                        state="normal"
                                        )
        self.select_output.pack(side="top")


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
                                        to=1,
                                        orient=tk.HORIZONTAL,
                                        variable=self.q_value,
                                        command=self.update_q_value,
                                        length=180,
                                        )
        self.select_q_value.set(0.1)
        self.select_q_value.pack(side="top")

        # ---- Select retention time threshold
        self.r_time_label = ttk.Label(left_frame,
                                       text=f'Retention time range: {round(self.r_time.get(), 3)}',
                                       )
        self.r_time_label.pack()
        self.select_r_time = ttk.Scale(left_frame,
                                        from_=0.1,
                                        to=3.0,
                                        orient=tk.HORIZONTAL,
                                        variable=self.r_time,
                                        command=self.update_r_time,
                                        length=180,
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

        # ---- Select smoothing ----
        self.smoothing_label = ttk.Label(left_frame,
                                        text=f'Smoothing:',
                                        )
        self.smoothing_label.pack()
        self.select_smoothing = ttk.OptionMenu(left_frame,
                                                    self.smoothing,
                                                    self.smoothing_options[0],
                                                    *self.smoothing_options,
                                                    )
        self.select_smoothing.pack(side="top")


        # link to RIANA integrate
        self.integrate_run_button = ttk.Button(left_frame,
                                     text="Run RIANA Integrate",
                                     command=self.handle_integrate_button,
                                     width=20,
                                     state="disabled"
                                     )
        self.integrate_run_button.pack(side=BOTTOM)

        self.cancel_button = ttk.Button(left_frame,
                                        text="Cancel",
                                        command=self.handle_cancel_button,
                                        width=20,
                                        state="disabled"
                                        )
        self.cancel_button.pack(side=BOTTOM)

        # sys.stdout = TextRedirector(self.master.master.console, "stdout")
        # sys.stderr = TextRedirector(self.master.master.console, "stderr")




        # ---- Section separator ----
        ttk.Label(right_frame, text="Results").pack()
        ttk.Separator(right_frame, orient='horizontal').pack(fill='x', pady=5, padx=5, anchor='w')

        # ---- Using pandastable ----

        self.integration_result_view = ttk.Frame(right_frame,
                                     width=600,
                                     height=300,
                                     )
        self.integration_result_view.pack(fill=BOTH, expand=1)
        # self.load_result_table()
        self.integration_result_view.pack_forget()

        # ---- Inspect view from selected row on the data ----
        self.inspect_view = ttk.Frame(right_frame,
                                     width=600,
                                     height=200,
                                     )
        # self.inspect_view.pack(fill=BOTH, expand=1)

        self.inspect_text = ttk.Label(self.inspect_view, text='')
        self.inspect_text.pack()





        # ---- Using tree view ----

        # self.result_tree = ttk.Treeview(right_frame,
        #                                 columns = df.columns,
        #                                 show='headings',
        #                                 height=300,
        #                                 )
        #
        # df_col = df.columns
        #
        # # all the column name are generated dynamically.
        # self.result_tree["columns"] = df_col
        # counter = len(df)
        #
        # # generating for loop to create columns and give heading to them through df_col var.
        # for x in range(len(df_col)):
        #     self.result_tree.column(x, width=100)
        #     self.result_tree.heading(x, text=df_col[x])
        #     # generating for loop to print values of dataframe in treeview column.
        # for i in range(counter):
        #     self.result_tree.insert('', 'end', values=df.iloc[i, :].tolist())
        #
        # self.tree_scroll = ttk.Scrollbar(right_frame,
        #
        #                                     orient="vertical",
        #                                     command=self.result_tree.yview)
        #
        # self.result_tree.configure(yscrollcommand=self.tree_scroll.set)
        # self.tree_scroll.pack(side="right", fill="y")
        # self.result_tree.pack()


    # ---- Functions ----
    def load_result_table(self):

        try:  # try to read in data from RIANA integrate
            output_file_path = os.path.join(self.output_file, self.sample.get() + '_riana.txt')
            # print(output_file_path)
            self.master.master.console.insert(END, f"Reading in data from {output_file_path}\n")
            df = pd.read_csv(output_file_path, sep='\t', index_col=0)

            # Open the intensities data frame so we can display the data traces
            intensity_file_path = os.path.join(self.output_file, self.sample.get() + '_riana_intensities_summarized.txt')
            intensities_df = pd.read_csv(intensity_file_path, sep='\t', index_col=0)


        except:  # if not, read in sample data
            df = TableModel.getSampleData()

        self.result_table = Table(self.integration_result_view,
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

        def handle_left_click(event):
            rowclicked_single = self.result_table.get_row_clicked(event)
            # print(rowclicked_single)

            # Clear the canvas

            try:
                self.canvas.get_tk_widget().pack_forget()
                for item in self.canvas.get_tk_widget().find_all():
                    self.canvas.get_tk_widget().delete(item)
                self.toolbar.pack_forget()
            except AttributeError:
                pass

            self.result_table.setSelectedRow(rowclicked_single)
            self.result_table.redraw()

            # Refresh the inspect view every time there is a new selection
            self.inspect_view.pack_forget()
            self.inspect_view.pack(fill=BOTH, expand=1)
            self.trace_plot = None
            self.trace_plot_2 = None

            # Get the selected row data
            selected_data = self.result_table.getSelectedRowData()
            # print(selected_data['pep_id'])

            # Filter the intensities data by the selected pep_id
            intensities_df_subset = intensities_df[intensities_df['pep_id'] == selected_data['pep_id'].values[0]]
            # print(intensities_df_subset['rt'].values[0])
            # print(intensities_df_subset['m0'].values[0])


            # Display the text on the interface
            out_text= f'Row: {self.result_table.getSelectedRow()}, ' \
                      f'pep_id: {selected_data["pep_id"].values[0]}, ' \
                      f'Peptide: {selected_data["sequence"].values[0]}+{selected_data["charge"].values[0]}, ' \
                      f'Protein: {selected_data["protein id"].values[0]} ' \
                      # f'Retention time: {intensities_df_subset["rt"].values[0]}, ' \
                      #  f'Mass: {intensities_df_subset["m0"].values[0]}'
            self.inspect_text.config(text=out_text)
            self.inspect_text.pack()

            # Display the graphs on the interface
            self.fig = Figure(figsize=(5, 3), dpi=100)
            self.fig.suptitle(f'{selected_data["sequence"].values[0]}+{selected_data["charge"].values[0]}')

            # Produce the chromatographic trace
            self.trace_plot = self.fig.add_subplot(131)
            self.trace_plot.set_xlabel('Retention time (min)')
            self.trace_plot.set_ylabel('Intensity')
            self.trace_plot.set_title(f'Chromatogram')
            self.trace_plot.grid()

            # Produce the chromatographic trace for peak finding
            self.trace_plot_2 = self.fig.add_subplot(132)
            self.trace_plot_2.set_xlabel('Retention time (min)')
            self.trace_plot_2.set_ylabel('Intensity')
            self.trace_plot_2.set_title(f'Chromatogram/Peak Finding')
            self.trace_plot_2.grid()

            # Turn the literal string into a list of floats
            rt = ast.literal_eval(intensities_df_subset['rt'].values[0])
            rt = [float(i) for i in rt]
            print(rt)

            # Get every isotopomer from the column names of the intensities df for plotting
            isotopomers = []
            for col in intensities_df_subset.columns:
                if col.startswith('m'):
                    isotopomers.append(col)

            colors = plt.cm.rainbow(np.linspace(0, 1, len(isotopomers)))

            # Loop through isotopomers and create a line plot for each
            for i, isotop in enumerate(isotopomers):
                m_trace = ast.literal_eval(intensities_df_subset[isotop].values[0])
                m_trace = [float(i) for i in m_trace]
                print(m_trace)

                self.trace_plot.scatter(x=rt,
                                        y=m_trace,
                                        s=5,
                                        label=isotop,
                                        color=colors[i],
                                        marker='o',
                                        linestyle='-',)
                self.trace_plot.plot(rt, m_trace, color=colors[i], linewidth=0.5)

                # Peak finding
                peaks, properties = find_peaks(m_trace, prominence=0.1, width=1)
                print(peaks)
                self.trace_plot_2.scatter(x=rt,
                                        y=m_trace,
                                        s=5,
                                        label=isotop,
                                        color=colors[i],
                                        marker='o',
                                        linestyle='-', )
                self.trace_plot_2.plot(np.array(rt)[peaks], np.array(m_trace)[peaks], "x", color=colors[i], linewidth=0.5)
                self.trace_plot_2.vlines(x=np.array(rt)[peaks], ymin= np.array(m_trace)[peaks] - properties["prominences"],
                                            ymax= np.array(m_trace)[peaks], color=colors[i], linewidth=0.5)
                self.trace_plot_2.plot(rt, m_trace, color=colors[i], linewidth=0.5)

            self.trace_plot.legend(loc='upper right')
            self.trace_plot_2.legend(loc='upper right')


            # Produce the isotopomer profile
            self.bar_plot = self.fig.add_subplot(133)
            self.bar_plot.set_xlabel('Isotopomer')
            self.bar_plot.set_ylabel('Intensity')
            self.bar_plot.set_title(f'Isotope Envelope')
            # self.bar_plot.grid()

            for i, isotop in enumerate(isotopomers):
                self.bar_plot.bar(x=isotop,
                                  height=selected_data[isotop].values[0],
                                  label=isotop,
                                  width=0.1,
                                  color=colors[i])

            self.bar_plot.legend(loc='upper right')

            # Display the plots on the interface
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.inspect_view)
            self.canvas.draw()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.inspect_view)
            self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            self.toolbar.update()
            self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.result_table.bind("<Button-1>", handle_left_click)




    # percolator file dialog
    def percolator_file_dialog(self):

        self.id_path = filedialog.askopenfile(initialdir='./tests/data/sample1',
                                             title="Select The Percolator File")
        # filetypes = #(("jpeg files","*.jpg"), ("all files","*.*")) )
        self.master.master.console.insert(END, self.id_path.name)

        # Enable integrate button
        if self.id_path is not None and self.mzml_path is not None:
            self.integrate_run_button["state"] = "normal"


    # mzml folder dialog
    def mzml_dialog(self):
        """ Returns a selected directoryname. """
        self.mzml_path = filedialog.askdirectory(initialdir='./tests/data/sample1')
        # self.label = ttk.Label(self)
        # self.label.pack()
        # self.label.configure(text = self.mzml)
        self.master.master.console.insert(END, self.mzml_path)
        # print(self.mzml_path)

        # enable button
        if self.id_path is not None and self.mzml_path is not None:
            self.integrate_run_button["state"] = "normal"

    # output folder dialog
    def output_folder_dialog(self):
        """ Returns a selected directoryname. """
        self.output_file = filedialog.askdirectory(initialdir='./out/ui_test/',
                                                   title='Select output folder',
                                                   )
        self.master.master.console.insert(END, self.output_file)
        print(self.output_file)

    # dummy progress bar
    def progress(self):
        self.prog_bar = ttk.Progressbar(
            self.master, orient="horizontal",
            length=400, mode="indeterminate"
        )
        self.prog_bar.pack(side=BOTTOM, fill=X)

    def handle_integrate_button(self):
        """
        Handle the model button click, run integration
        :return:
        """
        rx.empty().subscribe(
            on_completed=self.on_click,
            scheduler=self.master.master.pool_scheduler
        )

    def handle_cancel_button(self):
        """
        Cancel integration
        :return:
        """
        pass
        # # self.master.master.pool_scheduler
        # self.integrate_run_button["text"] = "Run"
        # self.integrate_run_button["state"] = "normal"
        # self.master.tab(1, state='normal')
        # self.master.master.status_label["text"] = "Status: Idle"


    def on_click(self):
        """
        Handle the model button click, run integration
        :return:
        """

        print("Console - Clicked")
        self.master.master.console.insert(END, "Clicked")

        self.progress()
        self.prog_bar.start()

        # Clear the result view when the button is pressed
        self.integration_result_view.pack_forget()
        self.inspect_view.pack_forget()

        # Disable all buttones
        self.integrate_run_button["text"] = "Running..."
        self.integrate_run_button["state"] = "disabled"
        self.cancel_button["state"] = "normal"

        self.master.tab(1, state='disabled')
        # self.master.Frame2.leftmodel_run_button["text"] = "Running..."
        # self.master.Frame2.model_run_button["state"] = "disabled"

        self.master.master.status_label["text"] = "Status: Running..."

        integration_vars = IntegrationVars(
            id_path=self.id_path,
            mzml_path=self.mzml_path,
            out=self.output_file,
            unique=self.unique.get(),
            q_value=self.q_value.get(),
            r_time=self.r_time.get(),
            mass_tol=self.mass_tol.get(),
            mass_defect=self.mass_defect.get(),
            thread=1,
            write_intensities=True,
            sample=self.sample.get(),
            iso=[int(i) for i in (self.iso.get().split(','))],
            smoothing=self.smoothing.get(),
            gui=True,
        )

        print("Integration variables: ", integration_vars)

        riana_integrate.integrate_all(integration_vars)

        print('Done')

        # if you want a callback on the main thread:
        self.master.master.after(5, self.on_task_complete)

    def on_task_complete(self):
        """
        After the model run is complete, stop the progress bar and update the status
        :return:
        """
        print("Task complete")

        self.prog_bar.stop()
        self.prog_bar.destroy()

        # self.queue.task_done()
        self.integrate_run_button["text"] = "Integrate"
        self.integrate_run_button["state"] = "normal"
        self.cancel_button["state"] = "disabled"
        self.master.tab(1, state='normal')

        self.master.master.status_label["text"] = "Status: Finished"

        # refresh the file handle in the filedialog
        self.id_path = open(self.id_path.name, 'r')

        # Update the result view
        self.load_result_table()
        self.integration_result_view.pack(fill=BOTH, expand=1)




