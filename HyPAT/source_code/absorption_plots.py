""" Code for the Absorption Plots tab in HyPAT """
import matplotlib
matplotlib.use("TkAgg")  # backend of matplotlib. needs to be here
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
import numpy as np
import pandas as pd
from tkinter.filedialog import askdirectory, asksaveasfilename
from tkinter.messagebox import showerror
import matplotlib.pyplot as plt
import os
from .data_storage import Widgets, FormatLabel, LoadingScreen
from scipy.optimize import curve_fit
import mplcursors  # adds the hover-over feature to labels. This library has good documentation
import openpyxl
import platform  # allows for Mac vs Windows adaptions
# make certain warnings appear as errors, allowing them to be caught using a try/except clause
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)  # used during curve_fit


class AbsorptionPlots(tk.Frame):

    def __init__(self, parent, storage, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        # container for important variables
        self.storage = storage
        self.datafiles = {}  # Store Excel files
        self.inputs = {}
        self.header = []  # Container for titles of each column in your data file that will be used
        self.needed_cols = []  # Container that will hold the column number for each needed data set
        self.converters_dict = {}  # Container for functions for converting data into correct units

        # store widgets
        self.widgets = Widgets()
        self.add_text0 = self.widgets.add_text0
        self.add_text2 = self.widgets.add_text2
        self.add_entry = self.widgets.add_entry
        self.add_entry3 = self.widgets.add_entry3

        # to get 2 rows over 3, we need two more frames
        self.top_frame = tk.Frame(self, bd=10)
        self.top_frame.grid(row=0, column=0, sticky="nsew")
        self.bottom_frame = tk.Frame(self, bd=10)
        self.bottom_frame.grid(row=1, column=0, sticky="nsew")

        self.frame = tk.LabelFrame(self.top_frame, bd=8)
        self.frame.grid(row=0, column=0, sticky="nsew")

        # arrange frames in tab
        self.rowconfigure((0, 1), weight=1)
        self.columnconfigure((0, 1), weight=1)
        self.top_frame.rowconfigure(0, weight=1)
        self.top_frame.columnconfigure((0, 1), weight=1)
        self.bottom_frame.rowconfigure(0, weight=1)
        self.bottom_frame.columnconfigure((0, 1, 2), weight=1)

        self.directory = ""
        self.loading_data = False  # Keep track of whether data is currently being loaded by this tab
        self.refreshing = False  # If false, ask for a directory when select_file gets called
        self.file_type = ""
        self.time_type = 0  # Used to determine which converter function to use for the time instrument.
        self.days = -1  # How many days have passed during the experiment (used for datetime conversion to seconds)
        self.this_date_time = 0  # Current datetime (used for datetime conversion to seconds)
        self.init_time = 0  # Initial datetime in seconds (used for datetime conversion to seconds)
        self.extra_time = 0  # Extra time to be accounted for (used for datetime conversion to seconds)
        self.error_texts = ""  # Text variable for storing the uncertainty texts that may come up when files are loaded

        self.num_GasT0 = 1  # Number of TCs measuring GasT initial
        self.num_GasT = 1  # Number of TCs measuring GasT after the isolation valve has opened

        # When determining uncertainty, the program will pick the max of the calculated statistical uncertainty, the
        # constant uncertainty, and the proportional uncertainty. These are containers for the constant and
        # proportional uncertainties of each instrument
        self.GasT_cerr = {}  # degrees C, i.e., +/- 2.2 degrees C (constant uncertainty)
        self.GasT_perr = {}  # percentage, i.e., +/- 0.75% of total measurement (in Celsius) (proportional uncertainty)
        self.SampT_cerr = {}
        self.SampT_perr = {}
        self.Pres_cerr = {}  # Pa (constant uncertainty)
        self.Pres_perr = {}  # percentage, (proportional uncertainty)

        # variables for input and uncertainty of inputs
        self.sample_thickness = tk.DoubleVar(value=0.1)
        self.sample_thickness_err = tk.DoubleVar(value=round(self.sample_thickness.get()*0.05, 13))
        self.Vs_cm3 = tk.DoubleVar(value=2.9e-2)   # Volume of initial container
        self.Vs_cm3_err = tk.DoubleVar(value=round(self.Vs_cm3.get() * 0.05, 13))
        self.Vic_cm3 = tk.DoubleVar(value=379.0)  # Volume of initial container
        self.Vic_cm3_err = tk.DoubleVar(value=round(self.Vic_cm3.get() * 0.05, 13))
        self.Vsc_cm3 = tk.DoubleVar(value=66.4)  # cm^3 Volume of sample container
        self.Vsc_cm3_err = tk.DoubleVar(value=round(self.Vsc_cm3.get() * 0.05, 13))
        self.ms_g = tk.DoubleVar(value=1)  # g Mass of sample
        self.ms_g_err = tk.DoubleVar(value=round(self.ms_g.get() * 0.05, 13))
        # todo Uncomment self.rhos_gcm3 and self.rho_gcm3_err here and elsewhere, then connect them to the code so they
        #      do things and update when needed
        # self.rhos_gcm3 = tk.DoubleVar(value=round(self.ms_g.get() / self.Vs_cm3.get(), 13))  # g/cm^3 density of sample
        # self.rhos_gcm3_err = tk.DoubleVar(value=round((self.ms_g.get() / self.Vs_cm3.get()) *
        #                                               np.sqrt((self.ms_g_err.get() / self.ms_g.get()) ** 2
        #                                                       + (self.Vs_cm3_err.get() / self.Vs_cm3.get())** 2), 13))
        self.exp_type = tk.StringVar(value="Single")
        self.molar_mass = tk.DoubleVar(value=50.9440)  # g/mol default molar_mass
        # Variables for determining initial and final values
        self.e_tol = tk.DoubleVar(value=self.storage.tol.get())  # Tolerance for equilibrium
        self.e_t_del = tk.DoubleVar(value=0)  # Minimum time after t0 before looking for eq.
        self.i_range = {}  # Dict for holding the range used in finding initial values for each file
        self.e_range = {}  # Dict for holding the range used in finding eq and eq value for each file
        self.gen_i_range = tk.IntVar(value=self.storage.gen_dp_range.get())  # Default range to find initial vals over
        self.gen_e_range = tk.IntVar(value=self.storage.gen_dp_range.get())  # Default range to find eq and eq vals over

        # absorption variables
        self.p_num = {}  # number of pressures
        self.t0 = {}  # time when isolation valve opens (s)
        self.t_e = {}  # time when equilibrium pressure is achieved (s)
        self.Phi = {}  # permeability
        self.Phi_err = {}
        self.Tg0 = {}  # Temperature of gas at t0 (using an average over time and instruments) (K)
        self.Tg_e = {}  # Temperature of gas at t_e (using an average over time and instruments) (K)
        self.Tg_e_err = {}
        self.artists = []  # store points for Perm/Dif/Sol vs Temp graph
        self.Ts0 = {}  # Sample temperature at t0 (using an average over time and instruments) (K)
        self.Ts_e = {}  # Sample temperature at t_e (using an average over time and instruments) (K)
        self.Ts_e_err = {}
        self.pr0_avg = {}  # Pressure at t0 (using an average) (Pa)
        self.pr_e_avg = {}  # Pressure at t_e (using an average) (Pa)
        self.pr_e_err = {}
        self.ns_e = {}  # Total moles absorbed by sample
        self.ns_e_err = {}
        self.ns_t = {}  # Number of moles in sample as a function of time
        self.HM = {}  # Hydrogen to metal atoms ratio
        self.HM_err = {}

        # diffusivity variables
        self.D = {}
        self.D_err = {}
        self.D_time = {}  # time over which D is calculated
        self.lhs = {}  # for the diffusivity optimization comparison
        self.rhs = {}  # analytical solution (non-fit)
        self.rhs_cf = {}  # curve_fit

        # solubility variables
        self.Ks = {}
        self.Ks_err = {}

        # text variables for labels
        self.pr_e_label = tk.DoubleVar(value=0)  # Equilibrium Pressure
        self.pr_e_err_label = tk.DoubleVar(value=0)
        self.ns_e_label = tk.DoubleVar(value=0)    # Moles absorbed
        self.ns_e_err_label = tk.DoubleVar(value=0)
        self.eq_t_label = tk.DoubleVar(value=0)  # Time to equilibrium
        self.eq_t_err_label = tk.DoubleVar(value=0)
        self.Phi_label = tk.DoubleVar(value=0)  # Permeability
        self.Phi_err_label = tk.DoubleVar(value=0)
        self.D_label = tk.DoubleVar(value=0)  # Diffusivity
        self.D_err_label = tk.DoubleVar(value=0)
        self.K_label = tk.DoubleVar(value=0)  # Solubility
        self.K_err_label = tk.DoubleVar(value=0)
        self.p_num_label = tk.IntVar(value=0)  # Number of pressures in an isotherm

        # add frames and buttons to view important variable

        entry_row = 0
        self.add_entry3(self, self.frame, variable=self.inputs, key="ls", text="Sample Thickness: l",
                        subscript="s", tvar1=self.sample_thickness, tvar2=self.sample_thickness_err,
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        units="[mm]", row=entry_row)
        entry_row += 1

        self.add_entry3(self, self.frame, variable=self.inputs, key="ms_g", text="Sample Mass: m",
                        subscript="s", tvar1=self.ms_g, tvar2=self.ms_g_err,
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        units="[g]", row=entry_row)
        entry_row += 1

        self.add_entry3(self, self.frame, variable=self.inputs, key="Vs_cm3", text="Sample Volume: V",
                        subscript="s", tvar1=self.Vs_cm3, tvar2=self.Vs_cm3_err, units="[cm\u00B3]",  # "[cm^3]"
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        row=entry_row)
        entry_row += 1

        self.add_entry3(self, self.frame, variable=self.inputs, key="Vic_cm3", text="Initial Container Volume: V",
                        subscript="ic", tvar1=self.Vic_cm3, tvar2=self.Vic_cm3_err, units="[cm\u00B3]",  # "[cm^3]"
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        row=entry_row)
        entry_row += 1

        self.add_entry3(self, self.frame, variable=self.inputs, key="Vsc_cm3", text="Sample Container Volume: V",
                        subscript="tc", tvar1=self.Vsc_cm3, tvar2=self.Vsc_cm3_err, units="[cm\u00B3]",  # "[cm^3]"
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        row=entry_row)
        entry_row += 1

        label_row = entry_row
        # todo Uncomment self.rhos_gcm3 and self.rho_gcm3_err stuff here and elsewhere, then connect them to the code so
        #      they do things and update when needed
        # self.add_text2(self.frame, text="Sample Density: \u03C1", subscript="s", tvar1=self.rhos_gcm3,
        #                tvar2=self.rhos_gcm3_err, units="[g/cm\u00B3]",
        #                row=label_row)  # No "_label" in rho's name because rho doesn't change based on file
        # label_row += 1
        self.add_text2(self.frame, text="Equilibrium Pressure: P", subscript="e", tvar1=self.pr_e_label,
                       tvar2=self.pr_e_err_label, units="[Pa]", row=label_row)
        label_row += 1
        self.add_text2(self.frame, text="Moles Absorbed by Sample: n", subscript="se", tvar1=self.ns_e_label,
                       tvar2=self.ns_e_err_label, units="[mol]", row=label_row)
        label_row += 1
        self.add_text2(self.frame, text="Time to Equilibrium: t", subscript="e", tvar1=self.eq_t_label,
                       tvar2=self.eq_t_err_label, units="[s]", row=label_row)
        label_row += 1
        self.add_text2(self.frame, text="Permeability: \u03A6", subscript="", tvar1=self.Phi_label,
                       tvar2=self.Phi_err_label,
                       units="[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]",   # "[mol/msPa^0.5]"
                       row=label_row)
        label_row += 1
        self.add_text2(self.frame, text="Diffusivity: D", subscript="", tvar1=self.D_label, tvar2=self.D_err_label,
                       units="[m\u00b2 s\u207b\u00b9]", row=label_row)
        label_row += 1
        self.add_text2(self.frame, text="Solubility: K", subscript="s", tvar1=self.K_label, tvar2=self.K_err_label,
                       units="[mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]", row=label_row)
        label_row += 1

        menu_row = label_row
        self.add_text0(self.frame, text="Current File:", subscript="", row=menu_row)
        # menu to choose which file to display
        self.current_file = tk.StringVar(value='No files yet')
        self.filemenu = tk.OptionMenu(self.frame, self.current_file, self.current_file.get())
        self.filemenu.grid(row=menu_row, column=1, columnspan=4, sticky='ew')
        self.current_file.trace_add("write", self.generate_plots)

        # menu to choose which measurement (P, D, K) to display in bottom left graph
        self.current_variable = tk.StringVar(value='Solubility')
        self.add_text0(self.frame, text="Current Measurement:", subscript="", row=menu_row + 1)
        self.PDK_menu = tk.OptionMenu(self.frame, self.current_variable, 'Solubility', 'Diffusivity', 'Permeability')
        self.PDK_menu.grid(row=menu_row + 1, column=1, columnspan=4, sticky='ew')
        self.current_variable.trace_add("write", self.update_PDK_plot)

        sidecol = 5
        self.add_text0(self.frame, text="Experiment Type:", subscript="", row=0, column=sidecol, sticky='w')
        self.experiment_type_menu = tk.OptionMenu(self.frame, self.exp_type, *list(["Single", "Isotherm"]))
        self.experiment_type_menu.config(bg="yellow")
        self.experiment_type_menu.grid(row=1, column=sidecol, columnspan=4,
                                       sticky="ew")

        self.add_text0(self.frame, text="Sample Molar Mass:", subscript="", row=2, column=sidecol, sticky="w")
        mm_frame = tk.Frame(self.frame)  # Frame for the molar mass
        mm_frame.grid(row=3, column=sidecol, sticky="w")
        self.add_entry(self, mm_frame, variable=self.inputs, key="molar_mass", text="", subscript="",
                       units="[g/mol]", tvar=self.molar_mass, row=1,
                       command=lambda tvar, variable, key: self.storage.check_for_number(tvar, variable, key))
        mm_frame.columnconfigure((0, 2), weight=1)
        self.inputs['molar_mass'].config(width=10)

        self.npframe = tk.Frame(self.frame)
        self.npframe.grid(row=4, column=sidecol, sticky="w")
        self.add_text0(self.npframe, text="Number of Pressures:", subscript="", row=0, sticky="w")
        FormatLabel(self.npframe, textvariable=self.p_num_label, borderwidth=1, relief="ridge",
                    ).grid(row=0, column=1, sticky="w")
        tk.Label(self.npframe, text="  ").grid(row=0, column=2, sticky='w')  # Give the label some right-side padding
        self.npframe.columnconfigure(0, weight=1)

        button_row = 7
        ttk.Style(self).configure('Tight.TButton', width="")  # Remove the extra space from the button

        self.b0 = ttk.Button(self.frame, text='Choose folder', command=self.select_file)
        self.b0.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1
        b1 = ttk.Button(self.frame, text='Refresh', command=self.refresh_graphs, style="Tight.TButton")
        b1.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1
        settings_b = ttk.Button(self.frame, text='Settings', command=self.adjust_persistent_vars, style="Tight.TButton")
        settings_b.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1
        b2 = ttk.Button(self.frame, text='Close popout plots', command=lambda: plt.close('all'))
        b2.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1
        b3 = ttk.Button(self.frame, text='Export to Excel', command=self.export_data)
        b3.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1
        b4 = ttk.Button(self.frame, text='Save current figures', command=self.save_figures)
        b4.grid(row=button_row, column=sidecol, sticky="w")
        button_row += 1

        # make all the rows in the top left frame have the same size and change size at the same rate
        self.frame.rowconfigure(tuple(range(menu_row + 2)), weight=1, uniform="row")

        # Configure how the columns squish so that only the text box columns squish
        self.frame.columnconfigure((0, 2, 4, 5), weight=1)
        # self.frame.columnconfigure(tuple(range(6)), weight=1)  # option for squishing all columns instead

        # create bottom left plot
        self.ax_title = "Solubility vs Temperature"
        self.ax_xlabel = "Temperature (\u00B0C)"
        self.ax_ylabel = "Solubility (mol m$^{-3}$ Pa$^{-0.5}$)"
        self.fig, self.ax, self.canvas, self.toolbar = self.add_plot(self.bottom_frame,
                                                                     xlabel=self.ax_xlabel,
                                                                     ylabel=self.ax_ylabel,
                                                                     title=self.ax_title,
                                                                     row=0, column=0, axes=[.3, .15, .65, .75])
        # create top right plot
        self.ax1_title = "Raw Data"
        self.ax1_xlabel = "Time (s)"
        self.ax1_ylabel = "Pressure (Pa)"
        self.ax12_ylabel = "Temperature (\u00B0C)"
        self.fig1, self.ax1, self.canvas1, self.toolbar1 = self.add_plot(self.top_frame,
                                                                         xlabel=self.ax1_xlabel,
                                                                         ylabel=self.ax1_ylabel,
                                                                         title=self.ax1_title,
                                                                         row=0, column=1, rowspan=1)
        self.ax12 = self.ax1.twinx()
        self.ax12.set_ylabel(self.ax12_ylabel)
        # The top right plot would flicker every time the order of magnitude of the cursor's location would change, so
        # the following line changes the format of the displayed x & y coordinates. The specific format chosen was the
        # result of trial and error in conjunction with editing the text of the entry boxes for e_tol and e_t_del
        # (which are now accessed by a button)
        self.ax12.format_coord = lambda x, y: "({:3g}, ".format(y) + "{:3g})".format(x)

        # Create bottom middle plot
        self.ax2_title = "Pressure-Composition-Temperature"
        self.ax2_xlabel = "Composition (H/M)"
        self.ax2_ylabel = "Pressure (Pa)"
        self.fig2, self.ax2, self.canvas2, self.toolbar2 = self.add_plot(self.bottom_frame,
                                                                         xlabel=self.ax2_xlabel,
                                                                         ylabel=self.ax2_ylabel,
                                                                         title=self.ax2_title,
                                                                         row=0, column=1, axes=[.15, .15, .7, .75])

        # Create bottom right plot
        self.ax3_title = "Diffusivity Optimization Comparison"
        self.ax3_xlabel = "Time (s)"
        self.ax3_ylabel = r"(n$_t$)/(n$_{inf}$)"
        self.fig3, self.ax3, self.canvas3, self.toolbar3 = self.add_plot(self.bottom_frame,
                                                                         xlabel=self.ax3_xlabel,
                                                                         ylabel=self.ax3_ylabel,
                                                                         title=self.ax3_title,
                                                                         row=0, column=2, axes=[.15, .15, .75, .75])

        # Counteracts a bug of constrained_layout=True which causes four plots to be generated but not shown until a
        # Popout Plot button is clicked:
        plt.close("all")

    def update_function(self, tvar, var_type, key):
        """ Checks if the user entry is a number (float or int). If not, revert the entered string to its prior state.
            Wanted to have this functionality in the widgets class, but it wasn't working.
            :param tvar:
            :param var_type:
            :param key:
            :return: """
        # Check to make sure the entry is a float
        try:
            svar = float(self.inputs[key + var_type].get())  # string var for formatting
        except ValueError:
            # i.e. user typed a word or left it blank
            tk.messagebox.showwarning("Invalid Entry", "Please enter a number.")
            svar = tvar.get()  # reset to what it was before the invalid entry was typed

        # update the entry to be nicely formatted and update the variable to be what was entered (or reset the variable)
        if key in self.inputs.keys() and key+"_err" in self.inputs.keys():
            if var_type == "var":
                self.inputs[key].delete(0, "end")
                self.inputs[key].insert(0, "{:.2e}".format(svar))
            else:
                self.inputs[key+"_err"].delete(0, "end")
                self.inputs[key+"_err"].insert(0, "{:.2e}".format(svar))
            tvar.set(svar)
        return True

    def add_plot(self, parent, xlabel='', ylabel='', title='', row=0, column=0, rowspan=1, axes=[.1, .15, .8, .75]):
        """ Create a plot according to variables passed in (parent frame, x-label, y-label, title, row, and column).
            axes are kept in case constrained_layout has problems. """
        # location of main plot
        frame = tk.Frame(parent)
        frame.grid(row=row, column=column, rowspan=rowspan, sticky="nsew")
        """ Note, the below line is to allow plots to behave better when not at optimal size. See details here:
        https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
        The above website says the following as of 10/16/2021: 
        "Currently Constrained Layout is experimental. The behaviour and API are subject to change, 
        or the whole functionality may be removed without a deprecation period..."
        """
        # todo Update this to the official layout feature. See:
        #      https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.5.0.html
        fig, ax = matplotlib.pyplot.subplots(constrained_layout=True)
        # If something goes wrong with constrained_layout, comment out the above line and uncomment the two lines below
        # fig = Figure()
        # ax = fig.add_axes(axes)  # [left, bottom, width, height]

        # Configure plots
        ax.tick_params(direction='in')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, frame)
        # add button to pop out the plot
        b = tk.Button(master=toolbar, text="Popout Plot", command=lambda x=ax: self.popout_plot(x))
        b.pack(side="left", padx=1)

        # Add a button that allows user input of various equilibrium variables
        if title == 'Raw Data':
            entry_frame = tk.Frame(toolbar)
            entry_frame.pack(side="left", padx=1)
            b = tk.Button(entry_frame, text='Equilibrium Variables', command=self.adjust_e_vars)
            b.grid(row=0, column=4, sticky="w")

        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        return fig, ax, canvas, toolbar

    def refresh_graphs(self, loading_warning=True):
        """ Sets self.refreshing = True so that, when select_file is called, HyPAT doesn't ask for a directory.
            Calls select_file, which reloads all the data using filenames previously uploaded and plots them """
        if self.loading_data and loading_warning:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return
        self.refreshing = True
        self.select_file()

    def select_file(self):
        """ Facilitates loading data and plotting the absorption tab plots """
        if self.loading_data and not self.refreshing:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return

        # If select_file was called because the user is choosing a new folder (if self.refreshing=False),
        # then ask them for a directory from which to get data. If not, set self.refreshing back to false and
        # continue with the previously loaded data
        if not self.refreshing:
            # Save to a temp variable in case user cancels loading a new folder
            temp_dir = askdirectory(initialdir=os.path.dirname(__file__))
            if temp_dir:
                self.directory = temp_dir
        else:
            self.refreshing = False
            temp_dir = True  # temp_dir is truthy unless the user cancels loading the data

        if self.directory and temp_dir:
            self.options = [file for file in os.listdir(self.directory) if (os.path.splitext(file)[1] == '.xls' or
                                                                            os.path.splitext(file)[1] == '.xlsx')]
            if not self.options:
                tk.messagebox.showwarning(title="Missing Files",
                                          message="Please select a folder with data in the XLS or XLSX format.")
                return

            self.loading_data = True
            try:
                self.get_persistent_vars()  # Loads variables from a data file for use in analysis
            except Exception as e:  # If there is a problem while trying to load the persistent variables
                # If you want to see the traceback as part of the error text, change the below value to true
                want_traceback = False
                if want_traceback:
                    import traceback
                    full_traceback = "\n" + traceback.format_exc()
                else:
                    full_traceback = ""

                showerror('Loading Error', 'Unknown Error while obtaining persistent variables. No files were' +
                                           ' processed. The following exception was raised: "' + str(e) + '".' +
                                           full_traceback + '\n\n')
                self.loading_data = False
                return

            self.error_texts = ""  # Reset error_texts each time a folder is loaded
            self.datafiles.clear()  # Clear the dictionary to avoid old data points appearing in PDK plot

            # Create the progress bar
            pb_len = len(self.options)
            pb = LoadingScreen(self)
            pb.add_progress_bar(pb_len)

            # loop through files, test to make sure they will work, then process them
            i = 0
            i2 = 0
            while i < len(self.options):  # While loop required so that files don't get skipped when a pop is needed
                filename = self.options[i]
                i += 1
                i2 += 1

                # check if the file is good
                status = self.check_file(os.path.join(self.directory, filename))

                pname = ""  # Text to be added to the filename when using isotherms as a delineator
                # Process the data by making the various calculations
                if status:
                    # Reset these variables for each loadable file
                    self.time_type = 0
                    self.days = -1
                    self.extra_time = 0

                    self.p_num[filename] = 1  # Default number of pressures for each file
                    pcount = 0  # Keep track of which pressure in an isotherm is being processed
                    old_pname = ""  # Store filename without the pressure number delineator
                    # Keep processing data until all the pressures in the isotherm are processed
                    while pcount < self.p_num[filename]:
                        try:  # Catch most other errors that can happen while processing the data
                            pcount += 1
                            # Keep track of pressure number delineators
                            if self.exp_type.get() == "Isotherm":
                                pname = " p#{}".format(pcount)
                                old_pname = " p#{}".format(pcount - 1)

                            # Load the data first time through; otherwise, truncate the old data
                            if pcount == 1:
                                self.datafiles[filename + pname] = \
                                    self.extract_data(os.path.join(self.directory, filename))
                                self.options.pop(self.options.index(filename))  # Remove old filename from the options
                            else:  # Takes only the data after t0
                                self.datafiles[filename + pname] = \
                                    (self.datafiles[filename + old_pname])[self.t0[filename + old_pname]:]
                                # Ensure the index is properly labeled
                                self.datafiles[filename + pname].reset_index(drop=True, inplace=True)
                                i += 1

                            self.options.insert(i - 1, filename + pname)  # Add in the properly named option

                            # Find row number where the Isolation Valve is 1st opened after being closed (aka, the start
                            # of the experiment). Done in a try/except clause in case no opening is detected
                            try:
                                _t0 = np.where((self.datafiles[filename + pname])['Isolation Valve'] == 0)[0]
                                closed = np.where(np.diff(_t0) > 1)[0]
                                if not closed.any():  # if the isolation valve wasn't closed again after being opened...
                                    # ...set t0 to be the data point after the last zero in the file
                                    self.t0[filename + pname] = _t0[-1] + 1
                                else:  # else, set t0 to the data point after the last zero before the valve is opened
                                    self.t0[filename + pname] = _t0[closed[0]] + 1

                                    # If the user is submitting an isotherm, and either there are data points after the
                                    # last data point where the valve is closed, or there more than one separate
                                    # places where the valve is closed after being opened (meaning there are multiple
                                    # places where the valve is open)
                                    if self.exp_type.get() == "Isotherm" and \
                                            (len((self.datafiles[filename + pname])['Isolation Valve']) > _t0[-1] + 1 or
                                             len(closed) > 1):
                                        self.p_num[filename] += 1
                            except IndexError:
                                self.error_texts += "Loading Error with file " + filename + pname + \
                                                    ". Incorrect Isolation Valve format.\n\n"
                                self.datafiles.pop(filename + pname)
                                status = False
                            if status:  # Ensure status hasn't changed since original definition
                                self.calculate_solubility(filename + pname, self.datafiles[filename + pname])
                                self.calculate_diffusivity(filename + pname, self.datafiles[filename + pname])
                                self.calculate_permeability(filename + pname, self.datafiles[filename + pname])
                        except Exception as e:
                            # If you want to see the traceback as part of the error text, change the below value to true
                            want_traceback = False
                            if want_traceback:
                                import traceback
                                full_traceback = "\n" + traceback.format_exc()
                            else:
                                full_traceback = ""

                            self.error_texts += "Unknown Error with file " + filename + pname + '. The following' + \
                                                ' exception was raised: "' + str(e) + '".' + full_traceback + '\n\n'
                            try:
                                self.datafiles.pop(filename + pname)
                            except KeyError:
                                self.error_texts += "Note: " + filename + pname + ' was never successfully loaded.' + \
                                                    '\n\n'
                            status = False
                            break
                if not status:  # If something is wrong with the file or with analyzing the file
                    self.options.pop(self.options.index(filename + pname))
                    i -= 1
                # Update the progress bar
                try:
                    pb.update_progress_bar(100 // pb_len, i2)
                except tk.TclError:  # if the user closed the progress bar window
                    # clear out remaining files
                    while i < len(self.options):
                        self.options.pop(i)
                    tk.messagebox.showwarning("Loading Canceled", "Loading has been canceled.")
                    self.error_texts += "Warning: Not all files in selected folder were analyzed.\n\n"
                    break
            pb.destroy()  # Close the progress bar

            self.update_dataframe()
            self.update_option_menu()

            try:
                self.current_file.set(self.options[0])
            except IndexError:
                self.error_texts += "No files were able to be processed.\n\n"
                self.current_file.set('No files yet')

            if self.error_texts:
                # Create a new popup window
                popup = tk.Tk()
                popup.wm_title("Error Messages")

                from tkinter import font  # Change the default font for text boxes to TkDefaultFont

                # Create a textbox with the error text in it
                errortext = tk.Text(popup, font=font.nametofont("TkFixedFont"),
                                    background=self.frame.cget("background"), padx=50, pady=50)
                errortext.insert("insert", self.error_texts)
                errortext.configure(state="disabled")  # Turn off editing the text
                errortext.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

                # Create a scrollbar for the text box
                err_scrollbar = tk.Scrollbar(popup, command=errortext.yview)
                err_scrollbar.grid(row=0, column=1, sticky="nsew")
                errortext.configure(wrap=tk.WORD, yscrollcommand=err_scrollbar.set)

        self.loading_data = False

    def check_file(self, filename):
        """ method to check Excel files for any problems, and to move on if there is an issue.
            :param filename:
            :return: """
        self.file_type = os.path.splitext(filename)[1]
        # Try to load data from the file. If an error is thrown, return false
        try:
            if self.file_type == ".xls":
                # This is faster for large files, but isn't working on .xlsx file formats
                data = pd.read_csv(filename, sep='\t', header=None)
            elif self.file_type == ".xlsx":
                # openpyxl supports .xlsx file formats. According to documentation, engine=None should
                # support xlsx and xls, but it doesn't seem to work for xls as of 11/10/2021
                data = pd.read_excel(filename, header=None, engine="openpyxl")
            else:  # This shouldn't ever get called if the program is working right
                self.error_texts += "Loading Error with file " + filename + ". This file type, " + self.file_type + \
                                    ", is unsupported.\n\n"
                return False
            return True
        except Exception as e:
            self.error_texts += "Loading Error with file " + filename + ". HyPAT was unable to read the file." + \
                                ' The following exception was raised: "' + str(e) + '".\n\n'
            return False

    def update_option_menu(self):
        """ update file selecting menu when data is loaded """
        menu = self.filemenu["menu"]
        menu.delete(0, "end")
        for string in self.options:
            menu.add_command(label=string, command=lambda value=string: self.current_file.set(value))
        self.b0.config(text='Choose new folder')

    def update_dataframe(self):
        """ Update the dataframe that is used for filling the Excel sheet and plotting the data in overview plots """ # . will loading data in the other tab overwrite this and casue problems?
        # this will recreate the dataframe each time
        df = pd.DataFrame()
        for filename in self.options:
            df = df.append(pd.DataFrame(
                {"Gas Temperature [K]": self.Tg_e[filename],
                 "Gas Temperature Uncertainty [K]": self.Tg_e_err[filename],
                 "Sample Temperature [K]": self.Ts_e[filename],
                 "Sample Temperature Uncertainty [K]": self.Ts_e_err[filename],
                 "Pressure [Pa]": self.pr_e_avg[filename],
                 "Pressure Uncertainty [Pa]": self.pr_e_err[filename],
                 "Permeability [mol m^-1 s^-1 Pa^-0.5]": self.Phi[filename],
                 "Permeability Uncertainty [mol m^-1 s^-1 Pa^-0.5]": self.Phi_err[filename],
                 "Diffusivity [m^2 s^-1]": self.D[filename],
                 "Diffusivity Uncertainty [m^2 s^-1]": self.D_err[filename],
                 "Solubility [mol m^-3 Pa^-0.5]": self.Ks[filename],
                 "Solubility Uncertainty [mol m^-3 Pa^-0.5]": self.Ks_err[filename],
                 "Sample Composition [H/M]": self.HM[filename],
                 "Sample Composition Uncertainty [H/M]": self.HM_err[filename]
                 }, index=[filename]
            ))
        self.storage.TransportParameters = df

    def export_data(self):
        """ Loads dataframe containing the information for each material into an Excel sheet """
        if self.loading_data:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return

        if self.current_file.get() != "No files yet":

            save_filename = asksaveasfilename(initialdir=os.path.dirname(__file__), initialfile="test",
                                              defaultextension="*.xlsx")
            if save_filename != '':
                self.storage.TransportParameters.to_excel(save_filename)

        else:
            showerror(message='Please load data first.')

    def l2n(self, col):
        """ Take Excel column names (such as J or AR) and return which column number they are (0-indexed)"""
        col_num = 0
        for i, l in enumerate(reversed(col)):
            col_num += (ord(l) - 64) * (26 ** i)
        return col_num - 1

    def unit_conversion(self, value, multi, add):  # todo This seems to do things cell by cell. Maybe find another way?
        """ Returns the linearly adjusted input according to multi and add """
        if value == "":
            return float("NaN")
        try:
            new_val = float(value) * multi + add
        except ValueError:
            self.error_texts += "Data Extraction Warning. The following was submitted to HyPAT for conversion" + \
                                ' into correct units: "' + str(value) + '". This triggered a ValueError. The' + \
                                " dataframe location corresponding to this term has been set to NaN. No further" \
                                " information is available.\n\n"
            return float("NaN")
        return new_val

    def get_seconds(self, date_time):
        """ Converts a datetime variable into seconds (excluding the date) """
        if self.days == -1:  # Equals -1 if this is the first term in the column being loaded. Resets with each file
            self.this_date_time = date_time
            self.init_time = date_time.hour * 3600 + date_time.minute * 60 + date_time.second + \
                             date_time.microsecond * 10 ** -6
        last_date_time = self.this_date_time
        self.this_date_time = date_time
        try:
            # Check to see if the day has changed between one cell and the next
            self.days = round(date_time.day - last_date_time.day, 5)  # round to avoid precision errors
            # If the day has changed, store the time that has passed in a new variable to be added to the new time
            if self.days != 0:
                self.extra_time += 86400
            value = date_time.hour * 3600 + date_time.minute * 60 + date_time.second + date_time.microsecond * 10 ** -6 + \
                self.extra_time - self.init_time
        except AttributeError:
            # This is expected if the cell is blank or otherwise not a datetime. NaNs sometimes trigger other errors
            # that HyPAT is more ready to handle than this AttributeError
            value = float("NaN")
        return value

    def convert_time(self, time_term, multi, add):
        """ Checks to see if the time-instrument's data is formatted in datetime or cumulative time passed and sets the
            converter function accordingly """
        if time_term == "":
            return float("NaN")
        if self.time_type == 0:  # Resets to zero before each file is loaded
            try:
                dummy = float(time_term)
                self.time_type = 1
            except TypeError:
                # TypeError is expected if time_term is a datetime
                self.time_type = 2
            except ValueError:
                self.error_texts += "Data Extraction Warning. The following was submitted to HyPAT for conversion" + \
                                ' into seconds: "' + str(time_term) + '". This triggered a ValueError.' + \
                                " The dataframe location corresponding to this term has been set to NaN." + \
                                " No further information is available.\n\n"
                return float("NaN")
        if self.time_type == 1:
            # If the time is in cumulative time passed format, do a linear transform according to multi and add
            return self.unit_conversion(time_term, multi, add)
        elif self.time_type == 2:
            # If the time is in datetime format, convert to cumulative seconds passed
            return self.get_seconds(time_term)

    def get_persistent_vars(self):
        """ Read data from persistent_absorption_input_variables.xlsx. This data is critical to processing the
            data files loaded in by the user. """
        # Open up the persistent variable file for reading
        pv_filename = os.path.join('data_files', 'persistent_absorption_input_variables.xlsx')
        pv_wb = openpyxl.load_workbook(pv_filename)

        self.num_GasT0 = pv_wb['Numbers']['C2'].value  # Number of TCs measuring GasT initial
        self.num_GasT = pv_wb['Numbers']['C3'].value  # Number of TCs measuring GasT after isolation valve is opened

        # Reset lists in case this is the second time they are loaded
        self.header = []
        self.needed_cols = []
        # Start the lists/dicts with time
        self.header.append('t')  # Generic name for data from this instrument
        self.needed_cols.append(self.l2n(pv_wb['MiscInfo']['A2'].value))  # Number corresponding to which column this instrument's data is in
        self.converters_dict[self.needed_cols[-1]] = \
            lambda input, multi=pv_wb['MiscInfo']['A3'].value, add=pv_wb['MiscInfo']['A4'].value: \
            self.convert_time(input, multi, add)  # Function to convert the instrument's data to correct units

        # Read data for each instrument that measures GasT into a dataframe, then save that info in convenient locations
        GasT_info = pd.read_excel(pv_filename, sheet_name="GasT", header=0)
        for GasT_inst in GasT_info.keys():  # Each column header corresponds with a TC measuring GasT
            self.header.append(GasT_inst)  # Generic name for the instrument
            self.needed_cols.append(self.l2n(GasT_info[GasT_inst][0]))  # Column number
            self.converters_dict[self.needed_cols[-1]] = \
                lambda input, multi=GasT_info[GasT_inst][1], add=GasT_info[GasT_inst][2]: \
                self.unit_conversion(input, multi, add)  # Function for unit conversion
            self.GasT_cerr[GasT_inst] = GasT_info[GasT_inst][3]  # constant uncertainty
            self.GasT_perr[GasT_inst] = GasT_info[GasT_inst][4]  # proportional uncertainty

        # Repeat the above for other variables
        # Sample temperature
        SampT_info = pd.read_excel(pv_filename, sheet_name="SampT", header=0)
        for SampT_inst in SampT_info.keys():
            self.header.append(SampT_inst)
            self.needed_cols.append(self.l2n(SampT_info[SampT_inst][0]))
            self.converters_dict[self.needed_cols[-1]] = \
                lambda input, multi=SampT_info[SampT_inst][1], add=SampT_info[SampT_inst][2]: \
                self.unit_conversion(input, multi, add)
            self.SampT_cerr[SampT_inst] = SampT_info[SampT_inst][3]
            self.SampT_perr[SampT_inst] = SampT_info[SampT_inst][4]

        # Pressure
        Pres_info = pd.read_excel(pv_filename, sheet_name="Pres", header=0)
        for Pres_inst in Pres_info.keys():
            self.header.append(Pres_inst)
            self.needed_cols.append(self.l2n(Pres_info[Pres_inst][0]))
            self.converters_dict[self.needed_cols[-1]] = \
                lambda input, multi=Pres_info[Pres_inst][1], add=Pres_info[Pres_inst][2]: \
                self.unit_conversion(input, multi, add)
            self.Pres_cerr[Pres_inst] = Pres_info[Pres_inst][3]
            self.Pres_perr[Pres_inst] = Pres_info[Pres_inst][4]

        # Isolation Valve shouldn't need more than one instrument, so don't need a for loop. Also don't need uncertainty
        self.header.append("Isolation Valve")
        self.needed_cols.append(self.l2n(pv_wb['MiscInfo']['B2'].value))
        self.converters_dict[self.needed_cols[-1]] = \
            lambda input, multi=pv_wb['MiscInfo']['B3'].value, add=pv_wb['MiscInfo']['B4'].value: \
            round(self.unit_conversion(input, multi, add))

        # Variable for which row in Excel sheet to start obtaining data from (0-indexed)
        self.starting_row = pv_wb['MiscInfo']['C2'].value
        # Variable for how many rows at the end of the Excel sheet to not read
        self.footer_rows = pv_wb['MiscInfo']['D2'].value

    def n_of_cols(self, filename):
        """ Given a .xls or .xlsx filename, return the number of columns in that file"""
        if self.file_type == ".xls":
            # This is faster for large files, but isn't working on .xlsx file formats
            df = pd.read_csv(filename, sep='\t', header=None, skiprows=self.starting_row, skipfooter=self.footer_rows)
        else:  # assume .xlsx file
            # openpyxl supports .xlsx file formats. According to documentation, engine=None should
            # support xlsx and xls, but it doesn't seem to work for xls as of 11/10/2021
            df = pd.read_excel(filename, header=None, engine="openpyxl", skiprows=self.starting_row,
                               skipfooter=self.footer_rows)
        return len(df.columns)

    def extract_data(self, filename):
        """ Extract data from given file """
        # Attempt to extract data from the file. If it fails because of a ValueError or IndexError, enter into a while
        # loop that can help the user fix the problem so the program can proceed
        try:
            if self.file_type == ".xls":
                # This is faster for large files, but isn't working on .xlsx file formats
                Data = pd.read_csv(filename, sep='\t', header=None, usecols=self.needed_cols,
                                   converters=self.converters_dict, skiprows=self.starting_row,
                                   skipfooter=self.footer_rows)
            else:  # assume .xlsx file
                # openpyxl supports .xlsx Excel file formats. According to documentation, engine=None should
                # support xlsx and xls, but it doesn't seem to work for xls as of 11/10/2021
                Data = pd.read_excel(filename, header=None, engine="openpyxl", usecols=self.needed_cols,
                                     converters=self.converters_dict, skiprows=self.starting_row,
                                     skipfooter=self.footer_rows)
            # https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
            new_header = [x for _, x in sorted(zip(self.needed_cols, self.header))]
            Data.columns = new_header
        except (ValueError, IndexError):
            # todo When there's a column name that references a column not in the data sheet, the following appears:
            #         FutureWarning: Defining usecols with out of bounds indices is deprecated and will raise a ParserError in a future version.
            #         return self._reader.parse(
            #      Unfortunately, ParserError doesn't seem to be recognized as an exception right now (late April /
            #      early May 2022), so this will need to be delt with later

            # Next line fixes some IndexErrors that occur when loading XLSX files if you give an incorrect column name
            # for the isolation valve or if you give a column name referencing a column outside the datasheet
            self.converters_dict.clear()
            new_pvars = False  # Whether new persistent variables (and self.converters_dict) have been obtained

            num_of_cols = self.n_of_cols(filename)
            while len(set(self.needed_cols)) < len(self.needed_cols) or \
                    not all(x < num_of_cols for x in self.needed_cols):
                # check to make sure all values in the list are less than the length of the Excel sheet
                # https://www.geeksforgeeks.org/python-check-if-all-the-values-in-a-list-are-less-than-a-given-value/
                if not all(x < num_of_cols for x in self.needed_cols):
                    tk.messagebox.showwarning("Data Reading Error",
                                              "Please ensure all instruments have column names that correspond to " +
                                              "columns in the Excel sheet.")
                    self.adjust_persistent_vars(retry=True)  # gives user a chance to fix the column names
                    new_pvars = True

                # check to make sure every column name is unique
                if len(set(self.needed_cols)) < len(self.needed_cols):
                    tk.messagebox.showwarning("Data Reading Error",
                                              "Please ensure all instruments have unique column names")
                    self.adjust_persistent_vars(retry=True)  # gives user a chance to fix the column names
                    new_pvars = True

            # Make sure self.converters_dict is filled since it was cleared above
            if not new_pvars:
                self.get_persistent_vars()

            # Now the errors are handled, extract the data from the Excel file and load it into a dataframe
            if self.file_type == ".xls":
                Data = pd.read_csv(filename, sep='\t', header=None, usecols=self.needed_cols,
                                   converters=self.converters_dict, skiprows=self.starting_row,
                                   skipfooter=self.footer_rows)
            else:  # assume .xlsx file
                Data = pd.read_excel(filename, header=None, engine="openpyxl", usecols=self.needed_cols,
                                     converters=self.converters_dict, skiprows=self.starting_row,
                                     skipfooter=self.footer_rows)
            # https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
            new_header = [x for _, x in sorted(zip(self.needed_cols, self.header))]
            Data.columns = new_header

        n = len(Data)
        # calculate numerical derivative of pressure
        deriv = np.zeros(n)
        for i in range(n - 1):
            deriv[i] = ((Data.loc[i + 1, 'Pres'] - Data.loc[i, 'Pres']) /
                        (Data.loc[i + 1, 't'] - Data.loc[i, 't']))
        Data['dPres'] = deriv.tolist()
        # rearrange last two columns so 'Pres' and 'dPres' are adjacent
        cols = Data.columns.tolist()
        cols = cols[:-2] + [cols[-1]] + [cols[-2]]
        Data = Data[cols]
        return Data

    def adjust_persistent_vars(self, retry=False):
        """ Function for calling the window for adjusting persistent variables for loading absorption data,
            then update everything if something was changed """
        loading_warn = True
        if self.loading_data and not retry:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return
            else:
                loading_warn = False

        self.update()
        # get location of gui
        x = self.winfo_rootx()
        y = self.winfo_rooty()
        # Call the "Adjust Persistent Variables" window and, if something was changed, refresh graphs and variables
        d = AdjustPVars(storage=self.storage, pos=(x, y)).show_pv_window()
        if d and not retry:
            self.refresh_graphs(loading_warning=loading_warn)  # if a warning has already been given don't give it again
        # retry = True if adjust_persistent_vars is called but only the persistent variables need to be updated
        # (as opposed to updating the graphs as well)
        if d and retry:
            self.get_persistent_vars()  # Loads data from the absorption input file-format file

    def adjust_e_vars(self):
        """ Function for calling the window for adjusting equilibrium variables for loading absorption data,
            then update everything if desired """
        popup = tk.Toplevel(self)
        popup.wm_title("Adjust Equilibrium Variables")

        # Place the popup window near the cursor
        pos_right = self.winfo_pointerx()
        pos_down = self.winfo_pointery()
        popup.geometry("+{}+{}".format(pos_right, pos_down))

        entry_frame = tk.Frame(popup)
        entry_frame.pack(side="left", padx=1)
        # Tolerance for determining equilibrium
        self.add_entry(popup, entry_frame, self.inputs, key="e_tol",
                       text="Tolerance for Equilibrium:", subscript='', units="Pa s\u207b\u00b2",
                       tvar=self.e_tol, formatting=True, ent_w=8, in_window=True,
                       command=lambda tvar, variable, key, formatting, pf:
                       self.storage.check_for_number(tvar, variable, key, formatting, pf))
        # Minimum time to let pass after t0 before checking for equilibrium
        self.add_entry(popup, entry_frame, self.inputs, key="e_t_del",
                       text="Delay until Equilibrium:", subscript='', units="s",
                       tvar=self.e_t_del, ent_w=8, row=1, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))  # "false" to say no formatting
        # Number of data points used in determining the initial values
        self.add_entry(popup, entry_frame, self.inputs, key="i_range",
                       text="Initial Values Range:", subscript='', units="data points",
                       tvar=self.gen_i_range, ent_w=8, row=3, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))
        # Number of data points used in determining the equilibrium and max time and final values
        self.add_entry(popup, entry_frame, self.inputs, key="e_range",
                       text="Equilibrium Range:", subscript='', units="data points",
                       tvar=self.gen_e_range, ent_w=8, row=4, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))

        button_row = 5
        b1 = ttk.Button(entry_frame, text='Close & Refresh', command=lambda: self.close_and_refresh(popup))
        b1.grid(row=button_row, column=0, sticky="ew")
        b2 = ttk.Button(entry_frame, text='Refresh', command=self.refresh_graphs)
        b2.grid(row=button_row, column=2, sticky="ew")

    def close_and_refresh(self, win):
        """ Accepts window argument, then closes that window after refreshing the graphs of the Absorption Plots tab """
        self.refresh_graphs()
        win.destroy()

    def calculate_solubility(self, filename, data):
        """ Calculate Solubility using loaded data """
        R = self.storage.R * 1000  # J/mol K
        Vs = self.Vs_cm3.get() * 10 ** (-6)  # m^3
        Vs_err = self.Vs_cm3_err.get() * 10 ** (-6)  # m^3
        V_ic = self.Vic_cm3.get() * 10 ** (-6)  # m^3
        V_ic_err = self.Vic_cm3_err.get() * 10 ** (-6)  # m^3
        V_sc = self.Vsc_cm3.get() * 10 ** (-6) - Vs  # m^3 Note the correction for the sample's volume
        V_sc_err = self.Vsc_cm3_err.get() * 10 ** (-6)  # m^3
        V_totc = V_ic + V_sc  # m^3, total container
        V_totc_err = np.sqrt(V_ic_err ** 2 + V_sc_err ** 2)  # m^3

        # Set the init range to the general init range (which is editable by the user) or to something that works better
        if self.gen_i_range.get() < self.t0[filename]:
            self.i_range[filename] = self.gen_i_range.get()
        else:
            self.i_range[filename] = self.t0[filename] - 1
            self.error_texts += "Initial Range Warning with file " + filename + ". The user-input initial range" + \
                                " exceeded the limit of " + str(self.t0[filename]) + " data points. Initial range" + \
                                " set to " + str(self.i_range[filename]) + " data points.\n\n"
        # Set the eq range to the general eq range (which is editable by the user) or to something that works better
        if self.gen_e_range.get() < len(data['Pres']) - self.t0[filename]:
            self.e_range[filename] = self.gen_e_range.get()
        else:
            self.e_range[filename] = len(data['Pres']) - self.t0[filename] - 2
            self.error_texts += "Equilibrium Warning with file " + filename + ". The user-input equilibrium" + \
                                " range exceeded the limit of " + str(len(data['Pres']) - self.t0[filename]) + \
                                " data points. Equilibrium range set to " + str(self.e_range[filename]) + \
                                " data points.\n\n"

        # determine row number when equilibrium pressure is achieved, checking to ensure it is before the end of
        # the file and after the minimum time delay + t0
        dPres = data['dPres'].rolling(window=self.e_range[filename], center=True, min_periods=1).mean()
        # List of all terms in dPres which are less than -5e-3 (allowing for a slight decline during equilibrium)
        neg_dPres = np.where(dPres < -5E-3)[0]
        # Minimum terms in a row required to determine if the pressure drop off was reached
        min_seq = max(self.e_range[filename] // 2 - 1, 1)
        # variable to store the last point before experiment ends (and/or pressure drops and/or isolation valve opens)
        e_time_max = len(dPres)

        # Check for valve being closed after it was opened and set the max data point accordingly
        IVclosed = np.where(((self.datafiles[filename])['Isolation Valve'])[self.t0[filename]:] == 0)[0]
        if IVclosed.any():
            e_time_max = (IVclosed[0] + self.t0[filename])

        past_drop = False
        #  Find where the pressure drops off because the valve was opened or fail to find such a point
        for count, value in enumerate(neg_dPres[:len(neg_dPres) - min_seq]):
            if e_time_max > value > (self.t0[filename]):  # if the value is in the range for looking for equilibrium,
                # ...check to see if we are past the initial drop in pressure
                if np.diff(neg_dPres[count:count + 2]) > 2 and not past_drop:
                    past_drop = True
                # ...check for a sequence after the initial drop that is min_seq long in which neg_dPres is less than 0
                if sum(np.diff(neg_dPres[count:count + min_seq])) <= min_seq and past_drop:
                    e_time_max = neg_dPres[count]  # set the max time to when the pressure drops
                    break

        e_del = np.where(data.loc[self.t0[filename]:, 't'] > self.e_t_del.get() + data.loc[self.t0[filename], 't']
                         + 0.000001)[0]  # Added .000001 to ensure precision errors don't lead to e_del being 0
        try:
            # Set the time delay to be the minimum number of data points away from t0 that still is further than
            # the user's input time delay
            e_del = e_del[0] - 1  # -1 to compensate for the > sign in e_del's def and in later checks
        except IndexError:  # This is expected if e_del is empty, generally because self.e_t_del is too long in some way
            # Set e_del to the largest delay that still works
            e_del = e_time_max - (self.e_range[filename] + self.t0[filename] + 2)
            self.error_texts += "Equilibrium Warning with file " + filename + ". The user-input time delay" + \
                                " exceeded the limit of " + \
                                str(round(data.loc[e_time_max - (self.e_range[filename] + 2), 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s. Time delay set to " + \
                                str(round(data.loc[e_del + self.t0[filename], 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s.\n\n"
        if e_del > e_time_max - (self.e_range[filename] + self.t0[filename] + 2):
            # In case e_del produces a usable data point but that data point is past an abrupt downturn in pressure or
            # is past the valve closing,
            e_del = e_time_max - (self.e_range[filename] + self.t0[filename] + 2)
            self.error_texts += "Equilibrium Warning with file " + filename + ". The user-input time delay" + \
                                " exceeded the limit of " + \
                                str(round(data.loc[e_time_max - (self.e_range[filename] + 2), 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s. Time delay set to " + \
                                str(round(data.loc[e_del + self.t0[filename], 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s.\n\n"
        # Find the data point at which an equilibrium is reached
        ddPres = pd.Series((dPres.loc[2:].to_numpy() - dPres.loc[:len(dPres) - 3].to_numpy()) /
                           (data.loc[2:, 't'].to_numpy() - data.loc[:len(dPres) - 3, 't'].to_numpy()))
        ddPres = ddPres.rolling(window=self.e_range[filename], center=True, min_periods=1).mean()
        zeros = np.where(abs(ddPres) < self.e_tol.get())[0]  # Find points where ddPres is approximately zero
        zeros = [z for z in zeros if self.t0[filename] + e_del < z < e_time_max - min_seq - 2]
        # If no equilibrium was found with the user-input tolerance, loop until find an equilibrium using new_e_tol
        new_e_tol = self.e_tol.get()
        while not zeros and new_e_tol < self.e_tol.get() * 100000:
            new_e_tol = float("{:.2e}".format(new_e_tol * 5))
            zeros = np.where(abs(ddPres) < new_e_tol)[0]
            zeros = [z for z in zeros if self.t0[filename] + e_del < z < e_time_max - min_seq - 2]
        if new_e_tol != self.e_tol.get() and zeros:
            self.error_texts += "Equilibrium Warning with file " + filename + ". HyPAT was unable to find an" + \
                                " equilibrium using the user-input tolerance " + str(self.e_tol.get()) + \
                                ". New tolerance was set to {:.2e}".format(new_e_tol) + \
                                ", which successfully determined a time for the beginning of the equilibrium.\n\n"
        elif new_e_tol != self.e_tol.get():
            zeros = [min(self.t0[filename] + e_del + 1, e_time_max - min_seq - 3)]  # Added a +1 because e_del can be zero now, so this forces t_e != t0
            self.error_texts += "Equilibrium Error with file " + filename + \
                                ". HyPAT was unable to find equilibrium using the user-input tolerance " + \
                                str(self.e_tol.get()) + " or the computationally set tolerance " + \
                                "{:.2e}".format(new_e_tol) + ". Equilibrium time was set to " + \
                                str(round(data.loc[zeros[0], 't'], 3)) + \
                                " s, which is equal to either the starting time plus the user-input" + \
                                " time delay until equilibrium or the maximum usable time, whichever is lower.\n\n"
        self.t_e[filename] = zeros[0]

        # Get the initial gas temperature using an average of arbitrarily many TCs
        Tg_ivals = pd.DataFrame()
        for i in range(1, self.num_GasT0 + 1):
            Tg_ivals['GasT' + str(i)] = data.loc[:, 'GasT' + str(i)]
        Tg_imean = Tg_ivals.mean(axis=1)
        self.Tg0[filename] = Tg_imean.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1].mean() + \
            self.storage.standard_temp  # K
        # Get the equilibrium gas temperature using an average of arbitrarily many TCs
        Tg_fvals = pd.DataFrame()
        for i in range(1, self.num_GasT + 1):
            Tg_fvals['GasT' + str(i)] = data.loc[:, 'GasT' + str(i)]
        Tg_fmean = Tg_fvals.mean(axis=1)
        self.Tg_e[filename] = Tg_fmean.loc[self.t_e[filename] + 1:self.t_e[filename] +
                                           self.e_range[filename]].mean() + self.storage.standard_temp  # K

        # calculate average pressures
        self.pr0_avg[filename] = \
            data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 'Pres'].mean()
        self.pr_e_avg[filename] = \
            data.loc[self.t_e[filename] + 1: self.t_e[filename] + self.e_range[filename], 'Pres'].mean()

        # calculate moles
        nc0 = self.pr0_avg[filename] * V_ic / R / self.Tg0[filename]
        nc_e = self.pr_e_avg[filename] * V_totc / R / self.Tg_e[filename]
        self.ns_e[filename] = nc0 - nc_e

        # Pressure as a function of time (Pa)
        pr_t = data.loc[self.t0[filename] + 1:self.t_e[filename] + self.e_range[filename], 'Pres'].to_numpy()
        # Moles in the sample as a function of time (mol)
        self.ns_t[filename] = nc0 - pr_t * V_totc / R / \
            (Tg_fmean.loc[self.t0[filename] + 1:self.t_e[filename] + self.e_range[filename]].to_numpy()
             + self.storage.standard_temp)

        # calculate solubility
        self.Ks[filename] = self.ns_e[filename] / Vs / np.sqrt(self.pr_e_avg[filename])

        # Calculate the uncertainty in solubility

        # Find uncertainty of Tg0 when there are arbitrarily many TCs
        Tg0_errs = []
        for i in range(1, self.num_GasT0 + 1):
            Tg0_errs.append((max(data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 'GasT' + str(i)].std(),
                                 self.GasT_cerr['GasT' + str(i)],
                                 (self.Tg0[filename] - self.storage.standard_temp) *
                                 self.GasT_perr['GasT' + str(i)] / 100)) ** 2)
        Tg0_err = np.sqrt(sum(Tg0_errs)) / self.num_GasT0
        # Find uncertainty of Tg_e when there are arbitrarily many TCs
        Tge_errs = []
        for i in range(1, self.num_GasT + 1):
            Tge_errs.append((max(data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 'GasT' + str(i)].std(),
                                 self.GasT_cerr['GasT' + str(i)],
                                 (self.Tg_e[filename] - self.storage.standard_temp) *
                                 self.GasT_perr['GasT' + str(i)] / 100)) ** 2)
        Tge_err = np.sqrt(sum(Tge_errs)) / self.num_GasT

        # Pressure uncertainties
        Pres0_err = max(data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 'Pres'].std(),
                        self.Pres_cerr["Pres"], self.pr0_avg[filename] * self.Pres_perr["Pres"] / 100)
        Pres_e_err = max(data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 'Pres'].std(),
                         self.Pres_cerr["Pres"], self.pr_e_avg[filename] * self.Pres_perr["Pres"] / 100)

        nc0_err = abs(nc0) * np.sqrt((Pres0_err / self.pr0_avg[filename]) ** 2 + (V_ic_err / V_ic) ** 2 +
                                     (Tg0_err / self.Tg0[filename]) ** 2)
        nce_err = abs(nc_e) * np.sqrt((Pres_e_err / self.pr_e_avg[filename]) ** 2 + (V_totc_err / V_totc) ** 2 +
                                      (Tge_err / self.Tg_e[filename]) ** 2)
        nse_err = np.sqrt(nc0_err ** 2 + nce_err ** 2)

        self.Ks_err[filename] = abs(self.Ks[filename]) * np.sqrt(
            (nse_err / self.ns_e[filename]) ** 2 + (Vs_err / Vs) ** 2 + (Pres_e_err / 2 / self.pr_e_avg[filename]) ** 2)

        # store uncertainties
        self.ns_e_err[filename] = nse_err
        self.Tg_e_err[filename] = Tge_err
        self.pr_e_err[filename] = Pres_e_err

        # Data for pressure vs composition graph
        NA = self.storage.Na  # number of atoms per mole, the Avogadro constant
        a_num = self.ms_g.get() * NA / self.molar_mass.get()  # number of atoms of metal
        H_num = 2 * self.ns_e[filename] * NA  # Number of Hydrogen atoms in the sample
        self.HM[filename] = H_num / a_num  # Hydrogen to metal atoms ratio
        self.HM_err[filename] = abs(self.HM[filename]) * np.sqrt((self.ms_g_err.get()/self.ms_g.get()) ** 2 +
                                                                 (nse_err / self.ns_e[filename]) ** 2)

        # get temp from sample TC
        self.Ts0[filename] = data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 'SampT'].mean() + \
            self.storage.standard_temp
        self.Ts_e[filename] = data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 'SampT'].mean() + \
            self.storage.standard_temp
        self.Ts_e_err[filename] = max(data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 'SampT'].std(),
                                      self.SampT_cerr["SampT"], self.Ts_e[filename] * self.SampT_perr["SampT"] / 100)

    def calculate_diffusivity(self, filename, data):
        """ Calculate Diffusivity using loaded data """
        debug = False

        # Prepare for analysis
        sl = self.sample_thickness.get() * 0.001  # converted to meters
        sl_err = self.sample_thickness_err.get() * 0.001  # converted to meters  # todo Find a way to use this
        self.D_time[filename] = data.loc[self.t0[filename] + 1:self.t_e[filename] + self.e_range[filename], 't'] - \
                                data.loc[self.t0[filename] + 1, 't']
        h = data.loc[1, 't'] - data.loc[0, 't']

        # Calculate D without fitting for an initial guess. Source: Crank, pg238-239
        halfway_point = np.argmin(abs(self.ns_t[filename] / self.ns_e[filename] - 0.5))
        halfway_t = self.D_time[filename][halfway_point + self.t0[filename]]  # +t0 because of indexing
        prop_const = -np.log(np.pi ** 2 / 16 - (1 / 9)*(np.pi ** 2 / 16) ** 9) / np.pi ** 2  # Proportionality constant
        D = prop_const * sl ** 2 / halfway_t  # Diffusivity
        self.rhs[filename] = 1 - sum([(8 / ((2 * n + 1) ** 2 * np.pi ** 2)) * np.exp(
            -D * (2 * n + 1) ** 2 * np.pi ** 2 * (self.D_time[filename]) /
            (4 * (sl / 2) ** 2)) for n in range(0, 20)])

        # Function for optimizing D using curve fit
        def f(xdata, D_opt, dt_opt, A_opt):
            rhs = 1 - sum([(8 / ((2 * n + 1) ** 2 * np.pi ** 2)) *
                           np.exp(-D_opt * (2 * n + 1) ** 2 * np.pi ** 2 * (xdata + dt_opt) / (4 * (sl / 2) ** 2))
                           for n in range(0, 20)])  # sl/2 because the book's equation goes from -l to l, not 0 to l
            return rhs / A_opt

        # Attempt to optimize D using curve fit and the above function. Show a warning if it fails
        self.lhs[filename] = self.ns_t[filename] / self.ns_e[filename]  # assumes ns0 in sample is zero
        try:
            popt, pcov = curve_fit(f, self.D_time[filename], self.lhs[filename], p0=[D, 0, 1], xtol=D*1e-3,
                                   bounds=([0, 0, -1000], [10, 10*h, 1000]))
        except RuntimeError:
            self.error_texts += "Curve Fit Error with file " + filename + \
                                ". Curve fit unable to find optimal diffusivity parameters.\n\n"
            NaN = float("NaN")
            popt = [D, 0, 1]
            pcov = np.array([[NaN, NaN, NaN], [NaN, NaN, NaN], [NaN, NaN, NaN]])
        except OptimizeWarning:
            self.error_texts += "Curve Fit Warning with file " + filename + \
                                ". Curve fit unable to find covariance of the diffusivity parameters.\n\n"
            NaN = float("NaN")
            popt = [D, 0, 1]
            pcov = np.array([[NaN, NaN, NaN], [NaN, NaN, NaN], [NaN, NaN, NaN]])
        # Get calculated uncertainty
        perr = np.sqrt(np.diag(pcov))
        self.D_err[filename] = perr[0]

        if debug:
            print('\n', filename)
            print("New D, dt, A:", popt)
            print('pcov:\n', pcov)
            print("perr:", perr)

        # Store D and rhs_cf for use elsewhere in program
        self.D[filename], dt, A = popt
        self.rhs_cf[filename] = 1 - sum([(8 / ((2 * n + 1) ** 2 * np.pi ** 2)) * np.exp(
            -self.D[filename] * (2 * n + 1) ** 2 * np.pi ** 2 * (self.D_time[filename] + dt) /
            (4 * (sl / 2) ** 2)) for n in range(0, 20)])

        if debug:
            plt.figure()
            plt.plot(self.lhs[filename][1:], label='lhs')
            plt.plot(np.array(self.rhs_cf[filename][1:]), label='rhs')
            plt.title(filename + " curve fit")
            plt.legend()
            plt.show()

    def calculate_permeability(self, filename, data):
        """ Calculate Permeability using calculated Solubility and Diffusivity """
        self.Phi[filename] = self.D[filename] * self.Ks[filename]
        # get the uncertainty
        self.Phi_err[filename] = self.Phi[filename] * np.sqrt(
            (self.Ks_err[filename] / self.Ks[filename]) ** 2 + (self.D_err[filename] / self.D[filename]) ** 2)

    def generate_plots(self, *args):
        """ Create the absorption tab plots """
        filename = self.current_file.get()
        if self.current_file.get() == "No files yet":
            return

        data = self.datafiles[filename]

        # clear plots
        self.ax.clear()
        self.ax1.clear()
        self.ax12.clear()
        self.ax2.clear()
        self.ax3.clear()

        # Permeability/Diffusivity/Solubility vs Temperature (bottom left graph)
        self.PDK_plot(self.ax)

        # Pressure vs Time (top right graph)
        self.pressure_time_plot(data, filename, self.fig1, self.ax1, self.ax12)

        # Pressure vs Composition (bottom middle graph) (color coded by temperature)
        self.pres_comp_temp_plot(self.fig2, self.ax2)

        # Diffusivity comparison plots for optimization (bottom right graph)
        self.comparison_plot(filename, self.ax3)

        self.canvas.draw()
        self.toolbar.update()
        self.canvas1.draw()
        self.toolbar1.update()
        self.canvas2.draw()
        self.toolbar2.update()
        self.canvas3.draw()
        self.toolbar3.update()

        # update labels
        self.pr_e_label.set(self.pr_e_avg[filename])
        self.ns_e_label.set(self.ns_e[filename])
        self.eq_t_label.set(data.loc[self.t_e[filename], 't'] - data.loc[self.t0[filename], 't'])
        self.Phi_label.set(self.Phi[filename])
        self.D_label.set(self.D[filename])
        self.K_label.set(self.Ks[filename])
        self.pr_e_err_label.set(self.pr_e_err[filename])
        self.ns_e_err_label.set(self.ns_e_err[filename])
        # todo The following error assumes 1 sigma of confidence that t0 and t_e are within one data point of the
        #      correct value for t0 or t_e, respectively. A better estimate may be needed.
        t0_s_err = (data.loc[self.t0[filename] + 1, 't'] - data.loc[self.t0[filename] - 1, 't']) / 2
        te_s_err = (data.loc[self.t_e[filename] + 1, 't'] - data.loc[self.t_e[filename] - 1, 't']) / 2
        self.eq_t_err_label.set(np.sqrt(t0_s_err ** 2 + te_s_err ** 2))
        self.Phi_err_label.set(self.Phi_err[filename])
        self.D_err_label.set(self.D_err[filename])
        self.K_err_label.set(self.Ks_err[filename])
        if self.exp_type.get() == "Isotherm":
            fname = (filename.rsplit(maxsplit=1))[0]  # Remove the " p#_" by removing info after a whitespace
        else:
            fname = filename
        self.p_num_label.set(self.p_num[fname])

    def update_PDK_plot(self, *args):
        """ Clear and replot the bottom left (PDK) plot """
        if self.current_file.get() == "No files yet":
            return
        self.ax.clear()
        self.PDK_plot(self.ax)
        self.canvas.draw()
        self.toolbar.update()

    def popout_plot(self, ax):
        """ creates the selected plot using the regular popout from matplotlib """

        if self.loading_data:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return

        filename = self.current_file.get()
        if filename != "No files yet":
            data = self.datafiles[filename]
            plot = ax.get_title()
            plt.ioff()  # this sometimes lets us use the navigation buttons in the toolbar  # todo is this line needed?

            # Do things depending on which plot is going to be the popout plot
            if plot == self.ax_title:
                # PDK plot
                fig, axis = plt.subplots()
                self.PDK_plot(axis)

            elif plot == self.ax1_title:
                # Pressure vs Time
                fig, ax1 = plt.subplots()
                ax12 = ax1.twinx()
                self.pressure_time_plot(data, filename, fig, ax1, ax12)

            elif plot == self.ax2_title:
                # Pressure vs Composition
                fig, ax2 = plt.subplots()
                self.pres_comp_temp_plot(fig, ax2)

            elif plot == self.ax3_title:
                # diffusivity comparison plots for optimization
                fig, ax3 = plt.subplots()
                self.comparison_plot(filename, ax3)

        else:
            showerror(message="Please select a data folder first.")
        plt.show()

    def pressure_time_plot(self, data, filename, fig1, ax1, ax12):
        """ Creates the pressure vs time (top right) plot """
        # Plot the point where the isolation valve opens and where equilibrium is determined
        ax1.plot(data.loc[self.t0[filename], 't'], data.loc[self.t0[filename], 'Pres'], 'yo', label="t$_0$")
        ax1.plot(data.loc[self.t_e[filename], 't'], data.loc[self.t_e[filename], 'Pres'], 'mo', label="t$_e$")

        # Plot pressure and set up axes for it
        ax1.plot(data['t'], data['Pres'], '.', label='Pressure')
        ax1.set_title(self.ax1_title)
        ax1.set_xlabel(self.ax1_xlabel)
        ax1.set_ylabel(self.ax1_ylabel)

        # Plot temperatures and set up axes for them
        base_colors = ['#1f77b4', '#ff7f0e', '#d62728', '#2ca02c', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']  # Red moved from fourth to third
        color_index = 2  # skip blue, orange, and red
        for i in range(1, max(self.num_GasT0, self.num_GasT) + 1):
            color_index += 1  # Move to next color in the color list
            color_index %= 10  # Cycle back to beginning if reached the end of the color list
            ax12.plot(data['t'], data['GasT' + str(i)], label='GasT' + str(i), color=base_colors[color_index])
        ax12.set_ylabel(self.ax12_ylabel)

        # Make sure ticks are in the right spot
        ax1.yaxis.set_ticks_position('left')
        ax12.yaxis.set_ticks_position('right')
        ax1.xaxis.set_ticks_position('bottom')

        # Generate and plot the red initial line (showing which values are used in initial calculations)
        # and green equilibrium line (showing which values are used in final calculations)
        xi = data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 't']
        yi = data.loc[self.t0[filename] - self.i_range[filename]:self.t0[filename] - 1, 'Pres']
        x_e = data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 't']
        y_e = data.loc[self.t_e[filename] + 1:self.t_e[filename] + self.e_range[filename], 'Pres']
        ax1.plot(xi, yi, color='red', label='Initial Points')
        ax1.plot(x_e, y_e, color='lime', label='Equilibrium Points')

        # Ensure all plots are on the legend. Simply using fig1.legend(...) doesn't get rid of old legends,
        # causing them to stack up in a way I can't figure out how to remove.
        # https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend/47370214#47370214
        pt_lines, pt_labels = ax1.get_legend_handles_labels()
        pt_lines2, pt_labels2 = ax12.get_legend_handles_labels()
        ax12.legend(pt_lines + pt_lines2, pt_labels + pt_labels2, loc="upper left", bbox_to_anchor=(0, 1),
                    bbox_transform=ax1.transAxes, framealpha=0.5)

    def pres_comp_temp_plot(self, fig2, ax2):
        """ Creates the pressure vs composition (bottom middle) plot and color codes using temperature """
        # Group data files according to their temperature to the tens place
        approx_temps = {}
        for filename in self.datafiles.keys():
            approx_temps[filename] = round(self.Ts_e[filename] - self.storage.standard_temp, -1)
        # Sort the grouped temperatures
        # https://www.geeksforgeeks.org/python-sort-python-dictionaries-by-key-or-value/
        sorted_temps = sorted(approx_temps.items(), key=lambda kv: (kv[1], kv[0]))
        sorted_names = []  # Store sorted filenames
        # Assign colors to data files according to their group
        old_temp = 0
        base_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        color_index = -1
        tempe_color = {}  # Temperature color
        legend_plots = []
        for i, tup in enumerate(sorted_temps):
            new_temp = tup[1]
            if new_temp != old_temp:  # if the temperature of the latest datafile is different from the last datafile
                color_index += 1  # Move to next color in the color list
                color_index %= 10  # Cycle back to beginning if reached the end of the color list
                legend_plots.append([i, str(new_temp) + " \u00B0C"])  # Keep track of which artists will be in the legend
            tempe_color[tup[0]] = base_colors[color_index]  # Set what color the datapoint will be
            old_temp = new_temp
            sorted_names.append(tup[0])

        artists = []
        for filename in sorted_names:
            y = self.pr_e_avg[filename]  # Pa
            x = self.HM[filename]  # Hydrogen to metal atoms ratio
            T = self.Ts_e[filename] - self.storage.standard_temp  # degrees C
            label = f'{filename}\nComposition={x:.3e} H/M' + \
                    f'\nTemperature={T:.3e} \u00B0C\nPressure={y:.3e} Pa'

            a = ax2.errorbar(x, y, yerr=self.pr_e_err[filename],
                             xerr=self.HM_err[filename], marker='o', label=label, color=tempe_color[filename])
            artists.append(a)

            # set labels on the points
        mplcursors.cursor(artists, hover=mplcursors.HoverMode.Transient).connect(
            "add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))

        # Set up axes
        ax2.set_title(self.ax2_title)
        ax2.set_xlabel(self.ax2_xlabel)
        ax2.set_ylabel(self.ax2_ylabel)
        ax2.semilogy()
        # Make sure ticks are in the right spot
        ax2.yaxis.set_ticks_position('left')
        ax2.xaxis.set_ticks_position('bottom')

        # [[x[0] for x in legend_plots]] gets the list of integers since np.array(legend_plots)[:, 0] produces an
        #      array of strings, and strings can't be used as indexes.
        # np.array(ax2.lines) turns the list of artists into an array so that a list can be used as indexes to pick out
        #      the first artist for each color/temperature and only use those in the legend.
        # np.array(legend_plots)[:, 1] produces the labels for the legend
        ax2.legend(np.array(ax2.lines)[[x[0] for x in legend_plots]], np.array(legend_plots)[:, 1],
                   loc="lower right", bbox_to_anchor=(1, 0), bbox_transform=ax2.transAxes, framealpha=0.5)

    def comparison_plot(self, filename, ax3):
        """ Creates the diffusivity optimization comparison (bottom right) plot """
        # Diffusivity comparison plots for optimization
        ax3.plot(self.D_time[filename][1:], self.lhs[filename][1:], '.', label='Experimental Data', )
        ax3.plot(self.D_time[filename][1:], self.rhs[filename][1:], label='D: Sorption')
        ax3.plot(self.D_time[filename][1:], self.rhs_cf[filename][1:], '--', label='Optimized Fit')
        ax3.set_title(self.ax3_title)
        # Set up axes
        ax3.set_xlabel(self.ax3_xlabel)
        ax3.set_ylabel(self.ax3_ylabel)
        ax3.yaxis.set_ticks_position('left')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.legend(framealpha=0.5)

    def PDK_plot(self, ax):
        """ Creates the PDK (bottom left) plot for calculated properties of different files """
        plot = self.current_variable.get()
        artists = []
        self.ax_title = plot + " vs Temperature"

        # set y labels to get units right
        if plot == "Permeability":
            self.ax_ylabel = plot + " (mol m$^{-1}$ s$^{-1}$ Pa$^{-0.5}$)"
        elif plot == "Diffusivity":
            self.ax_ylabel = plot + " (m$^2$ s$^{-1}$)"
        elif plot == "Solubility":
            self.ax_ylabel = plot + " (mol m$^{-3}$ Pa$^{-0.5}$)"

        for filename in self.datafiles.keys():
            x = self.Ts0[filename] - self.storage.standard_temp
            p = self.pr_e_avg[filename]  # pressure
            # y values change based on plot type
            if plot == "Permeability":
                y = self.Phi[filename]
                yerr = self.Phi_err[filename]
                label = f'{filename}\n{plot}={y:.3e}' + r' mol m$^{-1}\,$s$^{-1}\,$Pa$^{-0.5}$' + \
                        f'\nTemperature={x:.3e} \u00B0C\nPressure={p:.3e} Pa'
            elif plot == "Diffusivity":
                y = self.D[filename]
                yerr = self.D_err[filename]
                label = f'{filename}\n{plot}={y:.3e}' + r' m$^{2}\,$s$^{-1}$' + \
                        f'\nTemperature={x:.3e} \u00B0C\nPressure={p:.3e} Pa'
            elif plot == "Solubility":
                y = self.Ks[filename]
                yerr = self.Ks_err[filename]
                label = f'{filename}\n{plot}={y:.3e}' + r' mol m$^{-3}\,$Pa$^{-0.5}$' + \
                        f'\nTemperature={x:.3e} \u00B0C\nPressure={p:.3e} Pa'
            else:
                # This shouldn't be possible
                showerror("Selection Error", "Error, you should have selected Permeability, Diffusivity, Solubility.")
                y = self.Phi[filename]
                yerr = self.Phi_err[filename]
                label = f'{filename}\n{plot}={y:.3e}' + r' mol m$^{-1}\,$s$^{-1}\,$Pa$^{-0.5}$' + \
                        f'\nTemperature={x:.3e} K\nPressure={p:.3e} Pa'

            a = ax.errorbar(x, y, yerr=yerr, marker='o', label=label)
            artists.append(a)

        # set labels on the points
        mplcursors.cursor(artists, hover=mplcursors.HoverMode.Transient).connect(
            "add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))

        ax.semilogy()
        # ax.minorticks_off()

        ax.set_title(self.ax_title)
        ax.set_ylabel(self.ax_ylabel)
        ax.set_xlabel(self.ax_xlabel)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

    def save_figures(self):
        """ Saves figures into a folder """
        if self.loading_data:
            proceed = tk.messagebox.askyesno(title="Loading Data",
                                             message="Warning: Data is currently loading. Proceeding may cause " +
                                                     "unexpected or erroneous results. Do you wish to proceed?")
            if not proceed:
                return

        filename = asksaveasfilename(initialdir=os.path.dirname(__file__), initialfile="figures")
        if filename != '':
            # check for extension and remove it
            directory = os.path.splitext(filename)[0]  # if no extension the second element will be blank ('dir', '')
            # create the directory. Need to check if it already exists
            try:
                os.mkdir(directory)
                # save the figures
                self.fig.savefig(os.path.join(directory, self.ax_title))
                self.fig1.savefig(os.path.join(directory, self.ax1_title))
                self.fig2.savefig(os.path.join(directory, self.ax2_title))
                self.fig3.savefig(os.path.join(directory, self.ax3_title))
            except FileExistsError:
                ans = tk.messagebox.askyesno(
                    message="The directory already exists. Do you wish add the figures to this directory?")
                if ans:
                    # save the figures
                    self.fig.savefig(os.path.join(directory, self.ax_title))
                    self.fig1.savefig(os.path.join(directory, self.ax1_title))
                    self.fig2.savefig(os.path.join(directory, self.ax2_title))
                    self.fig3.savefig(os.path.join(directory, self.ax3_title))


class AdjustPVars(tk.Toplevel):
    """ popup window used to change many of the variables used in loading absorption data. """

    def __init__(self, storage, pos, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # import custom widgets
        widgets = Widgets()
        self.add_entry = widgets.add_entry

        self.title("Adjust Persistent Variables")
        self.resizable(width=False, height=False)
        self.minsize(400, 200)

        # gui_x/y values determined by running self.updateidletasks() at the end of self.__init__ and then printing size
        gui_x = 1429
        gui_y = 742
        if platform.system() == 'Darwin':
            width = 900  # width gets scaled again later if number of TCs measuring GasT has changed
            height = 424
            extra_width = 205  # width to add for each new TC beyond three
        else:
            width = 680
            height = 330
            extra_width = 165
        pos_right = int(pos[0] + (gui_x - width) / 2)
        pos_down = int(pos[1] + (gui_y - height) / 3)
        self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        # Treat closing the window via the "x" button the same as if "cancel" was clicked
        self.protocol("WM_DELETE_WINDOW", self.close_pv_win)

        self.storage = storage
        self.changed = False  # Boolean of whether edits have been made that require plot refreshing

        # Get the path for the persistent variable file
        self.pv_filename = os.path.join('data_files', 'persistent_absorption_input_variables.xlsx')

        # Read the Excel sheets into dataframes for easier access
        self.numbers_info = pd.read_excel(self.pv_filename, sheet_name="Numbers", header=0)
        self.misc_info = pd.read_excel(self.pv_filename, sheet_name="MiscInfo", header=0)
        self.GasT_info = pd.read_excel(self.pv_filename, sheet_name="GasT", header=0)
        self.SampT_info = pd.read_excel(self.pv_filename, sheet_name="SampT", header=0)
        self.Pres_info = pd.read_excel(self.pv_filename, sheet_name="Pres", header=0)

        # Pull out variables for display and potential editing

        # Number of TCs measuring gas temperature
        self.num_GasT0 = tk.DoubleVar(value=self.numbers_info['GasT'][0])
        self.num_GasT = tk.DoubleVar(value=self.numbers_info['GasT'][1])

        # Method of scaling the window width according to number of TCs measuring GasT. At num_GasT0 = 10,
        # things are still fine. Potentially add a scrollbar instead.
        if max(self.num_GasT0.get(), self.num_GasT.get()) > 2:
            width += 90 + int(extra_width * (max(self.num_GasT0.get(), self.num_GasT.get()) - 3))
            self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        # time
        self.col_t = tk.StringVar(value=self.misc_info['t'][0])  # column name
        self.m_t = tk.DoubleVar(value=self.misc_info['t'][1])  # "m" part of unit converter, as in mx + b
        self.b_t = tk.DoubleVar(value=self.misc_info['t'][2])  # "b" part of unit converter
        self.cerr_t = tk.DoubleVar(value=self.misc_info['t'][3])  # constant uncertainty
        self.perr_t = tk.DoubleVar(value=self.misc_info['t'][4])  # proportional uncertainty

        # TCs measuring gas temperature
        self.col_GasT = {}
        self.m_GasT = {}
        self.b_GasT = {}
        self.cerr_GasT = {}  # constant uncertainty
        self.perr_GasT = {}  # proportional uncertainty

        for GasT_inst in self.GasT_info.keys():
            self.col_GasT[GasT_inst] = tk.StringVar(value=self.GasT_info[GasT_inst][0])
            self.m_GasT[GasT_inst] = tk.DoubleVar(value=self.GasT_info[GasT_inst][1])
            self.b_GasT[GasT_inst] = tk.DoubleVar(value=self.GasT_info[GasT_inst][2])
            self.cerr_GasT[GasT_inst] = tk.DoubleVar(value=self.GasT_info[GasT_inst][3])
            self.perr_GasT[GasT_inst] = tk.DoubleVar(value=self.GasT_info[GasT_inst][4])

        # Sample temperature
        self.col_SampT = tk.StringVar(value=self.SampT_info['SampT'][0])
        self.m_SampT = tk.DoubleVar(value=self.SampT_info['SampT'][1])
        self.b_SampT = tk.DoubleVar(value=self.SampT_info['SampT'][2])
        self.cerr_SampT = tk.DoubleVar(value=self.SampT_info['SampT'][3])
        self.perr_SampT = tk.DoubleVar(value=self.SampT_info['SampT'][4])

        # Pressure
        self.col_Pres = tk.StringVar(value=self.Pres_info['Pres'][0])
        self.m_Pres = tk.DoubleVar(value=self.Pres_info['Pres'][1])
        self.b_Pres = tk.DoubleVar(value=self.Pres_info['Pres'][2])
        self.cerr_Pres = tk.DoubleVar(value=self.Pres_info['Pres'][3])
        self.perr_Pres = tk.DoubleVar(value=self.Pres_info['Pres'][4])

        # Isolation valve
        self.col_IV = tk.StringVar(value=self.misc_info['Isolation Valve'][0])
        self.m_IV = tk.DoubleVar(value=self.misc_info['Isolation Valve'][1])
        self.b_IV = tk.DoubleVar(value=self.misc_info['Isolation Valve'][2])

        # Row at which the program starts reading data
        self.starting_row = tk.IntVar(value=self.misc_info['Starting Row'][0])
        # Rows at the end of the file which the program doesn't read
        self.footer_rows = tk.IntVar(value=self.misc_info['Rows in Footer'][0])

        # store entries
        self.inputs = {}

        # create frames
        self.create_labels()
        self.create_inputs()

        submit_button = ttk.Button(self, text="Submit", command=self.submit)
        submit_button.grid(row=14, column=7, columnspan=2, sticky='e')

        quit_button = ttk.Button(self, text="Cancel", command=self.close_pv_win)
        quit_button.grid(row=14, column=9, columnspan=2, sticky='w')

        parent = tk.Frame(self)
        parent.grid(row=14, column=0, columnspan=7, sticky="nsew", pady=5)
        self.add_entry(self, parent, variable=self.inputs, key='num_GasT0', text="Number of TCs measuring GasT0:",
                       subscript='', ent_w=3, tvar=self.num_GasT0, units='', row=0, column=0, in_window=True,
                       command=lambda tvar, variable, key, pf: self.one_to_max(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key='num_GasT', text="Number of TCs measuring GasT:",
                       subscript='', ent_w=3, tvar=self.num_GasT, units='', row=0, column=4, in_window=True,
                       command=lambda tvar, variable, key, pf: self.one_to_max(tvar, variable, key, pf))

    def one_to_max(self, tvar, variable, key, pf=None):
        """ Ensure the user entry is greater than 1. If not, set it equal to 1. """
        self.storage.check_for_number(tvar, variable, key, parent_frame=pf)
        maxTCnum = 8  # Max number of TCs to allow
        # Either set the entry to 1 if less than 1, to int(maxTCnum) if more than int(maxTCnum), or to the user's entry
        if key in variable.keys() and tvar.get() < 1:
            tk.messagebox.showwarning("Invalid Entry", "Please enter a number greater than or equal to 1.")
            self.deiconify()
            variable[key].delete(0, "end")
            variable[key].insert(0, "1")
            tvar.set(1)
        elif key in variable.keys() and tvar.get() > int(maxTCnum):
            # Keeps the user from having to type out too many unique column names if something like "100" gets submitted
            tk.messagebox.showwarning("Invalid Entry", 'Please enter a number less than or equal to ' +
                                                       '{}.\n\n'.format(int(maxTCnum)) + 'If more than ' +
                                                       '{} TCs are needed,'.format(int(maxTCnum)) + ' simply search ' +
                                                       'for "maxTCnum" in the source code and change the ' +
                                                       '{} in its definition '.format(int(maxTCnum)) +
                                                       'to your desired number.')
            self.deiconify()
            variable[key].delete(0, "end")
            variable[key].insert(0, "{}".format(int(maxTCnum)))
            tvar.set(int(maxTCnum))
        else:
            variable[key].delete(0, "end")
            variable[key].insert(0, "{}".format(int(tvar.get())))
            tvar.set(int(tvar.get()))
        return True

    def create_labels(self):
        """ Creates intro text to Absorption Plots tab's Adjust Persistent Variables window """
        parent = tk.Frame(self)
        parent.grid(row=0, column=0, columnspan=10, sticky="nsew")
        tk.Label(parent, text='Choose how the program will handle your uploaded data. ' +
                              'Click the "?" for symbol information:').grid(sticky="nsew", columnspan=2)

        # Create a nicely formatted help button
        s = ttk.Style()
        s.configure('Bold.TButton', font=('Helvetica', 10, 'bold'))
        help_button = ttk.Button(parent, text=" ? ", width=3, style='Bold.TButton', command=self.show_help)
        help_button.grid(row=0, column=10)

    def show_help(self):
        """ Calls the class for creating the help window """
        x = self.winfo_rootx()
        y = self.winfo_rooty()
        apv_width = self.winfo_width()
        apv_height = self.winfo_height()
        APSettingsHelp(pos=(x, y), size=(apv_width, apv_height)).show_help_window()

    def check_col_names(self, tvar, variable, key):
        """ Checks to make sure the entered string is made up of only letters"""
        if variable[key].get().isalpha():
            svar = variable[key].get().upper()  # set the string to uppercase
        else:
            tk.messagebox.showwarning("Invalid Entry", "Please enter a string of letters.")
            self.deiconify()
            svar = tvar.get()  # reset to what it was before the invalid entry was typed

        # Either set the entry back to what it was (if not a string of letters) or to the string of (uppercase) letters
        if key in variable.keys():
            variable[key].delete(0, "end")
            variable[key].insert(0, "{}".format(svar))
            tvar.set(svar)
        return True

    def create_inputs(self):
        """Create the many entry boxes this window requires"""
        parent = self
        entry_w = 8  # Width of each entry box in this window
        entry_row = 1  # First row of the entries

        entry_col = 0
        tk.Label(parent, text="t [s]").grid(row=entry_row + 0, column=1)
        self.add_entry(self, parent, variable=self.inputs, key="col_t", text="col", innersubscript="t", innertext=":",
                       subscript="", tvar=self.col_t, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_t", text="m", innersubscript="t", innertext=":",
                       subscript="", tvar=self.m_t, units="[s/?]", ent_w=entry_w,
                       row=entry_row + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_t", text="b", innersubscript="t", innertext=":",
                       subscript="", tvar=self.b_t, units="[s]", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="cerr_t", text="cerr", innersubscript="t", innertext=":",
                       subscript="", tvar=self.cerr_t, units="[s]", ent_w=entry_w,
                       row=entry_row + 4, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="perr_t", text="perr", innersubscript="t", innertext=":",
                       subscript="", tvar=self.perr_t, units="[%]", ent_w=entry_w,
                       row=entry_row + 5, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        entry_col += 3
        tk.Label(parent, text="Pres [Pa]").grid(row=entry_row + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="col_Pres", text="col", innersubscript="Pres",
                       innertext=":", subscript="", tvar=self.col_Pres, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_Pres", text="m", innersubscript="Pres",
                       innertext=":", subscript="", tvar=self.m_Pres, units="[Pa/?]", ent_w=entry_w,
                       row=entry_row + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_Pres", text="b", innersubscript="Pres",
                       innertext=":", subscript="", tvar=self.b_Pres, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="cerr_Pres", text="cerr", innersubscript="Pres",
                       innertext=":", subscript="", tvar=self.cerr_Pres, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 4, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="perr_Pres", text="perr", innersubscript="Pres",
                       innertext=":", subscript="", tvar=self.perr_Pres, units="[%]", ent_w=entry_w,
                       row=entry_row + 5, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        entry_col += 3
        tk.Label(parent, text="Isolation Valve [1/0]").grid(row=entry_row + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="col_IV", text="col", innersubscript="IV",
                       innertext=":", subscript="", tvar=self.col_IV, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_IV", text="m", innersubscript="IV",
                       innertext=":", subscript="", tvar=self.m_IV, units="", ent_w=entry_w,
                       row=entry_row + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_IV", text="b", innersubscript="IV",
                       innertext=":", subscript="", tvar=self.b_IV, units="", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        entry_col += 3
        tk.Label(parent, text="Starting Row").grid(row=entry_row + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="Starting Row", text="", innersubscript="",
                       innertext="", subscript="", tvar=self.starting_row, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        tk.Label(parent, text="Rows in Footer").grid(row=entry_row + 2, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="Footer Rows", text="", innersubscript="",
                       innertext="", subscript="", tvar=self.footer_rows, units="", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        # Second row of instruments

        row2entries = entry_row + 7
        entry_col = -3  # Reset to far left (noting that it will start with a += 3)
        tk.Label(parent, text="").grid(row=row2entries-1, column=0)  # Add some space between the rows

        # Add entries for arbitrarily many instruments measuring GasT
        for GasT_inst in self.GasT_info.keys():
            entry_col += 3
            tk.Label(parent, text="{} [\u00B0C]".format(GasT_inst)).grid(row=row2entries, column=entry_col + 1)
            self.add_entry(self, parent, variable=self.inputs, key="col_{}".format(GasT_inst), text="col",
                           innersubscript="{}".format(GasT_inst), innertext=":",
                           subscript="", tvar=self.col_GasT[GasT_inst], units="", ent_w=entry_w,
                           row=row2entries + 1, column=entry_col,
                           command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
            self.add_entry(self, parent, variable=self.inputs, key="m_{}".format(GasT_inst), text="m",
                           innersubscript="{}".format(GasT_inst), innertext=":",
                           subscript="", tvar=self.m_GasT[GasT_inst], units="[\u00B0C/?]", ent_w=entry_w,
                           row=row2entries + 2, column=entry_col, in_window=True,
                           command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                                 parent_frame=pf))
            self.add_entry(self, parent, variable=self.inputs, key="b_{}".format(GasT_inst), text="b",
                           innersubscript="{}".format(GasT_inst), innertext=":",
                           subscript="", tvar=self.b_GasT[GasT_inst], units="[\u00B0C]", ent_w=entry_w,
                           row=row2entries + 3, column=entry_col, in_window=True,
                           command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                                 parent_frame=pf))
            self.add_entry(self, parent, variable=self.inputs, key="cerr_{}".format(GasT_inst), text="cerr",
                           innersubscript="{}".format(GasT_inst), innertext=":",
                           subscript="", tvar=self.cerr_GasT[GasT_inst], units="[\u00B0C]", ent_w=entry_w,
                           row=row2entries + 4, column=entry_col, in_window=True,
                           command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                                 parent_frame=pf))
            self.add_entry(self, parent, variable=self.inputs, key="perr_{}".format(GasT_inst), text="perr",
                           innersubscript="{}".format(GasT_inst), innertext=":",
                           subscript="", tvar=self.perr_GasT[GasT_inst], units="[%]", ent_w=entry_w,
                           row=row2entries + 5, column=entry_col, in_window=True,
                           command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                                 parent_frame=pf))

        # Sample temperature
        entry_col += 3
        tk.Label(parent, text="SampT [\u00B0C]").grid(row=row2entries, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="col_SampT", text="col", innersubscript="SampT",
                       innertext=":", subscript="", tvar=self.col_SampT, units="", ent_w=entry_w,
                       row=row2entries + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_SampT", text="m", innersubscript="SampT",
                       innertext=":", subscript="", tvar=self.m_SampT, units="[\u00B0C/?]", ent_w=entry_w,
                       row=row2entries + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_SampT", text="b", innersubscript="SampT",
                       innertext=":", subscript="", tvar=self.b_SampT, units="[\u00B0C]", ent_w=entry_w,
                       row=row2entries + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="cerr_SampT", text="cerr", innersubscript="SampT",
                       innertext=":", subscript="", tvar=self.cerr_SampT, units="[\u00B0C]", ent_w=entry_w,
                       row=row2entries + 4, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="perr_SampT", text="perr", innersubscript="SampT",
                       innertext=":", subscript="", tvar=self.perr_SampT, units="[%]", ent_w=entry_w,
                       row=row2entries + 5, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

    def submit(self):
        """Updates the persistent variables with user input"""

        # Create a list that contains all non-column name entries in GasT for checking for NaNs
        GasT_vals = list(self.m_GasT.values()) + list(self.b_GasT.values()) + \
            list(self.cerr_GasT.values()) + list(self.perr_GasT.values())

        # List of all column names for checking for uniqueness
        col_list = [self.col_t.get(), self.col_SampT.get(), self.col_Pres.get(), self.col_IV.get()]
        for val in self.col_GasT.values():
            col_list.append(val.get())
        # check to make sure no numerical variable is assigned "nan"
        if sum(map(np.isnan, [self.m_t.get(), self.b_t.get(), self.cerr_t.get(), self.perr_t.get(),
                              self.m_SampT.get(), self.b_SampT.get(), self.cerr_SampT.get(), self.perr_SampT.get(),
                              self.m_Pres.get(), self.b_Pres.get(), self.cerr_Pres.get(), self.perr_Pres.get(),
                              self.m_IV.get(), self.b_IV.get(), self.starting_row.get(), self.footer_rows.get(),
                              self.num_GasT0.get(), self.num_GasT.get()])) + \
                              any(np.isnan(val.get()) for val in GasT_vals) > 0:
            tk.messagebox.showwarning(title="Missing Entry",
                                      message="Please fill in remaining boxes before continuing.")
            self.deiconify()

        elif len(set(col_list)) < len(col_list):  # Check to make sure all column names are unique
            tk.messagebox.showwarning(title="Duplicate Column",
                                      message="Please ensure all column names are unique.")
            self.deiconify()

        else:
            ans = tk.messagebox.askokcancel(
                title="Confirmation",
                message="Do you wish to save the variables as they are to the source file?")
            if ans:
                # Redefine dataframes to account for changes
                self.misc_info['t'] = [self.col_t.get(), self.m_t.get(), self.b_t.get(),
                                       self.cerr_t.get(), self.perr_t.get()]
                for GasT_inst in self.GasT_info.keys():
                    self.GasT_info[GasT_inst] = [self.col_GasT[GasT_inst].get(), self.m_GasT[GasT_inst].get(),
                                                 self.b_GasT[GasT_inst].get(),
                                                 self.cerr_GasT[GasT_inst].get(), self.perr_GasT[GasT_inst].get()]
                self.SampT_info['SampT'] = [self.col_SampT.get(), self.m_SampT.get(), self.b_SampT.get(),
                                            self.cerr_SampT.get(), self.perr_SampT.get()]
                self.Pres_info['Pres'] = [self.col_Pres.get(), self.m_Pres.get(), self.b_Pres.get(),
                                          self.cerr_Pres.get(), self.perr_Pres.get()]
                self.misc_info['Isolation Valve'] = [self.col_IV.get(), self.m_IV.get(),
                                                     self.b_IV.get(), np.nan, np.nan]
                self.misc_info['Starting Row'] = [self.starting_row.get(), np.nan, np.nan, np.nan, np.nan]
                self.misc_info['Rows in Footer'] = [self.footer_rows.get(), np.nan, np.nan, np.nan, np.nan]
                # Create or remove data from dataframe so that it has the max number of instruments measuring GasT
                # out of self.num_GasT0 and self.num_GasT. Note that this can handle non-integer values
                # of self.num_GasT0 and self.num_GasT
                max_num = max(self.numbers_info['GasT'][0], self.numbers_info['GasT'][1])
                while max_num > max(self.num_GasT0.get(), self.num_GasT.get()):
                    del self.GasT_info['GasT' + str(max_num)]
                    max_num -= 1
                while max_num < max(self.num_GasT0.get(), self.num_GasT.get()):
                    self.GasT_info['GasT' + str(max_num + 1)] = [self.col_GasT['GasT1'].get(), 1.0, 0.0, 0.0, 0.0]
                    max_num += 1
                self.numbers_info.loc[0, 'GasT'] = self.num_GasT0.get()
                self.numbers_info.loc[1, 'GasT'] = self.num_GasT.get()

                # overwrite the old source file with updated dataframes
                with pd.ExcelWriter(self.pv_filename, engine="openpyxl", if_sheet_exists='new', mode='a') as writer:
                    workBook = writer.book

                    # remove's shouldn't be required, but if_sheet_exists='replace' is malfunctioning as of 11/18/21
                    workBook.remove(workBook['Numbers'])
                    workBook.remove(workBook['MiscInfo'])
                    workBook.remove(workBook['GasT'])
                    workBook.remove(workBook['SampT'])
                    workBook.remove(workBook['Pres'])

                    # Write updated dataframes to the Excel sheets
                    self.numbers_info.to_excel(writer, sheet_name='Numbers', index=False)
                    self.misc_info.to_excel(writer, sheet_name='MiscInfo', index=False)
                    self.GasT_info.to_excel(writer, sheet_name='GasT', index=False)
                    self.SampT_info.to_excel(writer, sheet_name='SampT', index=False)
                    self.Pres_info.to_excel(writer, sheet_name='Pres', index=False)
                    writer.save()

                self.changed = True
                self.destroy()  # close the pop-up window
            else:
                self.deiconify()  # bring back the pop-up window

    def close_pv_win(self):
        """ Checks for unsaved updates to persistent variables """

        # Create new dataframes to compare against
        temp_misc_info = self.misc_info.copy()
        temp_GasT_info = self.GasT_info.copy()
        temp_SampT_info = self.SampT_info.copy()
        temp_Pres_info = self.Pres_info.copy()
        temp_numbers_info = self.numbers_info.copy()

        # update new dataframes with current entries
        temp_misc_info['t'] = [self.col_t.get(), self.m_t.get(), self.b_t.get(),
                               self.cerr_t.get(), self.perr_t.get()]
        for GasT_inst in self.GasT_info.keys():
            temp_GasT_info[GasT_inst] = [self.col_GasT[GasT_inst].get(), self.m_GasT[GasT_inst].get(),
                                         self.b_GasT[GasT_inst].get(),
                                         self.cerr_GasT[GasT_inst].get(), self.perr_GasT[GasT_inst].get()]
        temp_SampT_info['SampT'] = [self.col_SampT.get(), self.m_SampT.get(), self.b_SampT.get(),
                                    self.cerr_SampT.get(), self.perr_SampT.get()]
        temp_Pres_info['Pres'] = [self.col_Pres.get(), self.m_Pres.get(), self.b_Pres.get(),
                                  self.cerr_Pres.get(), self.perr_Pres.get()]
        temp_misc_info['Isolation Valve'] = [self.col_IV.get(), self.m_IV.get(), self.b_IV.get(), np.nan, np.nan]
        temp_misc_info['Starting Row'] = [self.starting_row.get(), np.nan, np.nan, np.nan, np.nan]
        temp_misc_info['Rows in Footer'] = [self.footer_rows.get(), np.nan, np.nan, np.nan, np.nan]
        temp_numbers_info['GasT'] = [int(self.num_GasT0.get()), int(self.num_GasT.get())]  # "int" is so HyPAT recognizes when there was no change

        # List of all column names for checking for uniqueness
        col_list = [self.col_t.get(), self.col_SampT.get(), self.col_Pres.get(), self.col_IV.get()]
        for val in self.col_GasT.values():
            col_list.append(val.get())

        # Check to see if any changes have been made and if there are any duplicate column names
        if temp_misc_info.equals(self.misc_info) and temp_GasT_info.equals(self.GasT_info) and \
                temp_SampT_info.equals(self.SampT_info) and \
                temp_Pres_info.equals(self.Pres_info) and temp_numbers_info.equals(self.numbers_info)\
                and len(set(col_list)) == len(col_list):
            self.destroy()
        elif len(set(col_list)) < len(col_list):  # If there are duplicate column names
            tk.messagebox.showwarning(title="Duplicate Column",
                                      message="Please ensure all column names are unique.")
            self.deiconify()
        else:  # If changes have been made, check with user before closing
            ans = tk.messagebox.askokcancel(
                title="Confirmation",
                message="Do you wish close without saving changes?")
            if ans:
                self.destroy()
            else:
                self.deiconify()

    def show_pv_window(self):
        """ This allows the window to open and then, when closed, trigger the main gui to update the plots if
            changes were made """
        self.deiconify()
        self.wait_window()

        return self.changed


class APSettingsHelp(tk.Toplevel):
    """ popup box to display information about the variables that can be edited. """

    def __init__(self, pos, size, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title("Help for Absorption's Persistent Variables Settings")
        self.resizable(width=False, height=False)
        self.minsize(400, 200)

        if platform.system() == 'Darwin':
            width = 900
            height = 500
        else:
            width = 650
            height = 430
        pos_right = int(pos[0] + (size[0] - width) / 2)
        pos_down = int(pos[1] + (size[1] - height) / 2)
        self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        self.create_text_boxes()

    def create_text_boxes(self):
        """ Create the instructions that appear in this help window """
        # Get fonts ready for use
        from tkinter import font
        tfont = font.nametofont("TkDefaultFont").copy()
        btfont, bbtfont, bb2tfont = tfont.copy(), tfont.copy(), tfont.copy()
        btfont.config(weight="bold")
        bbtfont.config(weight="bold", underline=True)

        # Create a text widget for each symbol, along with one (" ") for the header
        td = {}
        tdt = {}
        for item in (" ", "t [s]: ", "Pres [Pa]: ", "Isolation Valve [1/0]: ", "Starting Row: ", "Rows in Footer: ",
                     "GasT [\u00B0C]: ", "SampT [\u00B0C]: ", "col: ", "m: ", "b: ", "cerr: ", "perr: "):
            td[item] = tk.Text(self, borderwidth=0, background=self.cget("background"), spacing3=1)
            tdt[item + "text"] = tk.Text(self, borderwidth=0, background=self.cget("background"), spacing3=1)

        # Add in formatting options
        for symbol in td:
            td[symbol].tag_configure("bold", justify="right", font=btfont)
            tdt[symbol + "text"].tag_configure("bold_underline", font=bbtfont)

        # Insert text into right-side widgets
        tdt[" text"].insert("insert", "Symbols", "bold_underline")  # Header
        tdt["t [s]: text"].insert("insert", "Time. HyPAT will treat this column either as a column of datetimes or " +
                                  "as cumulative time passed since start of data recording, " +
                                  "depending on whether the first entry in the column is a datetime or a float/int. " +
                                  "The variables 'm' and 'b,' described further down, " +
                                  "are only applied to the time data if the data is composed of floats.")
        tdt["Pres [Pa]: text"].insert("insert", "Pressure.")
        tdt["Isolation Valve [1/0]: text"].insert("insert", "The data describing whether the isolation valve is " +
                                                  "closed (~0) or open (not ~0). Start of absorption is " +
                                                  "measured from when the valve first goes from closed to " +
                                                  "open.")
        tdt["Starting Row: text"].insert("insert", "The row at which the program starts reading data from the Excel " +
                                                   "sheet (0-indexed). HyPAT assumes the data file has no header, so " +
                                                   "use this to start analysis after the header.")
        tdt["Rows in Footer: text"].insert("insert", "How many rows at the bottom of the file to ignore.")
        tdt["GasT [\u00B0C]: text"].insert("insert", "The thermocouples (TCs) that measure the gas's temperature. " +
                                                    "The number of TCs can be set arbitrarily via " +
                                                    "'Number of TCs measuring GasT0' and 'Number of TCs measuring " +
                                                    "GasT,' each of which start counting with GasT1.")
        tdt["SampT [\u00B0C]: text"].insert("insert", "The thermocouple that measures the sample temperature.")
        tdt["col: text"].insert("insert", "Reference column for corresponding quantity in raw data file.")
        tdt["m: text"].insert("insert", "Every entry in the column of data is converted according to m*x + b, " +
                                        "where x is the entry in the column. Use this to convert the column's data " +
                                        "to SI units.")
        tdt["b: text"].insert("insert", "See entry for 'm' above.")
        tdt["cerr: text"].insert("insert", "Constant uncertainty of the instrument, e.g., +/- 2\u00B0C. " +
                                           "When calculating the uncertainty caused by each instrument, the " +
                                           "application chooses the largest out of the constant uncertainty, the " +
                                           "proportional uncertainty (see below), and the calculated statistical " +
                                           "uncertainty.")
        tdt["perr: text"].insert("insert", "Proportional uncertainty of the instrument, e.g., " +
                                           "+/- 0.75% of the measurement (in \u00B0C).")

        # Configure each text widget
        for i, symbol in enumerate(td):
            # Set width (and other things) of left side
            if platform.system() == 'Darwin':
                divisor = 9
                addend = 2
                divisor2 = 9
                addend2 = 0
            else:
                divisor = 7
                addend = 6
                divisor2 = 6
                addend2 = 0
            # btfont.measure() says how many pixels text+subscript would take to write out.
            # For Windows, //7 takes that number and turns it into how many zeros it would take to be that long because
            #   one zero takes up 7 pixels in this font. The +6 adds one more character to make up for // rounding down
            #   and 5 more characters to make room.
            # For Mac, setting divisor to 9 and addend to 2 seems to look good, but it's unclear why those numbers
            #   work well
            # The unit of width is characters, specifically 0s, not pixels, which is why "Isolation Valve [1/0]" is
            #   needed. Those words are the longest text on the left side.
            td[symbol].configure(width=btfont.measure("Isolation Valve [1/0]") // divisor + addend,
                                 wrap=tk.WORD, height=1)

            # Set height of each right-side widget according to how many lines it will take to display
            text_width = btfont.measure(tdt[symbol + "text"].get("1.0", 'end-1c')) // divisor2 + addend2
            text_height = 1
            while text_width >= tdt[symbol + "text"].cget("width"):  # Check if the text width > widget's width.
                text_height += 1
                text_width -= tdt[symbol + "text"].cget("width")
            tdt[symbol + "text"].configure(wrap=tk.WORD, height=text_height)

            td[symbol].insert("insert", symbol, "bold")  # Insert text into the left-side widgets

            td[symbol].grid(row=i, column=0, sticky="nse")
            tdt[symbol + "text"].grid(row=i, column=1, sticky="nsew")

            # Make the text un-editable (to user and to program)
            td[symbol].configure(state="disabled")
            tdt[symbol + "text"].configure(state="disabled")

    def show_help_window(self):
        """ This allows the window to open correctly """
        self.deiconify()
        self.wait_window()
        return
