""" Code for the Permeation Plots tab in HyPAT """
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
from .data_storage import Widgets, LoadingScreen
from scipy.optimize import curve_fit
import mplcursors  # adds the hover-over feature to labels. This library has good documentation
import openpyxl
import platform  # allows for Mac vs. Windows adaptions
# make certain warnings appear as errors, allowing them to be caught using a try/except clause
import warnings
from scipy.optimize import OptimizeWarning
from scipy.signal import savgol_filter
warnings.simplefilter("error", OptimizeWarning)


class PermeationPlots(tk.Frame):

    def __init__(self, parent, storage, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        # container for important variables
        self.storage = storage
        self.datafiles = {}  # Store Excel files

        self.header = []  # Container for titles of each column in your data file that will be used
        self.needed_cols = []  # Container that will hold the column number for each needed data set
        self.converters_dict = {}  # Container for functions for converting data into correct units

        self.num_GasT = 1  # Number of TCs measuring GasT

        # When determining uncertainty, the program will pick the max of the calculated statistical uncertainty,
        # the constant uncertainty, and the proportional uncertainty. These are containers for the constant and
        # proportional uncertainties of each instrument
        self.GasT_cerr = {}  # degrees C, i.e., +/- 2.2 degrees C (constant uncertainty)
        self.GasT_perr = {}  # percentage, i.e., +/- 0.75% of total measurement (in Celsius) (proportional uncertainty)
        self.SampT_cerr = {}
        self.SampT_perr = {}
        self.SecP_cerr = {}  # Pa (constant uncertainty)
        self.SecP_perr = {}  # percentage, (proportional uncertainty)
        self.PrimP_cerr = {}
        self.PrimP_perr = {}

        # store widgets
        self.widgets = Widgets()
        self.add_text0 = self.widgets.add_text0
        self.add_text = self.widgets.add_text
        self.add_text2 = self.widgets.add_text2
        self.add_entry = self.widgets.add_entry
        self.add_entry3 = self.widgets.add_entry3

        # to get 2 rows over 3, we need two more frames
        self.top_frame = tk.Frame(self, bd=10)
        self.top_frame.grid(row=0, column=0, sticky="nsew")
        self.bottom_frame = tk.Frame(self, bd=10)
        self.bottom_frame.grid(row=1, column=0, sticky="nsew")

        # user inputs
        self.frame = tk.LabelFrame(self.top_frame, bd=8)
        self.frame.grid(row=0, column=0, sticky="nsew")
        self.inputs = {}

        self.directory = ""
        self.loading_data = False  # Keep track of whether data is currently being loaded by this tab
        self.refreshing = False  # If false, ask for a directory when select_file gets called
        self.file_type = ""
        self.time_type = 0  # Used to determine which converter function to use for the time instrument
        self.days = -1  # How many days have passed during the experiment (used for datetime conversion to seconds)
        self.this_date_time = 0  # Current datetime (used for datetime conversion to seconds)
        self.init_time = 0  # Initial datetime in seconds (used for datetime conversion to seconds)
        self.extra_time = 0  # Extra time to be accounted for (used for datetime conversion to seconds)
        self.error_texts = ""  # Text variable for storing the error texts that may come up when files are loaded

        # arrange frames in tab
        self.rowconfigure((0, 1), weight=1)
        self.columnconfigure((0, 1), weight=1)
        self.top_frame.rowconfigure(0, weight=1)
        self.top_frame.columnconfigure((0, 1), weight=1)
        self.bottom_frame.rowconfigure(0, weight=1)
        self.bottom_frame.columnconfigure((0, 1, 2), weight=1)

        # variables for input
        self.sample_thickness = tk.DoubleVar(value=1)
        self.volume = tk.DoubleVar(value=(self.storage.defaults_info["Secondary Side Volume [cc]"][0]) * 10 ** -6)
        self.ss_tol = tk.DoubleVar(value=self.storage.tol.get())  # Tolerance for finding the steady state
        self.ss_t_del = tk.DoubleVar(value=self.storage.t_del.get())  # Minimum time after t0 before looking for ss
        self.leak_range = {}  # Dict for holding the range used in finding leaks and initial values for each file
        self.ss_range = {}  # Dict for holding the range used in finding ss and ss values for each file
        self.gen_leak_range = tk.IntVar(value=self.storage.gen_dp_range.get())  # Default range to determine leak over
        self.gen_ss_range = tk.IntVar(value=self.storage.gen_dp_range.get())  # Default range to determine ss over
        # default area is the inner area of the default O-ring
        self.A_perm = tk.DoubleVar()
        inner_diameter = self.storage.oring_info.loc[self.storage.oring_info.index[0],
                                                     self.storage.oring_info.columns[1]]
        A_perm_str = "{:.2e}".format(
            np.pi * (inner_diameter / 2 * 0.001) ** 2)  # Provides a way to round to sig figs
        self.A_perm.set(float(A_perm_str))
        # Variables for storing uncertainty, rounded to avoid numerical errors
        self.volume_err = tk.DoubleVar(value=round(self.volume.get()*0.05, 13))
        self.A_perm_err = tk.DoubleVar(value=round(self.A_perm.get()*0.05, 13))
        self.sample_thickness_err = tk.DoubleVar(value=round(self.sample_thickness.get()*0.05, 13))

        # variables for filtering
        self.poly_deg = tk.IntVar(value=2) #number of polynomial for savitzky-golay filter
        self.filter_win = tk.IntVar(value=20) #window size for filter
        self.filter = True #boolean that controls the Savgol filter
        self.rolling = False #boolean that controls the rolling average toggle switch

        # permeation variables
        self.t0 = {}  # time when isolation valve opens (s)
        self.tss = {}  # time when steady state pressure is achieved (s)
        self.pleak_lf = {}  # linear fit for leak rate
        self.pss_lf = {}  # linear fit for steady state
        self.Prate = {}  # steady state permeation rates (pss_lf minus pleak_lf)
        self.Prate_err = {}
        self.F = {}  # flux
        self.F_err = {}
        self.Phi = {}  # permeability
        self.Phi_err = {}
        self.Perm_Plot = {}  # real time permeability (for plotting against time) adjusted for the pressure leak
        self.Tg0 = {}  # Temperature of gas at t0 (using an average over time and instruments) (K)
        self.Tgss = {}  # Temperature of gas at tss (using an average over time and instruments) (K)
        self.Tgss_err = {}
        self.artists = []  # store points for Perm/Dif/Sol vs. Temp graph
        self.Tsamp0 = {}  # Sample temperature at t0 (using an average over time and instruments) (K)
        self.Tsampss = {}  # Sample temperature at tss (using an average over time and instruments) (K)
        self.Tsampss_err = {}
        self.SecP_average = {}  # Secondary Side Pressure (using an average over time)
        self.SecP_err = {}
        self.PrimP_average = {}  # Primary Side Pressure (using an average over time)
        self.PrimP_err = {}

        # diffusivity variables
        self.intercept = {}
        self.tlag = {}
        self.D = {}  # diffusivity
        self.D_err = {}
        self.D_tlag = {}
        self.D_tlag_err = {}
        self.A = {}  # proportionality constant
        self.dt = {}  # additive time constant
        self.D_time = {}  # time over which D is calculated
        self.lhs = {}  # for the diffusivity optimization comparison
        self.rhs = {}  # time-lag method
        self.rhs_cf = {}  # curve_fit

        # solubility variables
        self.Ks = {}
        self.Ks_err = {}

        # text variables for labels
        self.Prate_label = tk.DoubleVar(value=0)
        self.F_label = tk.DoubleVar(value=0)
        self.Phi_label = tk.DoubleVar(value=0)
        self.D_label = tk.DoubleVar(value=0)
        self.D_tlag_label = tk.DoubleVar(value=0)
        self.K_label = tk.DoubleVar(value=0)
        self.tss_label = tk.DoubleVar(value=0)
        self.Prate_err_label = tk.DoubleVar(value=0)
        self.F_err_label = tk.DoubleVar(value=0)
        self.Phi_err_label = tk.DoubleVar(value=0)
        self.D_err_label = tk.DoubleVar(value=0)
        self.D_tlag_err_label = tk.DoubleVar(value=0)
        self.K_err_label = tk.DoubleVar(value=0)

        # add frames and buttons to view important variables

        entry_row = 1
        self.add_entry3(self, self.frame, variable=self.inputs, key="ls", text="Sample Thickness: l",
                        subscript="s", tvar1=self.sample_thickness, tvar2=self.sample_thickness_err,
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        units="[mm]", row=entry_row)

        self.add_entry3(self, self.frame, variable=self.inputs, key="A_perm", text="Permeation Surface Area: A",
                        subscript="perm", tvar1=self.A_perm, tvar2=self.A_perm_err,
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        units="[m\u00b2]", row=entry_row + 1)  # "[m^2]"

        self.add_entry3(self, self.frame, variable=self.inputs, key="V", text="Secondary Side Volume: V",
                        subscript="sec", tvar1=self.volume, tvar2=self.volume_err, units="[m\u00B3]",  # "[m^3]"
                        update_function=lambda tvar, var_type, key: self.update_function(tvar, var_type, key),
                        row=entry_row + 2)

        label_row = 4
        self.add_text2(self.frame, text="Permeation Rate: dP/dt", subscript="", tvar1=self.Prate_label,
                       tvar2=self.Prate_err_label, units="[Pa s\u207b\u00b9]", row=label_row)
        self.add_text2(self.frame, text="Molar Flux: J", subscript="inf", tvar1=self.F_label, tvar2=self.F_err_label,
                       units="[mol m\u207b\u00b2 s\u207b\u00b9]", row=label_row + 1)  # "[mol/m^2 s]"
        self.add_text2(self.frame, text="Permeability: \u03A6", subscript="", tvar1=self.Phi_label, tvar2=self.Phi_err_label,
                       units="[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]",   # "[mol/msPa^0.5]"
                       row=label_row + 2)
        self.add_text2(self.frame, text="Diffusivity (Optimized): D", subscript="", tvar1=self.D_label, tvar2=self.D_err_label,
                       units="[m\u00b2 s\u207b\u00b9]", row=label_row + 3)
        self.add_text2(self.frame, text="Diffusivity (timelag): D", subscript="", tvar1=self.D_tlag_label, tvar2=self.D_tlag_err_label,
                       units="[m\u00b2 s\u207b\u00b9]", row=label_row + 4)
        self.add_text2(self.frame, text="Solubility: K", subscript="s", tvar1=self.K_label, tvar2=self.K_err_label,
                       units="[mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]", row=label_row + 5)
        self.add_text(self.frame, text="Time to Steady State: t", subscript="ss", tvar=self.tss_label,
                       units="[s]", row=label_row + 6)

        button_row = 10
        self.b0 = ttk.Button(self.frame, text='Choose folder', command=self.select_file)
        self.b0.grid(row=button_row + 2, column=4, sticky="ew")
        b1 = ttk.Button(self.frame, text='Refresh', command=self.refresh_graphs)
        b1.grid(row=button_row, column=4, sticky="ew")

        settings_b = ttk.Button(self.frame, text='Settings', command=self.adjust_persistent_vars)
        settings_b.grid(row=button_row + 1, column=4, sticky="ew")

        b2 = ttk.Button(self.frame, text='Close popout plots', command=lambda: plt.close('all'))
        b2.grid(row=button_row + 1, column=1, sticky="ew")

        self.coord_b = ttk.Button(self.frame, text="Enable Coordinates", command=self.toggle_coordinates)
        self.coord_b.grid(row=button_row, column=3, sticky="ew")

        b3 = ttk.Button(self.frame, text='Save current figures', command=self.save_figures)
        b3.grid(row=button_row + 1, column=3, sticky="w")

        menu_row = button_row + 2
        self.add_text0(self.frame, text="Current File:", subscript="", row=menu_row)
        # menu to choose which file to display
        self.current_file = tk.StringVar(value='No files yet')
        self.filemenu = tk.OptionMenu(self.frame, self.current_file, self.current_file.get())
        self.filemenu.grid(row=menu_row, column=1, columnspan=3, sticky='ew')
        self.current_file.trace_add("write", self.generate_plots)

        # menu to choose which measurement (P, D, K, or flux) to display in bottom left graph
        self.current_variable = tk.StringVar(value='Permeability')
        self.add_text0(self.frame, text="Current Measurement:", subscript="", row=menu_row + 1)
        self.PDK_menu = tk.OptionMenu(self.frame, self.current_variable, 'Permeability', 'Diffusivity',
                                      'Solubility', 'Flux')
        self.PDK_menu.grid(row=menu_row + 1, column=1, columnspan=3, sticky='ew')
        self.current_variable.trace_add("write", self.update_PDK_plot)

        b4 = ttk.Button(self.frame, text='Export to Excel', command=self.export_data)
        b4.grid(row=menu_row + 1, column=4, sticky="ew")

        # Store the "set_message" functions
        self.message_function = {}

        # create bottom left plot
        self.ax_title = "Permeability vs. Temperature"
        self.ax_xlabel = " Temperature (\u00B0C) "  # Space at the beginning and end is b/c on Mac, T & e are too close
        self.ax_ylabel = "Permeability (mol m$^{-1}\, $s$^{-1}\, $Pa$^{-0.5}$)"
        self.fig, self.ax, self.canvas, self.toolbar = self.add_plot(self.bottom_frame,
                                                                     xlabel=self.ax_xlabel,
                                                                     ylabel=self.ax_ylabel,
                                                                     title=self.ax_title,
                                                                     row=0, column=0, axes=[.3, .15, .65, .75])
        # create top right plot
        self.ax1_title = "Pressure vs. Time"
        self.ax1_xlabel = "Time (s)"
        self.ax1_ylabel = "Secondary Pressure (Pa)"
        self.ax12_ylabel = "Primary Pressure (Pa)"
        self.fig1, self.ax1, self.canvas1, self.toolbar1 = self.add_plot(self.top_frame,
                                                                         xlabel=self.ax1_xlabel,
                                                                         ylabel=self.ax1_ylabel,
                                                                         title=self.ax1_title,
                                                                         row=0, column=1, rowspan=1)
        self.ax12 = self.ax1.twinx()
        self.ax12.set_ylabel(self.ax12_ylabel)
        # The top right plot would flicker every time the order of magnitude of the cursor's location would change, so
        # the following line changes the format of the displayed x & y coordinates. The specific format chosen was the
        # result of trial and error in conjunction with editing the text of the entry boxes for ss_tol and ss_t_del
        # (which are now accessed by a button). This is kept in case coordinates are changed to on by default again.
        self.ax12.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)

        # Create bottom middle plot
        self.ax2_title = "Permeability vs. Time"
        self.ax2_xlabel = "Time (s)"
        self.ax2_ylabel = r"Permeability (mol m$^{-1}\, $s$^{-1}\, $Pa$^{-0.5}$)"
        self.ax22_ylabel = "Sample Temperature (\u00B0C)"
        self.fig2, self.ax2, self.canvas2, self.toolbar2 = self.add_plot(self.bottom_frame,
                                                                         xlabel=self.ax2_xlabel,
                                                                         ylabel=self.ax2_ylabel,
                                                                         title=self.ax2_title,
                                                                         row=0, column=1, axes=[.15, .15, .7, .75])
        self.ax22 = self.ax2.twinx()
        self.ax22.set_ylabel(self.ax22_ylabel)

        # Create bottom right plot
        self.ax3_title = "Diffusivity Optimization Comparison"
        self.ax3_xlabel = "Time (s)"
        self.ax3_ylabel = r"$(~J_{\mathrm{t}} - J_{\mathrm{0}}~)~/~(~J_{\mathrm{inf}} - J_{\mathrm{0}}~)$"
        self.fig3, self.ax3, self.canvas3, self.toolbar3 = self.add_plot(self.bottom_frame,
                                                                         xlabel=self.ax3_xlabel,
                                                                         ylabel=self.ax3_ylabel,
                                                                         title=self.ax3_title,
                                                                         row=0, column=2, axes=[.15, .15, .75, .75])

        # Turn off coordinates to avoid layout="constrained" causing the plots to shift constantly
        self.ax.format_coord = lambda x, y: ''
        self.ax1.format_coord = lambda x, y: ''
        self.ax12.format_coord = lambda x, y: ''
        self.ax2.format_coord = lambda x, y: ''
        self.ax22.format_coord = lambda x, y: ''
        self.ax3.format_coord = lambda x, y: ''

        # Counteracts a bug of layout="constrained" which causes four plots to be generated but not shown until a
        # Popout Plot button is clicked
        plt.close("all")

    def toggle_coordinates(self):
        """ Toggles between being able to see coordinates while hovering over plots """
        if self.coord_b.config('text')[-1] == 'Enable Coordinates':
            # Turn on coordinates
            self.ax.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            self.ax1.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            self.ax12.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            self.ax2.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            self.ax22.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            self.ax3.format_coord = lambda x, y: "({:3g}, ".format(x) + "{:3g})".format(y)
            # Turn on status bar
            self.toolbar.set_message = self.message_function[self.ax]
            self.toolbar1.set_message = self.message_function[self.ax1]
            self.toolbar2.set_message = self.message_function[self.ax2]
            self.toolbar3.set_message = self.message_function[self.ax3]
            self.coord_b.config(text='Disable Coordinates')  # Toggle the text so next time, it goes to the "else" part
        else:
            # Turn off coordinates
            self.ax.format_coord = lambda x, y: ''
            self.ax1.format_coord = lambda x, y: ''
            self.ax12.format_coord = lambda x, y: ''
            self.ax2.format_coord = lambda x, y: ''
            self.ax22.format_coord = lambda x, y: ''
            self.ax3.format_coord = lambda x, y: ''
            # Turn off status bar
            self.toolbar.set_message("")
            self.toolbar1.set_message("")
            self.toolbar2.set_message("")
            self.toolbar3.set_message("")
            self.toolbar.set_message = lambda s: ""
            self.toolbar1.set_message = lambda s: ""
            self.toolbar2.set_message = lambda s: ""
            self.toolbar3.set_message = lambda s: ""
            self.coord_b.config(text='Enable Coordinates')  # Toggle the text so next time, it goes to the "if" part
        self.canvas.draw()
        self.toolbar.update()

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
            axes are kept in case layout="constrained" has problems. """
        # location of main plot
        frame = tk.Frame(parent)
        frame.grid(row=row, column=column, rowspan=rowspan, sticky="nsew")
        """ Note, the below line is to allow plots to behave better when not at optimal size. See details here:
        https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html
        The above website says the following as of 10/16/2021: 
        "Currently Constrained Layout is experimental. The behaviour and API are subject to change, 
        or the whole functionality may be removed without a deprecation period..."
        As of 9/6/22, according to https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.5.0.html there
        was an update. This source says to use layout="constrained": https://matplotlib.org/stable/api/figure_api.html
        """
        fig, ax = matplotlib.pyplot.subplots(layout="constrained")
        # If something's wrong with layout="constrained", comment out the above line and uncomment the two lines below
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

        # Add a button that allows user input of various steady state variables
        if title == 'Pressure vs. Time':
            entry_frame = tk.Frame(toolbar)
            entry_frame.pack(side="left", padx=1)
            b = tk.Button(entry_frame, text='Steady State Variables', command=self.adjust_ss_vars)
            b.grid(row=0, column=4, sticky="w")

        # Add a button that allows user input of filtering and smoothing variables
        if title == 'Pressure vs. Time':
            entry_frame = tk.Frame(toolbar)
            entry_frame.pack(side="left", padx=1)
            b = tk.Button(entry_frame, text='Filter and Smoothing', command=self.adjust_filter)
            b.grid(row=0, column=5, sticky="w")

        # Store and turn off the capability to update the status bar
        self.message_function[ax] = toolbar.set_message
        toolbar.set_message = lambda s: ""

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
        """ Facilitates loading data and plotting the permeation tab plots """
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
                                          message="Please select a folder with data.")
                return

            self.loading_data = True
            try:
                self.get_persistent_vars()  # Loads variables from a data file for use in analysis
            except Exception as e:
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

            self.datafiles.clear()  # Clear the dictionary to avoid old data points appearing in PDK plot
            self.error_texts = ""  # Reset error_texts each time a folder is loaded

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

                # Process the data by making the various calculations
                if status:
                    # Reset these variables for each loadable file
                    self.time_type = 0
                    self.days = -1
                    self.extra_time = 0

                    try:  # Catch most other errors that can happen while processing the data
                        self.datafiles[filename] = self.extract_data(os.path.join(self.directory, filename))

                        # Find row number where the Isolation Valve is 1st opened after being closed (aka, the start of
                        # the experiment). Done in a try/except clause in case no opening is registered
                        try:
                            _t0 = np.where((self.datafiles[filename])['Isolation Valve'] == 0)[0]
                            closed = np.where(np.diff(_t0) > 1)[0]
                            if not closed.any():  # if the isolation valve wasn't closed again after being opened...
                                # ...set t0 to be the data point after the last zero in the file
                                self.t0[filename] = _t0[-1] + 1
                            else:  # else, set t0 to be the data point after the last zero before the valve is opened
                                self.t0[filename] = _t0[closed[0]] + 1
                        except IndexError:
                            self.error_texts += "Loading Error with file " + filename + \
                                                ". Incorrect Isolation Valve format.\n\n"
                            self.datafiles.pop(filename)
                            status = False

                        if status:  # Ensure status hasn't changed since original definition
                            self.calculate_permeability(filename, self.datafiles[filename])
                            self.calculate_diffusivity(filename, self.datafiles[filename])
                            self.calculate_solubility(filename, self.datafiles[filename])
                    except Exception as e:
                        # If you want to see the traceback as part of the error text, change the below value to true
                        want_traceback = True
                        if want_traceback:
                            import traceback
                            full_traceback = "\n" + traceback.format_exc()
                        else:
                            full_traceback = ""

                        self.error_texts += "Unknown Error with file " + filename + '. The following exception was' + \
                                            ' raised: "' + str(e) + '".' + full_traceback + '\n\n'
                        try:
                            self.datafiles.pop(filename)
                        except KeyError:
                            self.error_texts += "Note: " + filename + ' was never successfully loaded.' + '\n\n'
                        status = False

                if not status:  # If something is wrong with the file or with analyzing the file
                    self.options.pop(self.options.index(filename))
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
        """ Update the dataframe that is used for filling the Excel sheet and plotting the data in overview plots """
        # this will recreate the dataframe each time
        df = pd.DataFrame()
        for filename in self.options:
            df = pd.concat([df, pd.DataFrame(
                {"Gas Temperature [K]": self.Tgss[filename],
                 "Gas Temperature Uncertainty [K]": self.Tgss_err[filename],
                 "Sample Temperature [K]": self.Tsampss[filename],
                 "Sample Temperature Uncertainty [K]": self.Tsampss_err[filename],
                 "Secondary Side Pressure [Pa]": self.SecP_average[filename],
                 "Secondary Side Pressure Uncertainty [Pa]": self.SecP_err[filename],
                 "Primary Side Pressure [Pa]": self.PrimP_average[filename],
                 "Primary Side Pressure Uncertainty [Pa]": self.PrimP_err[filename],
                 "Permeability [mol m^-1 s^-1 Pa^-0.5]": self.Phi[filename],
                 "Permeability Uncertainty [mol m^-1 s^-1 Pa^-0.5]": self.Phi_err[filename],
                 "Diffusivity [m^2 s^-1]": self.D[filename],
                 "Diffusivity Uncertainty [m^2 s^-1]": self.D_err[filename],
                 "Diffusivity (timelag) [m^2 s^-1]": self.D_tlag[filename],
                 "Diffusivity (timelag) Uncertainty [m^2 s^-1]": self.D_tlag_err[filename],
                 "Solubility [mol m^-3 Pa^-0.5]": self.Ks[filename],
                 "Solubility Uncertainty [mol m^-3 Pa^-0.5]": self.Ks_err[filename],
                 "Flux [mol m^-2 s^-1]": self.F[filename],
                 "Flux Uncertainty [mol m^-2 s^-1]": self.F_err[filename],
                 "Proportionality Constant A": self.A[filename],
                 "Additive Time Constant dt [s]": self.dt[filename]
                 }, index=[filename]
                )])
        self.storage.TransportParameters = df
        self.storage.PTransportParameters = df

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
                self.storage.PTransportParameters.to_excel(save_filename)

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
                            " into seconds: " + str(value) + ". This triggered a ValueError. The dataframe " + \
                            " location corresponding to this term has been set to NaN. No further information" \
                            " is available.\n\n"
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
                                " into seconds: " + str(time_term) + ". This triggered a ValueError. The dataframe " + \
                                " location corresponding to this term has been set to NaN. No further information" \
                                " is available.\n\n"
                return float("NaN")
        if self.time_type == 1:
            # If the time is in cumulative time passed format, do a linear transform according to multi and add
            return self.unit_conversion(time_term, multi, add)
        elif self.time_type == 2:
            # If the time is in datetime format, convert to cumulative seconds passed
            return self.get_seconds(time_term)

    def get_persistent_vars(self):
        """ Read data from persistent_permeation_input_variables.xlsx. This data is critical to processing the
            data files loaded in by the user. """
        # Open up the persistent variable file for reading
        pv_filename = os.path.join('data_files', 'persistent_permeation_input_variables.xlsx')
        pv_wb = openpyxl.load_workbook(pv_filename)

        self.num_GasT = pv_wb['Numbers']['C2'].value  # Number of TCs measuring GasT

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

        # Primary pressure
        PrimP_info = pd.read_excel(pv_filename, sheet_name="PrimP", header=0)
        for PrimP_inst in PrimP_info.keys():
            self.header.append(PrimP_inst)
            self.needed_cols.append(self.l2n(PrimP_info[PrimP_inst][0]))
            self.converters_dict[self.needed_cols[-1]] = \
                lambda input, multi=PrimP_info[PrimP_inst][1], add=PrimP_info[PrimP_inst][2]: \
                self.unit_conversion(input, multi, add)
            self.PrimP_cerr[PrimP_inst] = PrimP_info[PrimP_inst][3]
            self.PrimP_perr[PrimP_inst] = PrimP_info[PrimP_inst][4]

        # Secondary pressure
        SecP_info = pd.read_excel(pv_filename, sheet_name="SecP", header=0)
        for SecP_inst in SecP_info.keys():
            self.header.append(SecP_inst)
            self.needed_cols.append(self.l2n(SecP_info[SecP_inst][0]))
            self.converters_dict[self.needed_cols[-1]] = \
                lambda input, multi=SecP_info[SecP_inst][1], add=SecP_info[SecP_inst][2]: \
                self.unit_conversion(input, multi, add)
            self.SecP_cerr[SecP_inst] = SecP_info[SecP_inst][3]
            self.SecP_perr[SecP_inst] = SecP_info[SecP_inst][4]

        # Isolation Valve shouldn't need more than one instrument, so don't need a for loop. Also don't need uncertainty
        self.header.append("Isolation Valve")
        self.needed_cols.append(self.l2n(pv_wb['MiscInfo']['B2'].value))
        self.converters_dict[self.needed_cols[-1]] = \
            lambda input, multi=pv_wb['MiscInfo']['B3'].value, add=pv_wb['MiscInfo']['B4'].value: \
            round(self.unit_conversion(input, multi, add))

        # Variable for which row in Excel sheet to start obtaining data from (0-indexed)
        self.starting_row = pv_wb['MiscInfo']['C2'].value
        # Variable for how many rows at the end of the Excel sheet to not read (0-indexed)
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
                                   skipfooter=self.footer_rows, engine='python')
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

            # Now that the errors are handled, extract the data from the Excel file and load it into a dataframe
            if self.file_type == ".xls":
                Data = pd.read_csv(filename, sep='\t', header=None, usecols=self.needed_cols,
                                   converters=self.converters_dict, skiprows=self.starting_row,
                                   skipfooter=self.footer_rows, engine='python')
            else:  # assume .xlsx file
                Data = pd.read_excel(filename, header=None, engine="openpyxl", usecols=self.needed_cols,
                                     converters=self.converters_dict, skiprows=self.starting_row,
                                     skipfooter=self.footer_rows)
            # https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
            new_header = [x for _, x in sorted(zip(self.needed_cols, self.header))]
            Data.columns = new_header

        #update the time in column A
        t_0 = Data.loc[0,'t']
        tOld = np.array(Data['t'])
        tnew = tOld - t_0
        new_time = pd.Series(tnew, name='t')
        Data.update(new_time)


        if self.filter==True: #only turn on the filter if the switch is on
            y = np.array(Data['SecP']) #turn the secondary pressure data into an array
            window = self.filter_win.get() #get the window size (default is 20)
            degree = self.poly_deg.get() #get the polynomial degree size (default is 2)
            SecP_filtered = pd.Series(savgol_filter(y, window, degree), name='SecP') # Savitzky-Golay filter
            Data.update(SecP_filtered) #update the data frame with the filtered pressures

        n = len(Data)

        # calculate numerical derivative of secondary pressure
        deriv = np.zeros(n)
        for i in range(2):
            #the first two points don't matter, hold their place with just a basic right hand derivative
            deriv[i] = ((Data.loc[i + 1, 'SecP'] + Data.loc[i, 'SecP']) /
                        (Data.loc[i + 1, 't'] - Data.loc[i, 't']))
        for i in range(2, n - 2):
            #use the 5-point stencil to calculate the derivative
            deriv[i] = ((-Data.loc[i + 2, 'SecP'] + 8*Data.loc[i + 1, 'SecP'] - 8*Data.loc[i - 1, 'SecP'] + Data.loc[i - 2, 'SecP']) /
                        12*(Data.loc[i + 1, 't'] - Data.loc[i, 't']))


        Data['dSecP'] = deriv.tolist()

        # rearrange last two columns so 'SecP' and 'dSecP' are adjacent
        cols = Data.columns.tolist()
        cols = cols[:-2] + [cols[-1]] + [cols[-2]]
        Data = Data[cols]

        return Data

    def adjust_persistent_vars(self, retry=False):
        """ Function for calling the window for adjusting persistent variables for loading permeation data,
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
            self.get_persistent_vars()  # Loads data from the permeation input file-format file

    def adjust_filter(self):
        """ Function for calling the window for adjusting filtering variables for loading permeation data,
            then update everything if desired """
        popup = tk.Toplevel(self)
        popup.wm_title("Adjust filter Variables")

        # Place the popup window near the cursor
        pos_right = self.winfo_pointerx()
        pos_down = self.winfo_pointery()
        popup.geometry("+{}+{}".format(pos_right, pos_down))

        entry_frame = tk.Frame(popup)
        entry_frame.pack(side="left", padx=1)

        #Make a toggle switch for turning on and off the rolling average
        def Rolling_toggle():
            if rolling_button.config('text')[-1] == 'ON':
                rolling_button.config(text='OFF')
                self.rolling = False
            else:
                rolling_button.config(text='ON')
                self.rolling = True

        #put the rolling average toggle switch on the frame
        rolling_button = ttk.Button(entry_frame, text="OFF", command=Rolling_toggle)
        rolling_button.grid(row=1, column=1, sticky='ew')
        savgol_label = ttk.Label(entry_frame, text='Rolling Average')
        savgol_label.grid(row=1, column=0)

        #Make a toggle switch for filtering
        def Savgol_toggle():
            if savgol_toggle_button.config('text')[-1] == 'ON':
                savgol_toggle_button.config(text='OFF')
                self.filter = False
            else:
                savgol_toggle_button.config(text='ON')
                self.filter = True

        #put the filter toggel switch on the page
        savgol_toggle_button = ttk.Button(entry_frame, text="ON", command=Savgol_toggle)
        savgol_toggle_button.grid(row=2, column=1, sticky='ew')
        savgol_label = ttk.Label(entry_frame, text='Savitzky-Golay Filter')
        savgol_label.grid(row=2, column=0)

        # Degree of polynomial
        self.add_entry(popup, entry_frame, self.inputs, key="poly_deg",
                       text="Polynomial Number:", subscript='', units="",
                       tvar=self.poly_deg, ent_w=8, row=3, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))
        # Filter window size
        self.add_entry(popup, entry_frame, self.inputs, key="filter_win",
                       text="Filter Window Size:", subscript='', units="",
                       tvar=self.filter_win, ent_w=8, row=4, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))
        
        button_row = 7
        b1 = ttk.Button(entry_frame, text='Close & Refresh', command=lambda: self.close_and_refresh(popup))
        b1.grid(row=button_row, column=0, sticky="ew")
        b2 = ttk.Button(entry_frame, text='Refresh', command=self.refresh_graphs)
        b2.grid(row=button_row, column=2, sticky="ew")
        tk.mainloop()

    def adjust_ss_vars(self):
        """ Function for calling the window for adjusting steady state variables for loading permeation data,
            then update everything if desired """
        popup = tk.Toplevel(self)
        popup.wm_title("Adjust Steady State Variables")

        # Place the popup window near the cursor
        pos_right = self.winfo_pointerx()
        pos_down = self.winfo_pointery()
        popup.geometry("+{}+{}".format(pos_right, pos_down))

        entry_frame = tk.Frame(popup)
        entry_frame.pack(side="left", padx=1)
        # Tolerance for determining steady state
        self.add_entry(popup, entry_frame, self.inputs, key="ss_tol",
                       text="Tolerance for Steady State:", subscript='', units="Pa s\u207b\u00b2",
                       tvar=self.ss_tol, formatting=True, ent_w=8, in_window=True,
                       command=lambda tvar, variable, key, formatting, pf:
                       self.storage.check_for_number(tvar, variable, key, formatting, pf))
        # Minimum time to let pass after t0 before checking for steady state
        self.add_entry(popup, entry_frame, self.inputs, key="ss_t_del",
                       text="Delay until Steady State:", subscript='', units="s",
                       tvar=self.ss_t_del, ent_w=8, row=1, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))  # "false" to say no formatting
        # Number of data points used in determining the leak
        self.add_entry(popup, entry_frame, self.inputs, key="leak_range",
                       text="Leak Range:", subscript='', units="data points",
                       tvar=self.gen_leak_range, ent_w=8, row=3, in_window=True,
                       command=lambda tvar, variable, key, pf:
                       self.storage.check_for_number(tvar, variable, key, False, pf))
        # Number of data points used in determining the steady state and max time
        self.add_entry(popup, entry_frame, self.inputs, key="ss_range",
                    text="Steady State Range:", subscript='', units="data points",
                    tvar=self.gen_ss_range, ent_w=8, row=4, in_window=True,
                    command=lambda tvar, variable, key, pf:
                    self.storage.check_for_number(tvar, variable, key, False, pf))

        button_row = 5
        b1 = ttk.Button(entry_frame, text='Close & Refresh', command=lambda: self.close_and_refresh(popup))
        b1.grid(row=button_row, column=0, sticky="ew")
        b2 = ttk.Button(entry_frame, text='Refresh', command=self.refresh_graphs)
        b2.grid(row=button_row, column=2, sticky="ew")

    def close_and_refresh(self, win):
        """ Accepts window argument, then closes that window after refreshing the graphs of the Permeation Plots tab """
        self.refresh_graphs()
        win.destroy()

    def get_rates(self, X, Y):
        """ Calculate rate of change using linear algebra. Calculations taken from
        https://stackoverflow.com/questions/22381497/python-scikit-learn-linear-model-parameter-standard-error"""

        X = pd.DataFrame(X.values, columns=['X'])
        y = pd.DataFrame(Y.values, columns=['Y'])
        # reproduce with linear algebra - will probably just remove scikit after :P
        N = len(X)
        p = len(X.columns) + 1  # plus one because LinearRegression adds an intercept term

        X_with_intercept = np.empty(shape=(N, p), dtype='float')
        X_with_intercept[:, 0] = 1
        X_with_intercept[:, 1:p] = X.values

        beta_hat = np.linalg.inv(X_with_intercept.T @ X_with_intercept) @ X_with_intercept.T @ y.values
        # order is intercept, then coefficients
        m = beta_hat[1][0]
        b = beta_hat[0][0]

        # compute the standard error for the parameter estimates
        y_hat = np.array(m * X + b)

        residuals = y.values - y_hat
        residual_sum_of_squares = residuals.T @ residuals
        sigma_squared_hat = residual_sum_of_squares[0, 0] / (N - p)
        var_beta_hat = np.linalg.inv(X_with_intercept.T @ X_with_intercept) * sigma_squared_hat
        err = []
        for p_ in range(p):
            standard_error = var_beta_hat[p_, p_] ** 0.5
            err.append(standard_error)

        return {"slope": m, "slope error": err[1], "intercept": b, "intercept error": err[0]}

    def calculate_permeability(self, filename, data):
        """ Calculate Permeability using loaded data """
        R = self.storage.R * 1000  # J/mol K

        # Set the leak range to the general leak range (which is editable by the user) or to something that works better
        if self.gen_leak_range.get() < self.t0[filename]:
            self.leak_range[filename] = self.gen_leak_range.get()
        else:
            self.leak_range[filename] = self.t0[filename] - 1
            self.error_texts += "Leak Range Warning with file " + filename + ". The user-input leak range exceeded" + \
                                " the limit of " + str(self.t0[filename]) + " data points. Leak range set to " + \
                                str(self.leak_range[filename]) + " data points.\n\n"
        # Set the ss range to the general ss range (which is editable by the user) or to something that works better
        if self.gen_ss_range.get() < len(data['dSecP']) - self.t0[filename]:
            self.ss_range[filename] = self.gen_ss_range.get()
        else:
            self.ss_range[filename] = len(data['dSecP']) - self.t0[filename] - 2
            self.error_texts += "Steady State Warning with file " + filename + ". The user-input steady state" + \
                                " range exceeded the limit of " + str(len(data['dSecP']) - self.t0[filename]) + \
                                " data points. Steady state range set to " + str(self.ss_range[filename]) + \
                                " data points.\n\n"
            
        # determine the leak rate. self.leak_range points before opening, not including t0 (note that .loc is inclusive)
        self.pleak_lf[filename] = \
            self.get_rates(data.loc[self.t0[filename] - self.leak_range[filename]:self.t0[filename] - 1, 't'],
                           data.loc[self.t0[filename] - self.leak_range[filename]:self.t0[filename] - 1, 'SecP'])
        
        # determine row number when steady state pressure rate is achieved, checking to ensure it is before the end of
        # the file and after the minimum time delay + t0
        if self.rolling==False:
            dSecP = data['dSecP']
        else:
            dSecP = data['dSecP'].rolling(window=self.ss_range[filename], center=True, min_periods=1).mean()
        neg_dSecP = np.where(dSecP < 0)[0]  # List of all terms in dSecP which are less than 0
        # Minimum terms in a row required to determine if the pressure drop off was reached
        min_seq = max(self.ss_range[filename] // 2 - 1, 1)
        ss_time_max = len(dSecP)  # variable to store the last point before experiment ends (and/or pressure drops)
        #  Find where the pressure drops off because the valve was opened or fail to find such a point
        for count, value in enumerate(neg_dSecP[:len(neg_dSecP)-min_seq]):
            if value > self.t0[filename]:  # if the value is late enough
                # Check for a sequence min_seq long in which neg_dSecP is less than zero
                if sum(np.diff(neg_dSecP[count:count+min_seq])) <= min_seq:
                    ss_time_max = neg_dSecP[count]  # set the max time to when the pressure drops
                    break

        # Convert the user's input time delay into a location in the array of data that doesn't depend on delta t
        # This line creates an array containing only the locations of times that are eligible to be a steady state
        ss_del = np.where(data.loc[self.t0[filename]:, 't'] > self.ss_t_del.get() + data.loc[self.t0[filename], 't']
                          + 0.000001)[0]  # Added .00001 to ensure precision errors don't lead to e_del being 0
        try:
            # Set the time delay to be the minimum number of data points away from t0 that still is further than
            # the user's input time delay
            ss_del = ss_del[0] - 1  # -1 to compensate for the > sign in e_del's def and in later checks
        except IndexError:  # This is expected if ss_del is empty
            # Set ss_del to the largest delay that still works
            ss_del = ss_time_max - (self.ss_range[filename] + self.t0[filename] + 2)
            self.error_texts += "Steady State Warning with file " + filename + ". The user-input time delay" + \
                                " exceeded the limit of " + \
                                str(round(data.loc[ss_time_max - (self.ss_range[filename] + 2), 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s. Time delay set to " + \
                                str(round(data.loc[ss_del + self.t0[filename], 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s.\n\n"
        if ss_del > ss_time_max - (self.ss_range[filename] + self.t0[filename] + 2):
            # In case ss_del produces a usable data point but that data point is past an abrupt downturn in pressure
            ss_del = ss_time_max - (self.ss_range[filename] + self.t0[filename] + 2)
            self.error_texts += "Steady State Warning with file " + filename + ". The user-input time delay" + \
                                " exceeded the limit of " + \
                                str(round(data.loc[ss_time_max - (self.ss_range[filename] + 2), 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s. Time delay set to " + \
                                str(round(data.loc[ss_del + self.t0[filename], 't'] -
                                          data.loc[self.t0[filename], 't'], 3)) + " s.\n\n"

        # Find the data point at which an equilibrium is reached
        ddSecP = pd.Series((-dSecP.loc[4:].to_numpy() + 8*dSecP.loc[3:len(dSecP)-2].to_numpy() - 8*dSecP.loc[1:len(dSecP)-4].to_numpy() + dSecP.loc[:len(dSecP)-5].to_numpy()) /
                           12*(data.loc[4:, 't'].to_numpy() - data.loc[3:len(dSecP)-2, 't'].to_numpy()))
        ddSecP = ddSecP.rolling(window=self.ss_range[filename], center=True, min_periods=1).mean()
        zeros = np.where(abs(ddSecP) < self.ss_tol.get())[0]
        zeros = [z for z in zeros if self.t0[filename] + ss_del < z < ss_time_max - min_seq - 2]
        # If no steady state was found with the user-input tolerance, loop until find a steady state using new_ss_tol
        new_ss_tol = self.ss_tol.get()
        while not zeros and new_ss_tol < self.ss_tol.get() * 100000:
            new_ss_tol = float("{:.2e}".format(new_ss_tol * 5))
            zeros = np.where(abs(ddSecP) < new_ss_tol)[0]
            zeros = [z for z in zeros if self.t0[filename] + ss_del < z < ss_time_max - min_seq - 2]
        if new_ss_tol != self.ss_tol.get() and zeros:
            self.error_texts += "Steady State Warning with file " + filename + ". HyPAT was unable to find a steady" + \
                                " state using the user-input tolerance " + str(self.ss_tol.get()) + \
                                ". New tolerance was set to {:.2e}".format(new_ss_tol) + \
                                ", which successfully determined a time for the beginning of the steady state.\n\n"
        elif new_ss_tol != self.ss_tol.get():
            zeros = [min(self.t0[filename] + ss_del + 1, ss_time_max - min_seq - 3)]
            self.error_texts += "Steady State Error with file " + filename + \
                                ". HyPAT was unable to find steady state using the user-input tolerance " + \
                                str(self.ss_tol.get()) + " or the computationally set tolerance " + \
                                "{:.2e}".format(new_ss_tol) + ". Steady state time was set to " + \
                                str(round(data.loc[zeros[0], 't'], 3)) + \
                                " s, which is equal to either the starting time plus the user-input" + \
                                " time delay until steady state or the maximum usable time, whichever is lower.\n\n"
        # Set steady state time
        self.tss[filename] = zeros[0]

        # determine steady state permeation rates
        self.pss_lf[filename] = self.get_rates(data.loc[self.tss[filename] + 1:
                                                        self.tss[filename] + self.ss_range[filename], 't'],
                                               data.loc[self.tss[filename] + 1:
                                                        self.tss[filename] + self.ss_range[filename], 'SecP'])
        self.Prate[filename] = self.pss_lf[filename]["slope"] - self.pleak_lf[filename]["slope"]  # Pa/s

        # Get the temperature using an average of arbitrarily many TCs
        Tg_vals = pd.DataFrame()
        for i in range(1, self.num_GasT + 1):
            Tg_vals['GasT' + str(i)] = data.loc[:, 'GasT' + str(i)]
        self.Tg_mean = Tg_vals.mean(axis=1)
        # Average initial temperature and average steady state temperature
        self.Tg0[filename] = self.Tg_mean.loc[self.t0[filename] - self.leak_range[filename]:
                                              self.t0[filename] - 1].mean() + self.storage.standard_temp  # K
        self.Tgss[filename] = \
            self.Tg_mean.loc[self.tss[filename] + 1:self.tss[filename] + self.ss_range[filename]].mean() + \
            self.storage.standard_temp  # K

        # convert to molar flow rate mol/s using ideal gas law
        Nrate = (self.pss_lf[filename]["slope"] * self.volume.get()) / (R * self.Tgss[filename]) - \
                (self.pleak_lf[filename]["slope"] * self.volume.get()) / (R * self.Tg0[filename])
        # convert to molar flux
        self.F[filename] = Nrate / self.A_perm.get()

        # calculate permeability
        self.SecP_average[filename] = data.loc[self.tss[filename] + 1:
                                               self.tss[filename] + self.ss_range[filename], 'SecP'].mean()
        self.PrimP_average[filename] = data.loc[self.tss[filename] + 1:
                                                self.tss[filename] + self.ss_range[filename], 'PrimP'].mean()

        sl = self.sample_thickness.get() * 0.001  # converted to meters
        sl_err = self.sample_thickness_err.get() * 0.001  # converted to meters
        self.Phi[filename] = self.F[filename] * sl / (np.sqrt(self.PrimP_average[filename]) -
                                                      np.sqrt(self.SecP_average[filename]))

        # Calculate the uncertainty in permeability
        # We need uncertainty from the following: sl, SecP, PrimP, A_perm, pss_lf[slope], pleak_lf[slope], V, Tgss, Tg0

        # Find uncertainty of Tg0 when there are arbitrarily many TCs
        Tg0_errs = []
        for i in range(1, self.num_GasT + 1):
            Tg0_errs.append((max(data.loc[self.t0[filename] - self.leak_range[filename]:self.t0[filename] - 1,
                                 'GasT' + str(i)].std(), self.GasT_cerr['GasT' + str(i)],
                                 (self.Tg0[filename] - self.storage.standard_temp) *
                                 self.GasT_perr['GasT' + str(i)] / 100)) ** 2)
        Tg0_err = np.sqrt(sum(Tg0_errs)) / self.num_GasT
        # Find uncertainty of Tgss when there are arbitrarily many TCs
        Tgss_errs = []
        for i in range(1, self.num_GasT + 1):
            Tgss_errs.append((max(data.loc[self.tss[filename] + 1:self.tss[filename] + self.ss_range[filename],
                                  'GasT' + str(i)].std(), self.GasT_cerr['GasT' + str(i)],
                                  (self.Tgss[filename] - self.storage.standard_temp) *
                                  self.GasT_perr['GasT' + str(i)] / 100)) ** 2)
        Tgss_err = np.sqrt(sum(Tgss_errs)) / self.num_GasT

        # Pressure uncertainty
        SecP_err = max(data.loc[self.tss[filename] + 1:self.tss[filename] + self.ss_range[filename], 'SecP'].std(),
                       self.SecP_cerr["SecP"], self.SecP_average[filename] * self.SecP_perr["SecP"] / 100)
        PrimP_err = max(data.loc[self.tss[filename] + 1:self.tss[filename] + self.ss_range[filename], 'PrimP'].std(),
                        self.PrimP_cerr["PrimP"], self.PrimP_average[filename] * self.PrimP_perr["PrimP"] / 100)

        # Rate uncertainty
        Prate_err = np.sqrt(self.pss_lf[filename]['slope error'] ** 2 + self.pleak_lf[filename]['slope error'] ** 2)
        Nrate_err = 1 / R * np.sqrt((abs(self.pss_lf[filename]['slope'] * self.volume.get() / self.Tgss[filename]) *
                                     np.sqrt((self.pss_lf[filename]['slope error'] / self.pss_lf[filename]['slope']) ** 2 +
                                             (Tgss_err / self.Tgss[filename]) ** 2 +
                                             (self.volume_err.get() / self.volume.get()) ** 2)
                                     ) ** 2 +
                                    (abs(self.pleak_lf[filename]['slope'] * self.volume.get() / self.Tg0[filename]) *
                                     np.sqrt((self.pleak_lf[filename]['slope error'] / self.pleak_lf[filename]['slope']) ** 2 +
                                             (Tg0_err / self.Tg0[filename]) ** 2 +
                                             (self.volume_err.get() / self.volume.get()) ** 2)
                                     ) ** 2
                                    )  # self.volume could get pulled out front, but hasn't been to ease comprehension

        # Final results
        F_err = abs(Nrate/self.A_perm.get())*np.sqrt((Nrate_err/Nrate) ** 2 +
                                                     (self.A_perm_err.get() / self.A_perm.get()) ** 2)
        Phi_err = self.Phi[filename] * np.sqrt((F_err / self.F[filename]) ** 2 + (sl_err/sl) ** 2 +
                                               (np.sqrt(
                                                   (SecP_err / (2 * np.sqrt(self.SecP_average[filename]))) ** 2 +
                                                   (PrimP_err / (2 * np.sqrt(self.PrimP_average[filename]))) ** 2) /
                                                        (np.sqrt(self.PrimP_average[filename]) - np.sqrt(
                                                            self.SecP_average[filename]))
                                                ) ** 2
                                               )

        # store uncertainties
        self.Prate_err[filename] = Prate_err
        self.Phi_err[filename] = Phi_err
        self.F_err[filename] = F_err
        self.Tgss_err[filename] = Tgss_err
        self.SecP_err[filename] = SecP_err
        self.PrimP_err[filename] = PrimP_err

        # Data for permeability vs. time graph
        self.Perm_Plot[filename] = ((dSecP - self.pleak_lf[filename]["slope"]) * self.volume.get() * sl) / \
                                   (((data['PrimP']) ** 0.5 - (
                                                        data['SecP']) ** 0.5) * R * (
                                                        self.Tg_mean + self.storage.standard_temp) * self.A_perm.get())

        # get temp from sample TC
        self.Tsamp0[filename] = data.loc[self.t0[filename] - self.leak_range[filename]:
                                         self.t0[filename] - 1, 'SampT'].mean() + self.storage.standard_temp
        self.Tsampss[filename] = data.loc[self.tss[filename] + 1:
                                          self.tss[filename] + self.ss_range[filename], 'SampT'].mean() + \
            self.storage.standard_temp
        self.Tsampss_err[filename] = max(data.loc[self.tss[filename] + 1:
                                                  self.tss[filename] + self.ss_range[filename], 'SampT'].std(),
                                         self.SampT_cerr["SampT"],
                                         self.Tsampss[filename] * self.SampT_perr["SampT"] / 100)

    def calculate_diffusivity(self, filename, data):
        """ Calculate Diffusivity using loaded data """
        debug = False
        R = self.storage.R * 1000  # J/mol K

        # Prep for the time-lag method
        x_intercept = (self.pss_lf[filename]["intercept"] - self.pleak_lf[filename]["intercept"]) / \
                      (self.pleak_lf[filename]["slope"] - self.pss_lf[filename]["slope"])
        
        #find error from x intercept
        ss_slope_err = self.pss_lf[filename]['slope error']
        ss_int_err = self.pss_lf[filename]['intercept error']
        leak_slope_err = self.pleak_lf[filename]['slope error']
        leak_int_err = self.pleak_lf[filename]['intercept error']
        dif1 = self.pleak_lf[filename]["slope"]-self.pss_lf[filename]["slope"]
        dif2 = self.pss_lf[filename]["intercept"] - self.pleak_lf[filename]["intercept"]
        x_err = np.sqrt((ss_int_err / dif1)**2 + 
                        (leak_int_err / dif1)**2 + 
                        (leak_slope_err * dif2 / dif1**2)**2 + 
                        (ss_slope_err * dif2 / dif1**2)**2)
        
        #find the y intercept as well as the time lag
        y_intercept = self.pleak_lf[filename]["slope"] * x_intercept + self.pleak_lf[filename]["intercept"]
        self.intercept[filename] = (x_intercept, y_intercept)
        self.tlag[filename] = x_intercept - data.loc[self.t0[filename], 't']
        if self.tlag[filename] < 0:
            self.error_texts += "Time-Lag Error with file " + filename + ". A negative time-lag was calculated. " + \
                                "Diffusivity will not be calculated for this file.\n\n"
            skipD = True
        else:
            skipD = False

        # calculate diffusivity using time-lag method
        sl = self.sample_thickness.get() * 0.001  # converted to meters
        D = sl ** 2 / (6 * self.tlag[filename])
        self.D_tlag[filename] = D
        if skipD:  # If skipping the calculations of D
            D = float("NaN")

        #calculate the uncetainty of D using the time lag method
        #note: the error in the time lag is the same as the error in the x intercept
        sl_err = self.sample_thickness_err.get() * 0.001  # converted to meters
        self.D_tlag_err[filename] = np.sqrt((sl_err * sl / 3 / self.tlag[filename])**2 + 
                                  (x_err * sl**2 / 6 / self.tlag[filename]**2)**2)

        # Prepare for optimization
        Jt = data.loc[self.t0[filename] + 1:self.tss[filename] + self.ss_range[filename], 'dSecP'].to_numpy() * \
            self.volume.get() / \
            (R * (self.Tg_mean.loc[self.t0[filename] + 1:self.tss[filename] + self.ss_range[filename]].to_numpy() +
                  self.storage.standard_temp) * self.A_perm.get())
        self.D_time[filename] = data.loc[self.t0[filename] + 1:self.tss[filename] + self.ss_range[filename], 't'] - \
            data.loc[self.t0[filename], 't']
        J0 = np.mean(data.loc[self.t0[filename] - self.leak_range[filename]:self.t0[filename] - 1, 'dSecP']) * \
            self.volume.get() / (R * self.Tg0[filename] * self.A_perm.get())  # leak rate
        Jinf = np.mean(data.loc[self.tss[filename] + 1:
                                self.tss[filename] + self.ss_range[filename], 'dSecP'].to_numpy()) * \
            self.volume.get() / (R * self.Tgss[filename] * self.A_perm.get())  # steady state flow rate

        if debug:
            # checking some numbers
            print('\n', filename)
            print("mean:", np.mean(data.loc[self.t0[filename] - self.leak_range[filename]:
                                            self.t0[filename] - 1, 'dSecP']))
            print("leak rate:", self.pleak_lf[filename]["slope"])
            print("max:", np.max(Jt))  # large in 1505
            print("Steady State Rate:", self.pss_lf[filename]["slope"])
            print("Diffusivity guess:", D)

        # before optimizing
        self.lhs[filename] = (Jt - J0) / (Jinf - J0)
        self.rhs[filename] = 1 + 2 * sum(
            [(-1) ** n * np.exp(-D * n ** 2 * np.pi ** 2 * self.D_time[filename] / sl ** 2) for n in range(1, 20)])
        if debug:
            plt.figure()
            plt.plot(self.lhs[filename][1:], label='lhs')
            plt.plot(self.rhs[filename][1:], label='rhs')
            plt.title(filename+" before optimizing")
            plt.legend()

            plt.figure(filename)
            plt.plot(self.rhs[filename][1:], label='no optimization')
            plt.plot(self.lhs[filename][1:], label='lhs')
            plt.show()

        # optimizing using curve fit - should just be a wrapper around least squares

        h = data.loc[1, 't'] - data.loc[0, 't']

        def f(xdata, D_opt, dt_opt, A_opt):
            # Function for optimizing D using curve fit
            rhs = 1 + 2 * sum(
                [(-1) ** n * np.exp(-D_opt * n ** 2 * np.pi ** 2 * (xdata + dt_opt) / sl ** 2) for n in range(1, 20)])
            return rhs * A_opt
        # Attempt to optimize D using curve fit and the above function. Show a warning if it fails due to RuntimeError
        try:
            if not skipD:  # If not skipping the calculation of D
                # todo look into making this a weighted curve fit
                popt, pcov = curve_fit(f, self.D_time[filename], self.lhs[filename], p0=[D, 0, 1], xtol=D * 1e-3,
                                       bounds=([0, -min(self.D_time[filename]), -1000], [10, 10*h, 1000]))
            else:
                NaN = float("NaN")
                popt = [NaN, NaN, NaN]
                pcov = np.array([[NaN, NaN, NaN], [NaN, NaN, NaN], [NaN, NaN, NaN]])
        except RuntimeError:
            self.error_texts += "Curve Fit Error with file " + filename + \
                                ". Curve fit unable to find optimal diffusivity parameters.\n\n"
            NaN = float("NaN")
            popt = [D, 0, 1]
            pcov = np.array([[NaN, NaN, NaN], [NaN, NaN, NaN], [NaN, NaN, NaN]])
        except OptimizeWarning:
            self.error_texts += "Curve Fit Warning with file " + filename + \
                                ". Curve fit unable to find covariance of the diffusivity parameters. D has been" + \
                                " set to the diffusivity calculated by the time-lag method.\n\n"
            NaN = float("NaN")
            popt = [D, 0, 1]
            pcov = np.array([[NaN, NaN, NaN], [NaN, NaN, NaN], [NaN, NaN, NaN]])
        perr = np.sqrt(np.diag(pcov))
        self.D_err[filename] = perr[0]
        # Store D, dt, A, and rhs_cf for use elsewhere in program
        self.D[filename], self.dt[filename], self.A[filename] = popt
        self.rhs_cf[filename] = self.A[filename] * (1 + 2 * sum(
            [(-1) ** n * np.exp(-self.D[filename] * n ** 2 * np.pi ** 2 * (self.D_time[filename] + self.dt[filename]) /
                                sl ** 2) for n in range(1, 20)]))
        if debug:
            lhs = self.A[filename] * (Jt - J0) / (Jinf - J0)
            plt.figure()
            plt.plot(lhs[1:], label='lhs')
            plt.plot(self.rhs_cf[filename][1:], label='rhs')
            plt.title(filename + " curve fit")
            plt.legend()

            # compare the optimized lines
            plt.figure(filename)
            plt.plot(self.rhs_cf[filename][1:], '--', label='curve fit')
            plt.title("comparison")
            plt.legend()

    def calculate_solubility(self, filename, data):
        """ Calculate Solubility using calculated Diffusivity and Permeability """
        self.Ks[filename] = self.Phi[filename]/self.D[filename]
        # get the uncertainty
        self.Ks_err[filename] = self.Ks[filename] * np.sqrt(
            (self.Phi_err[filename]/self.Phi[filename])**2 + (self.D_err[filename]/self.D[filename])**2)

    def generate_plots(self, *args):
        """ Create the permeation tab plots """
        filename = self.current_file.get()
        if self.current_file.get() == "No files yet":
            return

        data = self.datafiles[filename]

        # clear plots
        self.ax.clear()
        self.ax1.clear()
        self.ax12.clear()
        self.ax2.clear()
        self.ax22.clear()
        self.ax3.clear()

        # Permeability/Diffusivity/Solubility vs. Temperature and Flux vs. Pressure (bottom left graph)
        self.PDK_plot(self.ax)

        # Pressure vs. Time (top right graph)
        self.pressure_time_plot(data, filename, self.fig1, self.ax1, self.ax12)

        # Permeation vs. Time (bottom middle graph)
        self.perm_time_plot(data, filename, self.fig2, self.ax2, self.ax22)

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
        self.Phi_label.set(self.Phi[filename])
        self.Prate_label.set(self.Prate[filename])
        self.F_label.set(self.F[filename])
        self.D_label.set(self.D[filename])
        self.D_tlag_label.set(self.D_tlag[filename])
        self.K_label.set(self.Ks[filename])
        self.tss_label.set(self.tss[filename])
        self.Phi_err_label.set(self.Phi_err[filename])
        self.Prate_err_label.set(self.Prate_err[filename])
        self.F_err_label.set(self.F_err[filename])
        self.D_err_label.set(self.D_err[filename])
        self.D_tlag_err_label.set(self.D_tlag_err[filename])
        self.K_err_label.set(self.Ks_err[filename])

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
                fig, axis = plt.subplots(layout="constrained")
                self.PDK_plot(axis)

            elif plot == self.ax1_title:
                # Pressure vs. Time
                fig, ax1 = plt.subplots(layout="constrained")
                ax12 = ax1.twinx()
                self.pressure_time_plot(data, filename, fig, ax1, ax12)

            elif plot == self.ax2_title:
                # Permeation vs. Time
                fig, ax2 = plt.subplots(layout="constrained")
                ax22 = ax2.twinx()
                self.perm_time_plot(data, filename, fig, ax2, ax22)

            elif plot == self.ax3_title:
                # diffusivity comparison plots for optimization
                fig, ax3 = plt.subplots(layout="constrained")
                self.comparison_plot(filename, ax3)

        else:
            showerror(message="Please select a data folder first.")
        plt.show()

    def pressure_time_plot(self, data, filename, fig1, ax1, ax12):
        """ Creates the pressure vs. time (top right) plot """
        # Plot the point where the isolation valve opens and where steady state is determined
        ax1.plot(data.loc[self.t0[filename], 't'], data.loc[self.t0[filename], 'SecP'], 'yo', label="t$_0$")
        ax1.plot(data.loc[self.tss[filename], 't'], data.loc[self.tss[filename], 'SecP'], 'mo', label="t$_{ss}$")

        # Plot secondary pressure and set up axes for it
        ax1.plot(data['t'], data['SecP'], '.', label='Secondary Pressure')
        ax1.set_title(self.ax1_title)
        ax1.set_xlabel(self.ax1_xlabel)
        ax1.set_ylabel(self.ax1_ylabel)

        # Plot primary pressure and set up axes for it
        ax12.plot(data['t'], data['PrimP'], color='orange', label='Primary Pressure')
        ax12.set_ylabel(self.ax12_ylabel)

        # Generate and plot the red leak line (showing which values are used in initial calculations)
        # and green steady state line (showing which values are used in final calculations)
        xleak = data.loc[self.t0[filename] - self.leak_range[filename]:self.t0[filename] - 1, 't']
        yleak = self.pleak_lf[filename]["slope"] * xleak + self.pleak_lf[filename]["intercept"]
        xss = data.loc[self.tss[filename] + 1:self.tss[filename] + self.ss_range[filename], 't']
        yss = self.pss_lf[filename]["slope"] * xss + self.pss_lf[filename]["intercept"]
        ax1.plot(xleak, yleak, color='red', label='Leak')
        ax1.plot(xss, yss, color='lime', label='Steady State')
        # Extrapolate the lines further
        xleak_ex = np.array([data.loc[self.t0[filename] - self.leak_range[filename], 't'],
                             (self.intercept[filename])[0]])
        yleak_ex = self.pleak_lf[filename]["slope"] * xleak_ex + self.pleak_lf[filename]["intercept"]
        xss_ex = np.array([(self.intercept[filename])[0],
                           data.loc[self.tss[filename] + self.ss_range[filename], 't']])
        yss_ex = self.pss_lf[filename]["slope"] * xss_ex + self.pss_lf[filename]["intercept"]
        ax1.plot(xleak_ex, yleak_ex, ":", color='red')
        ax1.plot(xss_ex, yss_ex, ":", color='lime')
        # Plot a point at the intersection of the leak line and steady state line
        ax1.plot(*self.intercept[filename], 'o')

        # Make sure ticks are in the right spot
        ax1.yaxis.set_ticks_position('left')
        ax12.yaxis.set_ticks_position('right')
        ax1.xaxis.set_ticks_position('bottom')

        # Ensure all plots are on the legend. Simply using fig1.legend(...) doesn't get rid of old legends,
        # causing them to stack up in a way I can't figure out how to remove.
        # https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend/47370214#47370214
        pt_lines, pt_labels = ax1.get_legend_handles_labels()
        pt_lines2, pt_labels2 = ax12.get_legend_handles_labels()
        ax12.legend(pt_lines + pt_lines2, pt_labels + pt_labels2, loc="lower right", bbox_to_anchor=(1, 0),
                    bbox_transform=ax1.transAxes, framealpha=0.5)

    def perm_time_plot(self, data, filename, fig2, ax2, ax22):
        """ Creates the permeability vs. time (bottom middle) plot """
        # Plot the permeability over time
        ax2.plot(data.loc[self.t0[filename]:, 't'], self.Perm_Plot[filename].iloc[self.t0[filename]:],
                 label="RT Perm")
        ax2.set_title(self.ax2_title)
        ax2.set_xlabel(self.ax2_xlabel)
        ax2.set_ylabel(self.ax2_ylabel)

        # Plot temperature vs. time and set relevant axes
        ax22.plot(data.loc[self.t0[filename]:, 't'], data.loc[self.t0[filename]:, 'SampT'], color='orange',
                  label='Sample T')
        ax22.set_ylabel(self.ax22_ylabel)

        # Make sure ticks are in the right spot
        ax2.yaxis.set_ticks_position('left')
        ax22.yaxis.set_ticks_position('right')
        ax2.xaxis.set_ticks_position('bottom')

        # Plot a line showing the final calculated permeability
        x = np.linspace(*self.ax2.get_xlim())
        y = self.Phi[filename] * np.ones_like(x)
        ax2.plot(x, y, 'k--', label='SS Perm')
        # Center the permeation plot where it can be easily seen. Added the bottom limit because the permeation was
        # dipping way below 0. Added the top limit was because autoscaling was failing once the bottom limit was added.
        ax2.set_ylim(bottom=0, top=1.25*self.Phi[filename])

        # Ensure all plots are on the legend. Simply using fig2.legend(...) doesn't get rid of old legends,
        # causing them to stack up in a way I can't figure out how to remove.
        # https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend/47370214#47370214
        pt_lines, pt_labels = ax2.get_legend_handles_labels()
        pt_lines2, pt_labels2 = ax22.get_legend_handles_labels()
        ax22.legend(pt_lines + pt_lines2, pt_labels + pt_labels2, loc="lower right", bbox_to_anchor=(1, 0),
                    bbox_transform=ax2.transAxes, framealpha=0.5)

    def comparison_plot(self, filename, ax3):
        """ Creates the diffusivity optimization comparison (bottom right) plot """
        # diffusivity comparison plots for optimization
        ax3.plot(self.D_time[filename][1:], self.lhs[filename][1:], '.', label='Experimental Data', )
        ax3.plot(self.D_time[filename][1:], self.rhs[filename][1:], label='D: Time-Lag')
        ax3.plot(self.D_time[filename][1:], self.rhs_cf[filename][1:], '--', label='D: Optimized')
        ax3.set_title(self.ax3_title)
        ax3.set_xlabel(self.ax3_xlabel)
        ax3.set_ylabel(self.ax3_ylabel)
        ax3.yaxis.set_ticks_position('left')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.legend(framealpha=0.5)

    def PDK_plot(self, ax):
        """ Creates the PDK (bottom left) plot for calculated properties of different files """
        plot = self.current_variable.get()
        artists = []

        # Set plot title
        if plot == "Flux":
            self.ax_title = plot + " vs. Pressure"
            self.ax.loglog()
        else:
            self.ax_title = plot + " vs. Temperature"

        # store as a local variable to allow for a different xlabel in the flux vs. pressure plot
        xlabel = self.ax_xlabel

        # set y labels to get units right
        if plot == "Permeability":
            self.ax_ylabel = plot + " (mol m$^{-1}$ s$^{-1}$ Pa$^{-0.5}$)"
        elif plot == "Diffusivity":
            self.ax_ylabel = plot + " (m$^2$ s$^{-1}$)"
        elif plot == "Solubility":
            self.ax_ylabel = plot + " (mol m$^{-3}$ Pa$^{-0.5}$)"
        elif plot == "Flux":
            self.ax_ylabel = plot + " (mol m$^{-2}$ s$^{-1}$)"
            xlabel = "Pressure (Pa)"

        for filename in self.datafiles.keys():
            x = self.Tsamp0[filename] - self.storage.standard_temp  # temperature
            p = self.PrimP_average[filename]  # pressure
            # y values (and for flux, x values) change based on plot type
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
            elif plot == "Flux":
                y = self.F[filename]
                yerr = self.F_err[filename]
                x = p  # change x-axis to pressure (the name was changed earlier to still work when there is no data)
                label = f'{filename}\n{plot}={y:.3e}' + r' mol m$^{-2}\,$s$^{-1}$' + \
                        f'\nPressure={p:.3e} Pa'
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
        ax.set_xlabel(xlabel)
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
    """ popup window used to change many of the variables used in loading permeation data. """

    def __init__(self, storage, pos, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # import custom widgets
        widgets = Widgets()
        self.add_entry = widgets.add_entry

        self.title("Adjust Persistent Variables")
        # self.resizable(width=False, height=False)
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
        self.pv_filename = os.path.join('data_files', 'persistent_permeation_input_variables.xlsx')

        # Read the Excel sheets into dataframes for easier access
        self.numbers_info = pd.read_excel(self.pv_filename, sheet_name="Numbers", header=0)
        self.misc_info = pd.read_excel(self.pv_filename, sheet_name="MiscInfo", header=0)
        self.GasT_info = pd.read_excel(self.pv_filename, sheet_name="GasT", header=0)
        self.SampT_info = pd.read_excel(self.pv_filename, sheet_name="SampT", header=0)
        self.PrimP_info = pd.read_excel(self.pv_filename, sheet_name="PrimP", header=0)
        self.SecP_info = pd.read_excel(self.pv_filename, sheet_name="SecP", header=0)

        # Pull out variables for display and potential editing
        self.num_GasT = tk.DoubleVar(value=self.numbers_info['GasT'][0])

        # Method of scaling the window width according to number of TCs measuring GasT. At num_GasT = 10,
        # things are still fine. Potentially add a scrollbar instead.
        if self.num_GasT.get() > 2:
            width += int(extra_width * (self.num_GasT.get() - 2))
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

        # Primary pressure
        self.col_PrimP = tk.StringVar(value=self.PrimP_info['PrimP'][0])
        self.m_PrimP = tk.DoubleVar(value=self.PrimP_info['PrimP'][1])
        self.b_PrimP = tk.DoubleVar(value=self.PrimP_info['PrimP'][2])
        self.cerr_PrimP = tk.DoubleVar(value=self.PrimP_info['PrimP'][3])
        self.perr_PrimP = tk.DoubleVar(value=self.PrimP_info['PrimP'][4])

        # Secondary pressure
        self.col_SecP = tk.StringVar(value=self.SecP_info['SecP'][0])
        self.m_SecP = tk.DoubleVar(value=self.SecP_info['SecP'][1])
        self.b_SecP = tk.DoubleVar(value=self.SecP_info['SecP'][2])
        self.cerr_SecP = tk.DoubleVar(value=self.SecP_info['SecP'][3])
        self.perr_SecP = tk.DoubleVar(value=self.SecP_info['SecP'][4])

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
        submit_button.grid(row=14, column=6, columnspan=2, sticky='e')

        quit_button = ttk.Button(self, text="Cancel", command=self.close_pv_win)
        quit_button.grid(row=14, column=8, columnspan=2, sticky='w')

        parent = tk.Frame(self)
        parent.grid(row=14, column=0, columnspan=6, sticky="nsew", pady=5)
        self.add_entry(self, parent, variable=self.inputs, key='num_GasT', text="Number of TCs measuring GasT:",
                       subscript='',
                       tvar=self.num_GasT, units='', row=0, column=0, in_window=True,
                       command=lambda tvar, variable, key, pf: self.one_or_more(tvar, variable, key, pf))

    def one_or_more(self, tvar, variable, key, pf=None):
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
        """ Creates intro text to Permeation Plots tab's Adjust Persistent Variables window """
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
        PPSettingsHelp(pos=(x, y), size=(apv_width, apv_height)).show_help_window()

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
        tk.Label(parent, text="PrimP [Pa]").grid(row=entry_row + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="col_PrimP", text="col", innersubscript="PrimP",
                       innertext=":", subscript="", tvar=self.col_PrimP, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_PrimP", text="m", innersubscript="PrimP",
                       innertext=":", subscript="", tvar=self.m_PrimP, units="[Pa/?]", ent_w=entry_w,
                       row=entry_row + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_PrimP", text="b", innersubscript="PrimP",
                       innertext=":", subscript="", tvar=self.b_PrimP, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="cerr_PrimP", text="cerr", innersubscript="PrimP",
                       innertext=":", subscript="", tvar=self.cerr_PrimP, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 4, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="perr_PrimP", text="perr", innersubscript="PrimP",
                       innertext=":", subscript="", tvar=self.perr_PrimP, units="[%]", ent_w=entry_w,
                       row=entry_row + 5, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        entry_col += 3
        tk.Label(parent, text="SecP [Pa]").grid(row=entry_row + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="col_SecP", text="col", innersubscript="SecP",
                       innertext=":", subscript="", tvar=self.col_SecP, units="", ent_w=entry_w,
                       row=entry_row + 1, column=entry_col,
                       command=lambda tvar, variable, key: self.check_col_names(tvar, variable, key))
        self.add_entry(self, parent, variable=self.inputs, key="m_SecP", text="m", innersubscript="SecP",
                       innertext=":", subscript="", tvar=self.m_SecP, units="[Pa/?]", ent_w=entry_w,
                       row=entry_row + 2, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="b_SecP", text="b", innersubscript="SecP",
                       innertext=":", subscript="", tvar=self.b_SecP, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="cerr_SecP", text="cerr", innersubscript="SecP",
                       innertext=":", subscript="", tvar=self.cerr_SecP, units="[Pa]", ent_w=entry_w,
                       row=entry_row + 4, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="perr_SecP", text="perr", innersubscript="SecP",
                       innertext=":", subscript="", tvar=self.perr_SecP, units="[%]", ent_w=entry_w,
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

        entry_col += 3
        tk.Label(parent, text="Starting Row").grid(row=row2entries + 0, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="Starting Row", text="", innersubscript="",
                       innertext="", subscript="", tvar=self.starting_row, units="", ent_w=entry_w,
                       row=row2entries + 1, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        tk.Label(parent, text="Rows in Footer").grid(row=row2entries + 2, column=entry_col + 1)
        self.add_entry(self, parent, variable=self.inputs, key="Footer Rows", text="", innersubscript="",
                       innertext="", subscript="", tvar=self.footer_rows, units="", ent_w=entry_w,
                       row=row2entries + 3, column=entry_col, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

    def submit(self):
        """Updates the persistent variables with user input"""

        # Create a list that contains all non-column name entries in GasT for checking for NaNs
        GasT_vals = list(self.m_GasT.values()) + list(self.b_GasT.values()) + \
            list(self.cerr_GasT.values()) + list(self.perr_GasT.values())

        # List of all column names for checking for uniqueness
        col_list = [self.col_t.get(), self.col_SampT.get(), self.col_PrimP.get(),
                    self.col_SecP.get(), self.col_IV.get()]
        for val in self.col_GasT.values():
            col_list.append(val.get())

        # check to make sure no numerical variable is assigned "nan"
        if sum(map(np.isnan, [self.m_t.get(), self.b_t.get(), self.cerr_t.get(), self.perr_t.get(),
                              self.m_SampT.get(), self.b_SampT.get(), self.cerr_SampT.get(), self.perr_SampT.get(),
                              self.m_PrimP.get(), self.b_PrimP.get(), self.cerr_PrimP.get(), self.perr_PrimP.get(),
                              self.m_SecP.get(), self.b_SecP.get(), self.cerr_SecP.get(), self.perr_SecP.get(),
                              self.m_IV.get(), self.b_IV.get(), self.starting_row.get(), self.footer_rows.get(),
                              self.num_GasT.get()])) + any(np.isnan(val.get()) for val in GasT_vals) > 0:
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
                self.PrimP_info['PrimP'] = [self.col_PrimP.get(), self.m_PrimP.get(), self.b_PrimP.get(),
                                            self.cerr_PrimP.get(), self.perr_PrimP.get()]
                self.SecP_info['SecP'] = [self.col_SecP.get(), self.m_SecP.get(), self.b_SecP.get(),
                                          self.cerr_SecP.get(), self.perr_SecP.get()]
                self.misc_info['Isolation Valve'] = [self.col_IV.get(), self.m_IV.get(),
                                                     self.b_IV.get(), np.nan, np.nan]
                self.misc_info['Starting Row'] = [self.starting_row.get(), np.nan, np.nan, np.nan, np.nan]
                self.misc_info['Rows in Footer'] = [self.footer_rows.get(), np.nan, np.nan, np.nan, np.nan]

                # Create or remove data from dataframe so that it has the same number of instruments measuring GasT
                # as self.num_GasT. Note that this can handle non-integer values of self.numGasT
                while self.numbers_info['GasT'][0] > self.num_GasT.get():
                    del self.GasT_info['GasT' + str(self.numbers_info['GasT'][0])]
                    self.numbers_info.loc[0, 'GasT'] = self.numbers_info.loc[0, 'GasT'] - 1
                while self.numbers_info['GasT'][0] < self.num_GasT.get():
                    self.GasT_info['GasT' + str(self.numbers_info['GasT'][0] + 1)] = \
                        [self.col_GasT['GasT1'].get(), 1.0, 0.0, 0.0, 0.0]
                    self.numbers_info.loc[0, 'GasT'] = self.numbers_info['GasT'][0] + 1

                # overwrite the old source file with updated dataframes
                with pd.ExcelWriter(self.pv_filename, engine="openpyxl", if_sheet_exists='new', mode='a') as writer:
                    workBook = writer.book

                    # remove's shouldn't be required, but if_sheet_exists='replace' is malfunctioning as of 11/18/21
                    workBook.remove(workBook['Numbers'])
                    workBook.remove(workBook['MiscInfo'])
                    workBook.remove(workBook['GasT'])
                    workBook.remove(workBook['SampT'])
                    workBook.remove(workBook['PrimP'])
                    workBook.remove(workBook['SecP'])

                    # Write updated dataframes to the Excel sheets
                    self.numbers_info.to_excel(writer, sheet_name='Numbers', index=False)
                    self.misc_info.to_excel(writer, sheet_name='MiscInfo', index=False)
                    self.GasT_info.to_excel(writer, sheet_name='GasT', index=False)
                    self.SampT_info.to_excel(writer, sheet_name='SampT', index=False)
                    self.PrimP_info.to_excel(writer, sheet_name='PrimP', index=False)
                    self.SecP_info.to_excel(writer, sheet_name='SecP', index=False)
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
        temp_PrimP_info = self.PrimP_info.copy()
        temp_SecP_info = self.SecP_info.copy()
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
        temp_PrimP_info['PrimP'] = [self.col_PrimP.get(), self.m_PrimP.get(), self.b_PrimP.get(),
                                    self.cerr_PrimP.get(), self.perr_PrimP.get()]
        temp_SecP_info['SecP'] = [self.col_SecP.get(), self.m_SecP.get(), self.b_SecP.get(),
                                  self.cerr_SecP.get(), self.perr_SecP.get()]
        temp_misc_info['Isolation Valve'] = [self.col_IV.get(), self.m_IV.get(), self.b_IV.get(), np.nan, np.nan]
        temp_misc_info['Starting Row'] = [self.starting_row.get(), np.nan, np.nan, np.nan, np.nan]
        temp_misc_info['Rows in Footer'] = [self.footer_rows.get(), np.nan, np.nan, np.nan, np.nan]
        temp_numbers_info['GasT'] = [int(self.num_GasT.get())]  # "int" is so HyPAT recognizes when there was no change

        # List of all column names for checking for uniqueness
        col_list = [self.col_t.get(), self.col_SampT.get(), self.col_PrimP.get(),
                    self.col_SecP.get(), self.col_IV.get()]
        for val in self.col_GasT.values():
            col_list.append(val.get())

        # Check to see if any changes have been made and if there are any duplicate column names
        if temp_misc_info.equals(self.misc_info) and temp_GasT_info.equals(self.GasT_info) and \
                temp_SampT_info.equals(self.SampT_info) and temp_PrimP_info.equals(self.PrimP_info) and \
                temp_SecP_info.equals(self.SecP_info) and temp_numbers_info.equals(self.numbers_info)\
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


class PPSettingsHelp(tk.Toplevel):
    """ popup box to display information about the variables that can be edited. """

    def __init__(self, pos, size, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title("Help for Permeation's Persistent Variables Settings")
        # self.resizable(width=False, height=False)
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
        for item in (" ", "t [s]: ", "PrimP [Pa]: ", "SecP [Pa]: ", "Isolation Valve [1/0]: ", "Starting Row: ",
                     "Rows in Footer: ", "GasT [\u00B0C]: ", "SampT [\u00B0C]: ", "col: ", "m: ", "b: ", "cerr: ",
                     "perr: "):
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
        tdt["PrimP [Pa]: text"].insert("insert", "Primary side pressure.")
        tdt["SecP [Pa]: text"].insert("insert", "Secondary side pressure.")
        tdt["Isolation Valve [1/0]: text"].insert("insert", "The data describing whether the isolation valve is " +
                                                  "closed (~0) or open (not ~0). Start of permeation is " +
                                                  "measured from when the valve first goes from closed to " +
                                                  "open.")
        tdt["Starting Row: text"].insert("insert", "The row at which the program starts reading data from the Excel " +
                                                   "sheet (0-indexed). HyPAT assumes the data file has no header, so " +
                                                   "use this to start analysis after the header.")
        tdt["Rows in Footer: text"].insert("insert", "How many rows at the bottom of the file to ignore.")
        tdt["GasT [\u00B0C]: text"].insert("insert", "The thermocouples (TCs) that measure the gas's temperature. " +
                                                    "The number of TCs can be set arbitrarily via " +
                                                    "'Number of TCs measuring GasT.'")
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

