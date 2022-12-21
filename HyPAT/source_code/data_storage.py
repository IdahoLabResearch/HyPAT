""" Keep track of important variables and calculations for HyPAT, especially the Permeation Estimates tab """
import tkinter as tk
from tkinter import ttk
import numpy as np
import pandas as pd
import os
import platform  # allows for Mac vs. Windows adaptions


class Storage:
    def __init__(self):
        """ Store the variables we need access to anywhere.
            This is initialized in the application class, then passed by reference to the other tabs. """

        self.data = self.load_data()

        # Create dataframes to store values for updating source spreadsheets.
        # First the material_data spreadsheet
        filename = os.path.join('data_files', 'material_data.xlsx')
        self.material_data = pd.read_excel(filename, header=[0, 1],
                                           engine="openpyxl")  # openpyxl supports .xlsx but not .xls
        # set up the index using material names
        self.material_data.set_index(('Unnamed: 0_level_0', "Material"), inplace=True)
        self.material_data.rename_axis("Material", axis="index", inplace=True)
        # Next for the melting_tempK file
        melt_filename = os.path.join('data_files', 'melting_tempK.xlsx')
        self.melting_tempK = pd.read_excel(melt_filename, engine="openpyxl")
        self.melting_tempK.set_index("Material", inplace=True)

        # Read the o-ring data off an Excel sheet into a convenient dataframe.
        # (1/2, 1/4 VCR data comes from https://www.swagelok.com/downloads/webcatalogs/en/ms-01-24.pdf,
        # page 17 Silver Plated Nonretained)
        self.oring_filename = os.path.join('data_files', 'o-ring_data.xlsx')
        self.oring_info = pd.read_excel(self.oring_filename, header=0, index_col=0)
        self.oring_info_4file = self.oring_info.copy()  # Duplicate used exclusively for saving to the o-ring file

        # used parameters
        self.standard_temp = 273.15  # Kelvin
        self.Na = 6.023e23  # Avogadro constant
        self.torr_to_pa = 133.3223684  # Pa/Torr
        self.cc_to_mole = 4.46429e-05  # mol/cc
        self.R = 8.314462618153e-03  # kJ mol^-1 K^-1

        # user parameter inputs
        # sample input
        self.x_samp = tk.DoubleVar(value=1.0)
        self.x_samp2 = tk.DoubleVar()
        self.sample_material = tk.StringVar()

        # o-ring input
        self.oring_type = tk.StringVar()
        self.oring_type.trace_add("write", self.update_oring_info)
        self.outer_diameter = tk.DoubleVar()
        self.inner_diameter = tk.DoubleVar()
        self.oring_thickness = tk.DoubleVar(value=1)
        self.A_perm = tk.DoubleVar()
        self.A_perm2 = tk.DoubleVar()
        self.oring_type.set(list(self.oring_info.index)[0])

        # H mass transport properties
        self.D0 = tk.DoubleVar()
        self.E_D = tk.DoubleVar()
        self.K0 = tk.DoubleVar()
        self.E_K = tk.DoubleVar()
        self.P0 = tk.DoubleVar()
        self.E_P = tk.DoubleVar()
        self.sample_material.trace_add("write", self.update_properties)

        # Read the file for some default values
        self.defaults_filename = os.path.join('data_files', 'default_entry_vals.xlsx')
        self.defaults_info = pd.read_excel(self.defaults_filename, header=0)

        # Calibrated Leak Rate
        self.cal_leak_rate = tk.DoubleVar(value=self.defaults_info["Calibrated Leak Rate [mol s^-1 Torr^-1]"][0])

        # Temperature input
        self.Tc = tk.DoubleVar(value=250)
        self.Tk = tk.DoubleVar()
        self.invTk = tk.DoubleVar()

        # primary input
        self.pP_T2 = tk.DoubleVar(value=750.0)
        self.pP_T2_Pa = tk.DoubleVar()

        # secondary input
        self.sV = tk.DoubleVar(value=self.defaults_info["Secondary Side Volume [cc]"][0])
        self.t_accum = tk.DoubleVar(value=1.00)
        self.t_accum2 = tk.DoubleVar()

        # estimated secondary pressure parameters
        self.t_L = tk.DoubleVar(value=0)  # Estimated time lag
        self.Phi = tk.DoubleVar(value=0)  # Estimated permeability
        self.flux = tk.DoubleVar(value=0)  # Estimated molar flux
        self.flux_atoms = tk.DoubleVar(value=0)  # Estimated  atomic flux
        self.Q = tk.DoubleVar(value=0)  # Estimated molecular permeation rate
        self.del_sP = tk.DoubleVar(value=0)  # Estimated rate of pressure increase, Torr/s
        self.del_sP_Pa = tk.DoubleVar(value=0)  # /\, Pa/s
        self.sP_final = tk.DoubleVar(value=0)  # Estimated final secondary side pressure, Torr
        self.sP_final2 = tk.DoubleVar(value=0)  # /\, Pa
        self.Det_CM = tk.StringVar(value="NOPE")  # Check if detectable via capacitance manometer
        self.p_QMS = tk.DoubleVar(value=0)  # Estimated pressure in QMS
        self.Det_QMS = tk.StringVar(value="NOPE")  # Check if detectable via QMS
        self.Sat_QMS = tk.StringVar(value="NOPE")  # Check if saturate the QMS

        # tolerance used for finding steady state during permeation data analysis and equilibrium during absorption d.a.
        self.tol = tk.DoubleVar(value=1e-6)
        # minimum number of seconds following t0 before finding a steady state or equilibrium
        self.t_del = tk.DoubleVar(value=100)  # Note, as of 5/28/22, this isn't used in absorption_plots.py
        # default number of data points used in determining leak, steady state, and/or equilibrium
        self.gen_dp_range = tk.IntVar(value=30)

        # store the transport properties calculations
        self.TransportParameters = pd.DataFrame()
        self.PTransportParameters = pd.DataFrame()
        self.ATransportParameters = pd.DataFrame()

    def connect_variables(self):
        """ Set up variables so that functions are triggered if they are changed"""
        self.x_samp2.trace_add("write", self.update_pressure_values)
        self.sample_material.trace_add("write", self.update_pressure_values)
        self.A_perm2.trace_add("write", self.update_pressure_values)
        self.invTk.trace_add("write", self.update_pressure_values)
        self.pP_T2_Pa.trace_add("write", self.update_pressure_values)
        self.sV.trace_add("write", self.update_pressure_values)
        self.t_accum2.trace_add("write", self.update_pressure_values)
        self.cal_leak_rate.trace_add("write", self.update_pressure_values)

        # check melting temperature
        self.Tk.trace_add("write", self.check_melt)

    def check_melt(self, *args):
        """ Only some Materials include melting temps.
            We want to verify that the materials used aren't going to melt in practice.
            This currently sets the warning at half the melting point """

        melt = self.data.loc[self.sample_material.get(), ("melting temp. [K]", "")]
        if melt:
            if self.Tk.get() > 1/2 * melt:
                tk.messagebox.showerror("Melting Caution", "The temperature you entered may result in soft metal." +
                                                           " You entered {:.1f} K.".format(self.Tk.get()) +
                                                           " Half the melting temp is {:.1f} K or {:.1f} \u00B0C."
                                                           .format(1/2*melt, 1/2*melt-self.standard_temp))

    @staticmethod
    def load_data():
        """ Read in data from materials files to get the diffusivity and solubility of each material, along with
            the activation energies, minimum temperatures, and maximum temperatures associated with each material's
            diffusivity and solubility. Use that information to calculate those respective number for permeability.
            Read in the melting temp if available. Combine this all into a dataframe and return that dataframe """
        # read in diffusivity and solubility data
        # todo It'd be nice to have the user select a file similar to "SaveFileExample.xlsx" in the data_files folder
        #      then if two of the three quantities (diffusivity, solubility, permeability) are present, the third
        #      will be calculated. This can be similar to how it is done in overview_plots.py - EditMaterials class.
        filename = os.path.join('data_files', 'material_data.xlsx')
        df = pd.read_excel(filename, header=[0, 1], engine="openpyxl")  # openpyxl supports .xlsx file format, not .xls
        # set up the index using material names
        df.set_index(('Unnamed: 0_level_0', "Material"), inplace=True)
        df.rename_axis("Material", axis="index", inplace=True)

        # Clean up numerical errors that creep in at the 14th decimal place
        def round_to_13(df_col):
            list_col = []
            for i in range(len(df_col)):
                list_col.append(round(df_col[i], 13))
            return list_col

        # calculate permeability data
        col = df.columns
        df[("Permeability", "P0 [mol m^-1 s^-1 Pa^-0.5]")] = round_to_13(df[col[0]] * df[col[4]])
        df[("Permeability", "E_P [kJ mol^-1]")] = round_to_13(df[col[1]] + df[col[5]])
        df[("Permeability", "min. temp. [K]")] = [max((df[col[2]][i], df[col[6]][i])) for i in range(len(df))]
        df[("Permeability", "max. temp. [K]")] = [min((df[col[3]][i], df[col[7]][i])) for i in range(len(df))]

        # read in melting temp data
        melt_filename = os.path.join('data_files', 'melting_tempK.xlsx')
        df2 = pd.read_excel(melt_filename, engine="openpyxl")
        df2.set_index("Material", inplace=True)

        # add melting temp data to the other dataframe
        df['melting temp. [K]'] = df2['melting temp. [K]']  # NaN present for Material without the melting points.

        return df

    def update_x_samp(self, tvar, variable, key):
        """ Ensures the new x_samp value is a number, then updates x_samp2"""
        self.check_for_number(tvar, variable, key)
        # updates value in sample input frame
        self.x_samp2.set(self.x_samp.get() * 0.001)
        return True  # need to tell the widget this worked. This is crucial for the validatecommand method.

    def update_oring_info(self, *args):
        """ changes the displayed O-ring values when the user uses the option menu """
        self.outer_diameter.set(self.oring_info.loc[self.oring_type.get(), self.oring_info.columns[0]])
        self.inner_diameter.set(self.oring_info.loc[self.oring_type.get(), self.oring_info.columns[1]])
        self.A_perm.set(np.pi*(self.inner_diameter.get()/2 * 0.001)**2)
        self.A_perm2.set(self.A_perm.get() * 10000)
        self.oring_thickness.set(self.oring_info.loc[self.oring_type.get(), self.oring_info.columns[2]])

    def update_properties(self, *args):
        """ H mass transport properties depend on sample material """
        self.D0.set(self.data.loc[self.sample_material.get(), ("Diffusivity", "D0 [m^2 s^-1]")])
        self.E_D.set(self.data.loc[self.sample_material.get(), ("Diffusivity", "E_D [kJ mol^-1]")])
        self.K0.set(self.data.loc[self.sample_material.get(), ("Solubility", "K0 [mol m^-3 Pa^-0.5]")])
        self.E_K.set(self.data.loc[self.sample_material.get(), ("Solubility", "E_K [kJ mol^-1]")])
        self.P0.set(self.data.loc[self.sample_material.get(), ("Permeability", "P0 [mol m^-1 s^-1 Pa^-0.5]")])
        self.E_P.set(self.data.loc[self.sample_material.get(), ("Permeability", "E_P [kJ mol^-1]")])

        # Update other values using new transport properties
        self.update_pressure_values()

    def update_temperature(self, tvar, variable, key):
        """ compute the temperature in Kelvin and the inverse based on the user-defined temperature in Celsius """
        self.check_for_number(tvar, variable, key)
        self.Tk.set(self.Tc.get() + self.standard_temp)
        self.invTk.set(1/self.Tk.get())
        return True

    def update_pP_T2(self, tvar, variable, key):
        """ Ensures the new pressure value is a number, then updates the Pascal value"""
        self.check_for_number(tvar, variable, key)
        self.pP_T2_Pa.set(self.torr_to_pa * self.pP_T2.get())
        return True

    def update_accumulation_time(self, tvar, variable, key):
        """ Ensures the new time value is a number, then updates the hour value"""
        self.check_for_number(tvar, variable, key)
        self.t_accum2.set(self.t_accum.get() * 3600)
        return True

    def check_for_number(self, tvar, variable, key, scientific_format=False, parent_frame=None):
        """ Checks if the user entry is a number (float or int). If not, revert the entered string to its prior state"""
        try:
            svar = float(variable[key].get())  # string var for formatting
        except ValueError:
            # i.e. user typed a word or left it blank
            tk.messagebox.showwarning("Invalid Entry", "Please enter a number.")
            # If this entry is in a window (i.e. the variable in_window has been passed to the add_entry function as
            # True), then make that window pop back up after the warning is done.
            if parent_frame:
                parent_frame.deiconify()
            svar = tvar.get()
        # Either set the entry to what it was prior to user entry (if user entry is not a number) or to the user entry
        if key in variable.keys():
            variable[key].delete(0, "end")
            if not scientific_format:
                variable[key].insert(0, "{}".format(svar))
            else:
                variable[key].insert(0, "{:.2e}".format(svar))
            tvar.set(svar)
        return True

    def update_pressure_values(self, *args):
        """ On any user input, update the values of variables in
            the estimated secondary pressure section.
            Variable names and equations copied from Excel. This method reduced extra calls to the tk variables. """
        self.t_L.set(self.x_samp2.get()**2 / 6 / (
                self.D0.get() * np.exp(-1 * self.E_D.get() * self.invTk.get() / self.R)))  # s, Time-lag
        self.Phi.set(self.P0.get() *
                     np.exp(-1 * self.E_P.get() * self.invTk.get() / self.R))  # mol(Q2)/m/s/sqrt(Pa), Molecular permeability
        self.flux.set(self.Phi.get() *
                      np.sqrt(self.pP_T2_Pa.get()) / self.x_samp2.get())  # mol(Q2)/s/m^2, Molecular permeation flux
        self.flux_atoms.set(self.flux.get() * self.Na * 2)  # atoms/s/m^2, Atomic permeation flux
        self.Q.set(self.flux.get() * self.A_perm.get())  # mol(Q2)/s, Molecular permeation rate
        self.del_sP.set((self.Q.get() * self.R * (self.standard_temp + 20) / self.sV.get()) *
                        760)  # Torr/s, rate of pressure increase, assumes room temperature gas
        self.del_sP_Pa.set(self.del_sP.get() * self.torr_to_pa)  # Pa/s, /\
        self.sP_final.set(self.del_sP.get() * self.t_accum2.get())  # Torr, final secondary pressure
        self.sP_final2.set(self.del_sP_Pa.get() * self.t_accum2.get())  # Pa, /\
        # Detectable with capacitance manometer (>1E-4 Torr)
        if self.sP_final.get() > 0.0001:
            self.Det_CM.set("YES")
        else:
            self.Det_CM.set("NO")

        self.p_QMS.set(self.Q.get() / self.cal_leak_rate.get())
        # Detectable with QMS (>1E-10 Torr)
        if self.p_QMS.get() > 0.0000000001:
            self.Det_QMS.set("YES")
        else:
            self.Det_QMS.set("NO")
        # Saturate with QMS  (>1E-6 Torr)
        if self.p_QMS.get() > 0.000001:
            self.Sat_QMS.set("SATURATE")
        else:
            self.Sat_QMS.set("OK")


class Widgets:

    """ Container for custom widgets: entry, label. """

    def __init__(self):

        # design (these font variables are also used in other tabs)
        self.font = ("TkDefaultFont", 15, "bold")
        self.font2 = ("TkDefaultFont", 18, "bold")
        self.bold = ("TkDefaultFont", 12, "bold")

    def view_value(self, event, tvar, variable, key):
        """ Allows user to see full value of a formatted value when entry is clicked """
        if float(variable[key].get()) != tvar.get():
            variable[key].delete(0, "end")
            variable[key].insert(0, "{}".format(tvar.get()))

    def add_entry(self, parent_frame, frame, variable, key, text, subscript, tvar, units, command, column=0, row=0,
                  ent_w=None, formatting=False, subsubscript='', innertext='', innersubscript='', in_window=False):
        """ Creates an entry widget along with accompanying text. Configures that widget.
            key, vals, and units can all be lists, if we want the value expressed in more than one unit """
        # Add opening text
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript, row, column)

        # Create user entry
        variable[key + "var"] = tk.StringVar(value=tvar.get())  # Variable that stores user entry
        variable[key] = tk.Entry(frame, textvariable=variable[key + "var"], justify="center")
        # Establish the validate command depending on whether formatting and/or in_window are/is true
        if not in_window and not formatting:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar, y=variable, z=key: command(x, y, z))
            command(tvar, variable, key)
        elif in_window and not formatting:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar, y=variable, z=key, pf=parent_frame:
                                 command(x, y, z, pf=pf))
            command(tvar, variable, key, pf=parent_frame)
        elif in_window and formatting:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar, y=variable, z=key, format=formatting, pf=parent_frame:
                                 command(x, y, z, format, pf=pf))
            command(tvar, variable, key, formatting, pf=parent_frame)
            variable[key].bind("<FocusIn>", lambda ev, x=tvar, y=variable, z=key: self.view_value(ev, x, y, z))
        else:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar, y=variable, z=key, format=formatting: command(x, y, z,
                                                                                                              format))
            command(tvar, variable, key, formatting)
            variable[key].bind("<FocusIn>", lambda ev, x=tvar, y=variable, z=key: self.view_value(ev, x, y, z))

        # Finish configuring entry
        variable[key].bind("<Return>", lambda x: parent_frame.focus_set())
        variable[key].bind("<KP_Enter>", lambda x: parent_frame.focus_set())
        variable[key].config(bg="yellow", width=ent_w)
        variable[key].grid(row=row, column=column + 1)
        tk.Label(frame, text=units).grid(row=row, column=column + 2, sticky='w')

    def add_entry2(self, parent_frame, frame, variable, key, text, subscript, tvar1, tvar2, units, conversion,
                   row=0, column=0, ent_w=None, formatting=False, subsubscript='', innertext='', innersubscript=''):
        """ adds a row with a box for user entry, and displays that value in
            two sets of units, i.e. [mm] and [m]

            variable: dictionary that stores entry widgets and the strings of the entries
            key: specific name for widget
            text: entry label for the user
            subscript: subscript of the text (set to "" if unwanted)
            tvar1/tvar2: floats that are used for storing the numerical values of the entries and
                         during calculations elsewhere
            units: list of units, in order of appearance
            conversion: function that defines the unit transformation
            - the function passed into conversion must return true
            Optional:
            row, col: initial row and column of entry sequence. Column will increase by 1 for each part of sequence
            ent_w: width of the entry widget
            formatting: Boolean that says whether this variable should be formatted in scientific notation
            subsubscript: subscript to the subscript
            innersubscript, innertext: label of the entry is formatted as follows:
                                        text + inner subscript + inner text + subscript + subscript to the subscript
        """
        # Add opening text
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript, row, column)

        # Create user entry
        variable[key + "var"] = tk.StringVar(value=tvar1.get())  # Variable that stores user entry
        variable[key] = tk.Entry(frame, textvariable=variable[key + "var"], justify="center")
        # Establish the validate command depending on whether formatting is True
        if not formatting:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar1, y=variable, z=key: conversion(x, y, z))
            conversion(tvar1, variable, key)
        else:
            variable[key].config(validate="focusout",
                                 validatecommand=lambda x=tvar1, y=variable, z=key, format=formatting:
                                 conversion(x, y, z, format))
            conversion(tvar1, variable, key, formatting)
            # Enable viewing full variable if entry is focused on
            variable[key].bind("<FocusIn>", lambda ev, x=tvar1, y=variable, z=key: self.view_value(ev, x, y, z))
        conversion(tvar1, variable, key)

        # Redirect focus to general window if either enter key (keyboard or keypad) is pressed
        variable[key].bind("<Return>", lambda x: parent_frame.focus_set())
        variable[key].bind("<KP_Enter>", lambda x: parent_frame.focus_set())
        # for active Entry to be highlighted
        variable[key].config(bg="yellow", width=ent_w)
        variable[key].grid(row=row, column=column + 1)
        # first units
        tk.Label(frame, text=units[0]).grid(row=row, column=column + 2)
        # second units
        FormatLabel(frame, textvariable=tvar2, borderwidth=1, relief="ridge", format="{:.2e}"
                    ).grid(row=row, column=column + 3)
        tk.Label(frame, text=units[1]).grid(row=row, column=column + 4, sticky='w')

    def add_entry3(self, parent_frame, frame, variable, key, text, subscript, tvar1, tvar2, units, update_function,
                   row=0, column=0, ent_w=10, subsubscript='', innertext='', innersubscript=''):
        """ add a row with two entry boxes, one for the value and the other for the error.
         expects the two units to be the same, so only provide one.

          Naming scheme:
            variable[key] is the entry widget for the data
            variable[key+"_err"] is the second entry widget for the error
            variable[key + "var"] is the stringVar for formatting the data
            variable[key + "err"] is the stringVar for formatting the error
            tvar 1 and 2 are DoubleVar for storing the data and error

        the set of string variables is necessary to make the formatting work.
        -note, this was really finicky, so be cautious making changes.
          """
        # need a temporary storage for the tvar's. tvar1 and tvar2 will store the data in DoubleVar, while
        # another variable will store a StringVar for display

        # store scientific notation value for variable and its error
        variable[key + "var"] = tk.StringVar(value=tvar1.get())
        variable[key + "err"] = tk.StringVar(value=tvar2.get())

        # text preceding the entry
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript, row, column)
        # entry for the value here
        variable[key] = tk.Entry(frame, textvariable=variable[key + "var"], validate="focusout", justify="center",
                                 validatecommand=lambda x=tvar1, y='var', z=key: update_function(x, y, z), width=ent_w)

        # Enable viewing full variable if entry is focused on
        variable[key].bind("<FocusIn>", lambda ev, x=tvar1, y=variable, z=key: self.view_value(ev, x, y, z))
        # Redirect focus to general window if either enter key (keyboard or keypad) is pressed
        variable[key].bind("<Return>", lambda x: parent_frame.focus_set())
        variable[key].bind("<KP_Enter>", lambda x: parent_frame.focus_set())
        # active Entry will be highlighted
        variable[key].config(bg="yellow")
        variable[key].grid(row=row, column=column+1, sticky="ew")
        # separator to indicate the error is the next box
        tk.Label(frame, text=" +/- ").grid(row=row, column=column+2, sticky="ew")
        # second entry for the error
        key2 = key + "_err"
        variable[key2] = tk.Entry(frame, textvariable=variable[key + "err"], width=ent_w,
                                  validate="focusout", justify="center",
                                  validatecommand=lambda x=tvar2, y='err', z=key: update_function(x, y, z))
        variable[key2].bind("<Return>", lambda x: parent_frame.focus_set())
        variable[key2].bind("<KP_Enter>", lambda x: parent_frame.focus_set())

        # Enable viewing full variable if entry is focused on
        variable[key2].bind("<FocusIn>", lambda ev, x=tvar2, y=variable, z=key2: self.view_value(ev, x, y, z))
        # for active Entry to be highlighted
        variable[key2].config(bg="yellow")
        variable[key2].grid(row=row, column=column+3, sticky="ew")
        # first units
        tk.Label(frame, text=units).grid(row=row, column=column+4, sticky="w")

        # display in scientific notation
        variable[key].delete(0, "end")
        variable[key].insert(0, "{:.2e}".format(tvar1.get()))
        variable[key + "_err"].delete(0, "end")
        variable[key + "_err"].insert(0, "{:.2e}".format(tvar2.get()))

    @staticmethod
    def add_text0(frame, text, subscript, row=0, column=0, sticky='e',
                  subsubscript='', innertext='', innersubscript=''):
        """ Create text without creating additional labels and/or entry boxes """
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript,
                                       row, column, sticky)

    def add_text(self, frame, text, subscript, tvar, units, row=0, column=0,
                 subsubscript='', innertext='', innersubscript=''):
        """ similar to add_entry, but for displaying important variables the user can't change """
        # text preceding the entry
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript, row, column)

        # Add the little box with the number. If the number in that box is changeable elsewhere in the program via
        # being a DoubleVar variable, then use the class FormatLabel to deal with it.
        if isinstance(tvar, tk.DoubleVar):
            FormatLabel(frame, textvariable=tvar, borderwidth=1, relief="ridge", width=10, format="{:.2e}"
                        ).grid(row=row, column=column + 1, sticky="ew")
        else:
            tk.Label(frame, textvariable=tvar, borderwidth=1, relief="ridge", width=10, font=self.bold
                     ).grid(row=row, column=column + 1, sticky="ew")
        # Add the units
        tk.Label(frame, text=units).grid(row=row, column=column + 2, columnspan=2, sticky='w')

    @staticmethod
    def add_text2(frame, text, subscript, tvar1, tvar2, units, row=0, column=0,
                  subsubscript='', innertext='', innersubscript=''):
        """ Create a pair of labels with text before, after, and in between them.
            Expects tuple or string for units. If string is passed the unit is duplicated so
            both labels will display the same units and a +/- text widget will be placed between the labels"""
        if type(units) == str:
            # this function is to show a value and its error
            units = (" +/- ", units)
            mid_sticky = "ew"
        else:
            mid_sticky = "w"

        # text preceding the label
        Widgets.add_text_and_subscript(frame, text, subscript, subsubscript, innertext, innersubscript, row, column)

        # use the class FormatLabel to deal with the DoubleVar variables while creating these labels.
        FormatLabel(frame, textvariable=tvar1, borderwidth=1, relief="ridge", width=10, format="{:.2e}"
                    ).grid(row=row, column=column + 1, sticky="ew")
        tk.Label(frame, text=units[0]).grid(row=row, column=column + 2, sticky=mid_sticky)
        FormatLabel(frame, textvariable=tvar2, borderwidth=1, relief="ridge", width=10, format="{:.2e}"
                    ).grid(row=row, column=column + 3, sticky="ew")
        tk.Label(frame, text=units[1]).grid(row=row, column=column + 4, sticky='w')

    def add_text_and_subscript(self, text, subscript, subsubscript='', innertext='', innersubscript='',
                               row=0, column=0, sticky='e'):
        from tkinter import font
        t = tk.Text(self, height=1, borderwidth=0, background=self.cget("background"), spacing3=1)

        # todo This method creates a lot of fonts, but allows the subscripts to be in the same font as the parents.
        #      Find out if the many fonts is a problem.
        # Either set the text box's font to the parent frame's font or set it to TkDefaultFont
        try:
            # frame.cget("font") gets the name of the font the frame is using, then font.nametofont() takes that name
            # and turns it into an item that python recognizes.
            tfont = font.nametofont(self.cget("font")).copy()
        except tk.TclError:
            # Expected error if the parent frame doesn't have a font, such as the pop-up windows
            tfont = font.nametofont("TkDefaultFont").copy()
        t.config(font=tfont)

        # Create fonts for subscripts
        stfont, sstfont = tfont.copy(), tfont.copy()
        stfont.config(size=tfont.cget("size") - 2)  # subscript font
        sstfont.config(size=tfont.cget("size") - 4)  # sub-subscript font

        # Set width of text box while calculating required size and adapting according to os
        if platform.system() == 'Darwin':
            divisor = 8
            addend = 0
        else:
            divisor = 6
            addend = 2
        # tfont.measure() says how many pixels text+subscript would take to write out.
        # For Windows, //6 takes that number and turns it into how many zeros it would take to be that long because one
        #   zero takes up 6 pixels in this font. The +2 adds one more character to make up for // rounding down and
        #   another character to make room for a space.
        # For Mac, setting divisor to 8 and addend to 0 seems to look good, but it's unclear why those numbers work well
        # The unit of width is characters, specifically 0s, not pixels, which is why text + subscript + ... is needed
        t.configure(width=tfont.measure(text + subscript + subsubscript +
                                        innertext + innersubscript) // divisor + addend)

        # Set up tags for the text box
        t.tag_configure("right", justify="right")  # Right aligns the text
        # Create the subscripts. Offset determines its height, font gives it the correct font and a smaller font size
        t.tag_configure("subscript", offset=-2, justify="right", font=stfont)
        t.tag_configure("sub-subscript", offset=-4, justify="right", font=sstfont)
        t.insert("insert", text, "right", innersubscript, "subscript", innertext, "right",
                 subscript, "subscript", subsubscript + " ", "sub-subscript")

        t.grid(row=row, column=column, sticky=sticky)

        t.configure(state="disabled")  # Turn off editing the text


class FormatLabel(tk.Label):
    """
    https://stackoverflow.com/questions/59408126/tkinter-formatting-doublevars-for-display
    """

    def __init__(self, master=None, cnf={}, **kw):

        # default values
        self._format = '{}'
        self._textvariable = None

    # get new format and remove it from `kw` so later `super().__init__` doesn't use them (it would get error message)
        if 'format' in kw:
            self._format = kw['format']
            del kw['format']

        # get `textvariable` to assign own function which set formatted text in Label when variable change value
        if 'textvariable' in kw:
            self._textvariable = kw['textvariable']
            self._textvariable.trace_add("write", self._update_text)
            del kw['textvariable']

        # run `Label.__init__` without `format` and `textvariable`
        super().__init__(master, cnf={}, **kw)

        # update text after running `Label.__init__`
        if self._textvariable:
            # self._update_text(None, None, None)
            self._update_text(self._textvariable, '', 'w')

    def _update_text(self, a, b, c):
        """update text in label when variable change value"""
        self["text"] = self._format.format(self._textvariable.get())


class LoadingScreen(tk.Toplevel):
    """  This class creates and adjusts a progress bar.
         Majority of the text for this class obtained from
         https://github.com/IdahoLabResearch/CDB_Analysis_Program/blob/main/CDB-AP_1.0/src/PlotModule.py """
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.title('Loading')
        self.update()

    def add_progress_bar(self, num_of_files):
        """ Create the progress bar and accompanying label"""
        self.progress_bar = ttk.Progressbar(self, orient=tk.HORIZONTAL,
                                            length=250, mode='determinate')
        self.progress_bar.pack(pady=20)
        self.file_num = num_of_files
        self.label = tk.Label(self, text="0/{} files".format(self.file_num))
        self.label.pack()
        self.update()

    def update_progress_bar(self, amount, i):
        """ Update the progress bar and accompanying label"""
        self.progress_bar['value'] += amount
        self.label['text'] = "{}/{} files".format(i, self.file_num)
        self.update()