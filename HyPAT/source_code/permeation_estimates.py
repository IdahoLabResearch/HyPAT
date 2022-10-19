""" Main page to accept input from the user. This is the Permeation Estimates tab in HyPAT """
import tkinter as tk
from tkinter import ttk, font
from .data_storage import Widgets
# Imports for the ORingsAndDefaultVals class
import numpy as np
import os
import platform


class InputForm(tk.Frame):

    def __init__(self, parent, storage, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # container for important variables
        self.storage = storage

        # store widgets
        self.widgets = Widgets()
        self.add_entry = self.widgets.add_entry
        self.add_entry2 = self.widgets.add_entry2
        self.add_text0 = self.widgets.add_text0
        self.add_text = self.widgets.add_text
        self.add_text2 = self.widgets.add_text2

        self.input = {}  # stores all entry and text variables

        # design
        self.font = self.widgets.font
        self.named_font = font.Font(font=self.font, name="named_font")
        self.font2 = self.widgets.font2

        # generate the sub frames
        self.create_parameter_input_frame(row=0, column=0)
        self.create_secondary_est_frame(row=0, column=1)
        self.create_initial_conditions_frame(row=1, column=0)
        self.create_final_output_frame(row=1, column=1)
        # configure to allow resizing
        self.rowconfigure((0, 1), weight=1)
        self.columnconfigure((0, 1), weight=1)

        # perform the initial calculation of Secondary Pressure estimates
        self.storage.update_pressure_values()
        self.storage.connect_variables()

        # update drop down menus upon clicking on the tab
        self.bind("<FocusIn>", self.update_option_menu)

    def create_parameter_input_frame(self, row, column):
        """
        top left frame containing Sample Input, H Mass Transport Properties, O-Ring Input, and Calibrated Leak Rate
        :param row:
        :param column:
        :return:
        """
        parent = tk.LabelFrame(self, text="PARAMETER INPUTS", font=self.font2, bd=10)
        tk.Label(parent, text="PLEASE INPUT ALL IN YELLOW.", background="yellow").grid(row=0, column=0)
        self.create_sample_input_frame(parent, row=1, column=0)
        self.create_properties_frame(parent, row=2, column=0)
        self.create_o_ring_input_frame(parent, row=3, column=0)
        self.create_cal_leak_rate_frame(parent, row=4, column=0)
        parent.grid(row=row, column=column, columnspan=1, sticky="nsew")

        for n in range(4):
            parent.rowconfigure(n, weight=1)
        parent.columnconfigure(0, weight=1)

    def create_sample_input_frame(self, parent, row, column):
        """ Sets up frame for sample specific information """
        box1 = tk.LabelFrame(parent, text="SAMPLE INPUT", fg="blue")

        # Sample thickness user entry
        self.add_entry2(self, box1, variable=self.input, key="x_samp", text="Sample thickness: x", subscript="samp",
                        tvar1=self.storage.x_samp, tvar2=self.storage.x_samp2, units=("[mm]", "[m]"),
                        conversion=lambda tvar1, variable, key: self.storage.update_x_samp(tvar1, variable, key), row=1)

        # Drop down menu for sample material
        self.add_text0(box1, text="Sample material:", subscript="", row=3)
        self.storage.sample_material.set(self.storage.data.index[2])
        self.input["sample_material"] = tk.OptionMenu(box1, self.storage.sample_material,
                                                      *list(self.storage.data.index))
        self.input["sample_material"].config(bg="yellow", indicatoron=False)
        self.input["sample_material"].grid(row=3, column=1, sticky="ew")

        box1.grid(row=row, column=column, sticky="nsew")

    def update_option_menu(self, *args):
        """ refresh the list of materials """
        # delete current options
        menu = self.input["sample_material"]["menu"]
        menu.delete(0, tk.END)
        # add new options
        for material in list(self.storage.data.index):
            menu.add_command(label=material, command=lambda value=material: self.storage.sample_material.set(value))

    def create_o_ring_input_frame(self, parent, row, column):
        """ Sets up O-ring input frame """
        # O-ring input
        box3 = tk.LabelFrame(parent, text="O-RING INPUT", fg="blue")
        # drop down menu for O-ring type
        self.add_text0(box3, text="O-ring type:", subscript="", row=0)
        self.input["o-ring type"] = tk.OptionMenu(box3, self.storage.oring_type,
                                                  *list(self.storage.oring_info.index))
        self.input["o-ring type"].config(bg="yellow", indicatoron=False)
        self.input["o-ring type"].grid(row=0, column=1, columnspan=4, sticky="ew")

        # O-ring characteristics
        self.add_text(box3, text="O-ring outer diameter: D", subscript="samp",
                      tvar=self.storage.outer_diameter, units="[mm]", row=1)
        self.add_text(box3, text="O-ring inner diameter: d", subscript="samp",
                      tvar=self.storage.inner_diameter, units="[mm]", row=2)
        self.add_text(box3, text="O-ring thickness: x", subscript="samp",
                      tvar=self.storage.oring_thickness, units="[mm]", row=3)
        self.add_text2(box3, text="Permeation surface area: A", subscript="perm",
                       tvar1=self.storage.A_perm, tvar2=self.storage.A_perm2,
                       units=("[m\u00b2]", "[cm\u00b2]"), row=4)  # "[m^2]", "[cm^2]"
        box3.grid(row=row, column=column, rowspan=1, sticky="nsew")

    def create_properties_frame(self, parent, row, column):
        """ Sets up H mass transport properties frame """
        box2 = tk.LabelFrame(parent, text="H MASS TRANSPORT PROPERTIES", fg="blue")
        self.add_text(box2, text="Pre-exponential factor for diffusivity: D", subscript="0",
                      units="[m\u00b2 s\u207b\u00b9]", tvar=self.storage.D0, row=0)  # "[m^2 s^-1]"
        self.add_text(box2, text="Activation energy for diffusivity: E", subscript="D",
                      units="[kJ mol\u207b\u00b9]", tvar=self.storage.E_D, row=1)  # "[kJ mol^-1]",
        self.add_text(box2, text="Pre-exponential factor for solubility: K", subscript="0",
                      units="[mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]",  # "[mol m^-3 Pa^-0.5]"
                      tvar=self.storage.K0, row=2)
        self.add_text(box2, text="Activation energy for solubility: E", subscript="K",
                      units="[kJ mol\u207b\u00b9]", tvar=self.storage.E_K, row=3)  # "[kJ mol^-1]"
        self.add_text(box2, text="Pre-exponential factor for permeability: \u03A6", subscript="0",
                      # units="[mol m^-1 s^-1 Pa^-0.5]"
                      units="[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]",
                      tvar=self.storage.P0, row=4)
        self.add_text(box2, text="Activation energy for permeability: E", subscript="\u03A6",
                      units="[kJ mol\u207b\u00b9]", tvar=self.storage.E_P, row=5)  # "[kJ mol^-1]",

        box2.grid(row=row, column=column, rowspan=1, sticky="nsew")

    def create_cal_leak_rate_frame(self, parent, row, column):
        """ Sets up the calibrated leak rate frame """
        box4 = tk.LabelFrame(parent, text="CALIBRATED LEAK RATE", fg="blue")
        self.add_entry(self, box4, self.input, key="cal_leak_r",
                       text="Calibrated leak rate:", subscript='',
                       units="[mol s\u207b\u00b9 Torr\u207b\u00b9]",
                       tvar=self.storage.cal_leak_rate, row=0, formatting=True,
                       command=lambda tvar, variable, key, formatting: self.storage.check_for_number(tvar, variable,
                                                                                                     key, formatting))
        box4.grid(row=row, column=column, sticky="nsew")

    def create_secondary_est_frame(self, row, column):
        """ Sets up frame for the right column of frames. These are related to Secondary Permeation Flux Estimation """
        parent = tk.LabelFrame(self, text="SECONDARY PERMEATION FLUX ESTIMATION", font=self.font2, bd=10)
        tk.Label(parent, text="PLEASE INPUT ALL IN YELLOW.", bg="yellow").grid(row=0)
        self.create_temperature_input_frame(parent, row=1, column=0)
        self.create_primary_input_frame(parent, row=2, column=0)
        self.create_secondary_input_frame(parent, row=3, column=0)
        self.create_pressure_frame(parent, row=4, column=0)

        parent.grid(row=row, column=column, columnspan=1, sticky="nsew")

        # Button for editing default values and O-ring list
        settings_b = ttk.Button(parent, text='Settings', command=self.edit_orings_and_default_vals)
        settings_b.grid(row=0, sticky="e")

        for n in range(4):
            parent.rowconfigure(n, weight=1)  # Make every row expand/contract based on available space
        parent.columnconfigure(0, weight=1)

    def edit_orings_and_default_vals(self):
        """ Calls the window that allows for editing the list of O-rings and the default values of some variables.
            Updates things as necessary if changes were made. """
        self.update()
        # get location of gui
        x = self.winfo_rootx()
        y = self.winfo_rooty()
        # Call the O-Rings and Default Values window and, if something was changed, refresh appropriate values
        o_change, d_change = ORingsAndDefaultVals(storage=self.storage, pos=(x, y)).show_pv_window()
        if o_change:
            # Reset the stored values
            self.storage.oring_type.set(list(self.storage.oring_info.index)[0])
            self.storage.update_oring_info()

            # Reset the drop-down menu
            self.input["o-ring type"]['menu'].delete(0, "end")
            for string in self.storage.oring_info.index:
                self.input["o-ring type"]['menu'].\
                    add_command(label=string, command=lambda value=string: self.storage.oring_type.set(value))

        if d_change:
            # As in the populate function, both the entry and the variable seem to need to be updated separately
            self.input['cal_leak_r'].delete(0, "end")
            self.input['cal_leak_r'].insert(0, "{}".format(
                self.storage.defaults_info["Calibrated Leak Rate [mol s^-1 Torr^-1]"][0]))
            self.storage.cal_leak_rate.set(self.storage.defaults_info["Calibrated Leak Rate [mol s^-1 Torr^-1]"][0])

            self.input['sV'].delete(0, "end")
            self.input['sV'].insert(0, "{}".format(
                self.storage.defaults_info["Secondary Side Volume [cc]"][0]))
            self.storage.sV.set(self.storage.defaults_info["Secondary Side Volume [cc]"][0])

    def create_temperature_input_frame(self, parent, row, column):
        """ Sets up the temperature input frame"""
        box1 = tk.LabelFrame(parent, text="TEMPERATURE INPUT", fg="blue")
        self.add_entry2(self, box1, self.input, key="temp", text="Temperature: T", innersubscript="C",
                        innertext=", T", subscript="K",
                        units=("[\u00B0C]", "[K]"), tvar1=self.storage.Tc, tvar2=self.storage.Tk, row=0,
                        conversion=lambda tvar1, variable, key: self.storage.update_temperature(tvar1, variable, key))
        self.add_text(box1, text="Reciprocal temperature: 1/T", subscript="K",
                      units="[1/K]", tvar=self.storage.invTk, row=1)
        box1.grid(row=row, column=column, sticky="nsew")

    def create_primary_input_frame(self, parent, row, column):
        """ Sets up frame for primary pressure input """
        box2 = tk.LabelFrame(parent, text="PRIMARY INPUT", fg="blue")

        self.add_entry2(self, box2, variable=self.input, key="pP_T2",
                        text="Primary calculated partial pressure: pP", subscript="T", subsubscript="2",
                        units=("[Torr]", "[Pa]"), tvar1=self.storage.pP_T2, tvar2=self.storage.pP_T2_Pa, row=0,
                        conversion=lambda tvar1, variable, key: self.storage.update_pP_T2(tvar1, variable, key))
        box2.grid(row=row, column=column, sticky="nsew")

    def create_secondary_input_frame(self, parent, row, column):
        """ Sets up frame for secondary input """
        box3 = tk.LabelFrame(parent, text="SECONDARY INPUT", fg="blue")
        self.add_entry(self, box3, self.input, key="sV", text="Secondary side volume: V", subscript="sec",
                       units="[cc]", tvar=self.storage.sV, row=0,
                       command=lambda tvar, variable, key: self.storage.check_for_number(tvar, variable, key))
        self.add_entry2(self, box3, self.input, key="t_accum", text="Accumulation time: t", subscript="accum",
                        units=("[hr]", "[sec]"), tvar1=self.storage.t_accum,
                        tvar2=self.storage.t_accum2, row=1,
                        conversion=lambda tvar1, variable, key:
                        self.storage.update_accumulation_time(tvar1, variable, key))
        box3.grid(row=row, column=column, sticky="nsew")

    def create_pressure_frame(self, parent, row, column):
        """ Sets up frame for displaying predicted values """
        box4 = tk.LabelFrame(parent, text="ESTIMATED SECONDARY PRESSURE BY REFERENCE PERMEABILITY DATA", fg="blue")
        self.add_text(box4, text="Estimated time-lag: t", subscript="L",
                      units="[sec]", tvar=self.storage.t_L, row=0)
        self.add_text(box4, text="Estimated permeability: \u03A6", subscript="",
                      # units="[mol m^-1 s^-1 Pa^-0.5]"
                      units="[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]",
                      tvar=self.storage.Phi, row=1)
        self.add_text2(box4, text="Estimated permeation flux: J", subscript="inf", subsubscript="",
                       # units=("[mol m^-2 s^-1]", "[atoms m^-2 s^-1]")
                       units=("[mol m\u207b\u00b2 s\u207b\u00b9]", "[atoms m\u207b\u00b2 s\u207b\u00b9]"),
                       tvar1=self.storage.flux, tvar2=self.storage.flux_atoms, row=2)
        self.add_text(box4, text="Estimated permeation rate: Q", subscript="", subsubscript="",
                      units="[mol s\u207b\u00b9]",  # "[mol s^-1]",
                      tvar=self.storage.Q, row=4)
        self.add_text2(box4, text="Estimated rate of pressure increase: dP/dt", subscript="",
                       units=("[Torr s\u207b\u00b9]", "[Pa s\u207b\u00b9]"),  # ("[Torr s^-1]", "[Pa s^-1]")
                       tvar1=self.storage.del_sP, tvar2=self.storage.del_sP_Pa, row=5)
        self.add_text2(box4, text="Estimated final secondary pressure: sP", subscript="final",
                       units=("[Torr]", "[Pa]"), tvar1=self.storage.sP_final,
                       tvar2=self.storage.sP_final2, row=6)
        self.add_text(box4, text="Detectable with capacitance manometer (>1E-4 Torr):", subscript="",
                      units="YES=detectable, NO=undetectable", tvar=self.storage.Det_CM, row=7)
        self.add_text(box4, text="Estimated pressure in QMS: P", subscript="", subsubscript="QMS",
                      units="[Torr]", tvar=self.storage.p_QMS, row=8)
        self.add_text(box4, text="Detectable with QMS (>1E-10 Torr):", subscript="",
                      units="YES=detectable, NO=undetectable", tvar=self.storage.Det_QMS, row=9)
        self.add_text(box4, text="Saturate with QMS (>1E-6 Torr):", subscript="",
                      units="SATURATE or OK", tvar=self.storage.Sat_QMS, row=10)
        box4.grid(row=row, column=column, sticky="nsew", rowspan=3)

    def create_initial_conditions_frame(self, row, column):
        """ Sets up frame for displaying initial conditions """
        parent = tk.LabelFrame(self, text="INITIAL CONDITIONS", fg="blue", bd=10)
        self.add_text0(parent, text="Reference database column:", subscript="", row=0)
        tk.Label(parent, textvariable=self.storage.sample_material, borderwidth=1,
                 relief="ridge", width=10).grid(row=0, column=1)

        self.add_text2(parent, text="Sample thickness: x", subscript="samp",
                       units=("[mm]", "[m]"), tvar1=self.storage.x_samp,
                       tvar2=self.storage.x_samp2, row=1)
        self.add_text(parent, text="Temperature: T", subscript="K",
                      units="[K]", tvar=self.storage.Tk, row=2)
        self.add_text2(parent, text="Accumulation time: t", subscript="accum",
                       units=("[hr]", "[sec]"), tvar1=self.storage.t_accum, tvar2=self.storage.t_accum2, row=3)
        parent.config(font=self.named_font)
        parent.grid(row=row, column=column, sticky="nsew")

    def create_final_output_frame(self, row, column):
        """ Sets up frame for displaying final output """
        parent = tk.LabelFrame(self, text="FINAL OUTPUT", fg="blue", bd=10)
        self.add_text(parent, text="Estimated time-lag: t", subscript="L",
                      units="[sec]", tvar=self.storage.t_L, row=0)
        self.add_text(parent, text="Estimated permeation rate: Q", subscript="", subsubscript="",
                      units="[mol s\u207b\u00b9]",  # "[mol s^-1]",
                      tvar=self.storage.Q, row=1)
        self.add_text(parent, text="Estimated pressure in QMS: P", subscript="", subsubscript="QMS",
                      units="[Torr]", tvar=self.storage.p_QMS, row=2)
        self.add_text(parent, text="Estimated final secondary pressure: sP", subscript="final",
                      units="[Torr]", tvar=self.storage.sP_final, row=3)

        parent.config(font=self.named_font)
        parent.grid(row=row, column=column, sticky="nsew")


class ORingsAndDefaultVals(tk.Toplevel):
    """ popup box used to add, edit, and remove O-rings in the source spreadsheet, along with editing the default
        values of some variables """

    def __init__(self, storage, pos, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # import custom widgets
        widgets = Widgets()
        self.add_entry = widgets.add_entry

        self.title("Add or Edit O-Rings / Edit Default Values")
        # self.resizable(width=False, height=False)
        self.minsize(400, 170)

        # gui_x/y values determined by running self.updateidletasks() at the end of self.__init__ and then printing size
        gui_x = 1429
        gui_y = 742

        # Set window size and position
        if platform.system() == 'Darwin':
            width = 1100
            height = 212
        else:
            width = 875
            height = 170
        pos_right = int(pos[0] + (gui_x - width) / 2)
        pos_down = int(pos[1] + (gui_y - height) / 3)
        self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        self.storage = storage
        self.orings_changed = False  # True if anything related to O-rings was changed
        self.defaults_changed = False  # True if anything related to the default values of other variables was changed

        # store information
        self.ring = tk.StringVar()

        self.D_ring = tk.DoubleVar(value=np.NaN)  # Outer diameter of O-ring
        self.d_ring = tk.DoubleVar(value=np.NaN)  # Inner diameter of O-ring
        self.x_ring = tk.DoubleVar(value=np.NaN)  # Thickness of O-ring

        self.defaults_info = self.storage.defaults_info
        # Calibrated leak rate
        self.cal_leak_r = tk.DoubleVar(value=self.defaults_info["Calibrated Leak Rate [mol s^-1 Torr^-1]"][0])
        # secondary volume
        self.sV = tk.DoubleVar(value=self.defaults_info["Secondary Side Volume [cc]"][0])

        # store entry values and Combobox widget
        self.inputs = {}

        # populate the window
        self.create_labels()
        self.create_inputs()

        quit_button = ttk.Button(self, text="Cancel", command=self.destroy)
        quit_button.grid(row=0, column=6, sticky='e')

    def create_labels(self):
        """ Creates intro text to Add or Edit O-Rings Edit Default Values window """
        parent = tk.Frame(self)
        parent.grid(row=0, column=0, columnspan=3, sticky="nsew")

        tk.Label(parent, text="Please enter the O-ring's information or edit the default values (but not both at once)."
                 ).grid(sticky="nsew", columnspan=2)

    def create_inputs(self):
        """ Creates a frame for O-rings and a frame for Default Values and populates those frames """

        # Create O-ring frame

        parent = tk.Frame(self, borderwidth=10, relief="groove")
        parent.grid(row=1, column=0, columnspan=3, sticky="nsew")

        tk.Label(parent, text="O-Ring:").grid(row=0, column=0, sticky="e")

        # allow user to add new O-ring, edit an existing O-ring, or remove an existing O-ring
        self.inputs['ring'] = ttk.Combobox(parent, textvariable=self.ring, width=17)
        self.inputs['ring']['values'] = [""] + list(self.storage.oring_info.index)
        self.inputs['ring'].bind('<<ComboboxSelected>>', self.populate_labels)
        self.inputs['ring'].grid(row=0, column=1)

        self.add_entry(self, parent, variable=self.inputs, key="D_ring", text="Outer Diameter:",
                       subscript="", tvar=self.D_ring, units="[mm]", row=1, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="d_ring", text="Inner Diameter:",
                       subscript="", tvar=self.d_ring, units="[mm]", row=2, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.add_entry(self, parent, variable=self.inputs, key="x_ring", text="Thickness:",
                       subscript="", tvar=self.x_ring, units="[mm]", row=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        submit_button = ttk.Button(parent, text="Submit O-Ring Edits/Additions", command=self.submit_4orings)
        submit_button.grid(row=4, column=0, sticky="e")

        remove_button = ttk.Button(parent, text="Remove O-Ring", command=self.remove_oring)
        remove_button.grid(row=4, column=1)

        # Create default values frame

        parent2 = tk.Frame(self, borderwidth=10, relief='groove')
        parent2.grid(row=1, column=4, columnspan=3, sticky="nsew")
        tk.Label(parent2, text="Default Values").grid(row=0, column=4)
        self.add_entry(self, parent2, variable=self.inputs, key="cal_leak_rate", text="Calibrated Leak Rate:",
                       subscript="", tvar=self.cal_leak_r, units="[mol s^-1 Torr^-1]",
                       row=1, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        tk.Label(parent2, text="").grid(row=2, column=0)

        self.add_entry(self, parent2, variable=self.inputs, key="sV", text="Secondary Side Volume:",
                       subscript="", tvar=self.sV, units="[cc]", row=3, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))

        submit_defaults_button = ttk.Button(parent2, text="Submit Default Values", command=self.submit_4default)
        submit_defaults_button.grid(row=4, column=5, sticky="w")

    def populate_labels(self, *args):
        """ load the data based on user choice of O-ring """
        ring = self.ring.get()
        if ring != "":  # If an actual O-ring was selected...
            # disable editing in combobox
            self.inputs['ring']['state'] = 'readonly'

            # fill in entry boxes
            ''' Deleting the current contents of of the entry box and updating it also updates 
                self.inputs['D_ring' + 'var'], the variable that contains the entry box's string. However, it does this
                without calling the verification function, which is why self.D_ring needs to get updated separately.'''
            self.inputs['D_ring'].delete(0, "end")
            self.inputs['D_ring'].insert(0, "{}".format(
                self.storage.oring_info.loc[ring, self.storage.oring_info.columns[0]]))
            self.D_ring.set(self.storage.oring_info.loc[ring, self.storage.oring_info.columns[0]])

            self.inputs['d_ring'].delete(0, "end")
            self.inputs['d_ring'].insert(0, "{}".format(
                self.storage.oring_info.loc[ring, self.storage.oring_info.columns[1]]))
            self.d_ring.set(self.storage.oring_info.loc[ring, self.storage.oring_info.columns[1]])

            self.inputs['x_ring'].delete(0, "end")
            self.inputs['x_ring'].insert(0, "{}".format(
                self.storage.oring_info.loc[ring, self.storage.oring_info.columns[2]]))
            self.x_ring.set(self.storage.oring_info.loc[ring, self.storage.oring_info.columns[2]])

        else:  # If the user selects the blank box...
            # enable editing
            self.inputs['ring']['state'] = 'normal'
            # clear entries
            ''' Deleting the current contents of of the entry box and updating it also updates 
                self.inputs['D_ring' + 'var'], the variable that contains the entry box's string. However, it does this
                without calling the verification function, which is why self.D_ring needs to get updated separately.'''
            self.inputs['D_ring'].delete(0, "end")
            self.inputs['D_ring'].insert(0, "{}".format(np.NaN))
            self.D_ring.set(np.NaN)

            self.inputs['d_ring'].delete(0, "end")
            self.inputs['d_ring'].insert(0, "{}".format(np.NaN))
            self.d_ring.set(np.NaN)

            self.inputs['x_ring'].delete(0, "end")
            self.inputs['x_ring'].insert(0, "{}".format(np.NaN))
            self.x_ring.set(np.NaN)

    def submit_4orings(self):
        """ Handles potential new O-rings and the editing of old O-rings.
            Note: submitting via this function means that changes to Default Values aren't saved"""
        # Check if all the entries have been filled in. If not, don't submit.
        if sum(map(np.isnan, [self.D_ring.get(), self.d_ring.get(), self.x_ring.get()])) > 0:
            tk.messagebox.showwarning(title="Missing Entry",
                                      message="Please fill in remaining boxes before continuing.")
            self.deiconify()
        else:  # if all entries have been filled in...
            # popup message dependent on if it's an edit or an addition
            if self.ring.get() in self.storage.oring_info.index:
                c_message = "Edit the recorded properties of {} in the source file for the O-rings?".format(
                    self.ring.get()) + \
                            '\n \u2022 Selecting "yes" will change them here and in the original Excel sheet.' + \
                            '\n \u2022 Selecting "no" will change them here but not in the original spreadsheet.' + \
                            '\n \u2022 Selecting "cancel" will return you to the Add or Edit O-Rings window.'
            else:
                c_message = "Add {} to the source file for the O-rings?".format(self.ring.get()) + \
                            '\n \u2022 Selecting "yes" will add it here and in the original Excel sheet.' + \
                            '\n \u2022 Selecting "no" will add it here but not in the original spreadsheet.' + \
                            '\n \u2022 Selecting "cancel" will return you to the Add or Edit O-Rings window.'
            ans = tk.messagebox.askyesnocancel(
                title="Confirmation",
                message=c_message)
            if ans:
                # Change the selected O-ring to have the new properties (or create that O-ring if it didn't exist)
                self.storage.oring_info.loc[self.ring.get()] = [self.D_ring.get(), self.d_ring.get(), self.x_ring.get()]

                # Change the O-ring dataframe reserved for editing the source file
                self.storage.oring_info_4file.loc[self.ring.get()] = \
                    [self.D_ring.get(), self.d_ring.get(), self.x_ring.get()]

                # overwrite the old source file with the updated dataframe
                oring_filename = os.path.join('datafiles', 'o-ring_data.xlsx')
                self.storage.oring_info_4file.to_excel(oring_filename)

                self.orings_changed = True
                self.destroy()
            elif ans is False:
                # Change the selected O-ring to have the new properties (or create that O-ring if it didn't exist)
                self.storage.oring_info.loc[self.ring.get()] = [self.D_ring.get(), self.d_ring.get(), self.x_ring.get()]

                self.orings_changed = True
                self.destroy()
            else:
                self.deiconify()  # bring back the pop-up window

    def remove_oring(self):
        """ Handles removal of old O-rings.
            Note: removing via this function means that changes to Default Values aren't saved"""
        ans = tk.messagebox.askyesnocancel(
            title="Confirmation",
            message='Do you want to also remove {} from the source file for O-rings?'.format(self.ring.get()) +
                    '\n \u2022 Selecting "yes" will remove it from here and the original Excel sheet.' +
                    '\n \u2022 Selecting "no" will remove it from here but not the original spreadsheet.' +
                    '\n \u2022 Selecting "cancel" will return you to the Add or Edit O-Rings window.')

        if ans:
            try:
                # Remove the O-ring from the dataframe that influences the program (storage.oring_info) and from the
                # dataframe that will be used to replace the source file (which only gets edited when if specifically
                # selected to be edited)
                self.storage.oring_info.drop(index=self.ring.get(), inplace=True)
                if self.ring.get() in self.storage.oring_info_4file.index:  # Check that the O-ring's in this dataframe
                    self.storage.oring_info_4file.drop(index=self.ring.get(), inplace=True)
                    # overwrite the old source file with the updated dataframe
                    oring_filename = os.path.join('datafiles', 'o-ring_data.xlsx')
                    self.storage.oring_info_4file.to_excel(oring_filename)

                self.orings_changed = True
                self.destroy()  # Closes O-rings/Default Values window
            except KeyError:
                # If no O-ring has been selected
                tk.messagebox.showwarning(title="Missing Selection",
                                          message="Please select an O-ring to remove.")
                self.deiconify()
        elif ans is False:
            try:
                # Remove the O-ring from the dataframe that influences the program, storage.oring_info
                self.storage.oring_info.drop(index=self.ring.get(), inplace=True)

                self.orings_changed = True
                self.destroy()  # Closes O-rings/Default Values window
            except KeyError:
                # If no O-ring has been selected
                tk.messagebox.showwarning(title="Missing Selection",
                                          message="Please select an O-ring to remove.")
                self.deiconify()
        else:
            self.deiconify()  # bring back the pop-up window

    def submit_4default(self):
        """ Handles changes to default values.
            Note: submitting via this function means that changes to O-rings aren't saved"""
        #  Ensure all default value entries have been filled
        if sum(map(np.isnan, [self.cal_leak_r.get(), self.sV.get()])) > 0:
            tk.messagebox.showwarning(title="Missing Entry",
                                      message="Please fill in remaining boxes before continuing.")
            self.deiconify()
        else:
            ans = tk.messagebox.askokcancel(
                title="Confirmation",
                message="Do you wish the save the default values as they are to the source file?")
            if ans:
                # Set the dataframe with the new values
                self.storage.defaults_info.loc[0, "Calibrated Leak Rate [mol s^-1 Torr^-1]"] = self.cal_leak_r.get()
                self.storage.defaults_info.loc[0, "Secondary Side Volume [cc]"] = self.sV.get()

                # overwrite the old source file with the updated dataframe
                default_entry_vars_filename = os.path.join('datafiles', 'default_entry_vals.xlsx')
                self.storage.defaults_info.to_excel(default_entry_vars_filename, index=False)

                self.defaults_changed = True
                self.destroy()
            else:
                self.deiconify()  # bring back the pop-up window

    def show_pv_window(self):
        """ This allows the window to open and then, when closed, trigger the main gui to update the O-rings or
            default values (if something was changed)"""
        self.deiconify()
        self.wait_window()

        return self.orings_changed, self.defaults_changed
