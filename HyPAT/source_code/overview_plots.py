""" Code for the Overview Plots tab in HyPAT """
import matplotlib
matplotlib.use("TkAgg")  # backend of matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
import tkinter as tk
from tkinter import ttk, colorchooser
import mplcursors
import platform
if platform.system() == 'Darwin':
    from tkmacosx import Button as MacButton  # As of August 2021, tk.Button didn't work right for Macs
    BUTTON_WIDTH = 50
    SCROLL_SCALE = int(1)
else:
    MacButton = tk.Button
    BUTTON_WIDTH = 5  # account for the different button
    SCROLL_SCALE = 120  # Speed of scrolling in Materials frame
from tkinter.filedialog import asksaveasfilename, askopenfilename
import os
from .data_storage import Widgets
from pandas import DataFrame, read_excel


class Plots(tk.Frame):

    def __init__(self, parent, storage, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # container for important variables
        self.storage = storage
        self.index = list(self.storage.data.index)
        self.columns = list(self.storage.data.columns)
        self.item = {}  # Container for storing whether a material is selected for plotting
        self.exp_data = DataFrame()  # Container for data received from permeation_plots
        self.start = -1

        self.plotted_once = False
        # Keep track of changes in materials for use in setting the Materials' frame scroll region
        self.added_material = False
        self.subtracted_material = False

        # Containers
        self.check_buttons = []
        self.colors = {}  # for choosing plot line color
        self.color_buttons = {}
        # list of matplotlib colors
        self.base_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.current_base_color_idx = 0  # for default plotting
        self.artists1, self.artists2 = [], []
        self.experimental_data0 = {}  # Container for the diffusivity from each file loaded into permeation_plots
        self.experimental_data1 = {}  # Container for the solubility from each file loaded into permeation_plots
        self.experimental_data2 = {}  # Container for the permeability from each file loaded into permeation_plots

        # make the frames
        self.create_materials_frame()
        self.label = tk.Label(self)
        # self.label.grid(row=1, column=0)
        self.create_plot_frame(row=0, column=1)
        self.plot()

    def create_materials_frame(self):
        """ Create the materials frame during initial set up """

        self.mats_frame = tk.LabelFrame(self, text="Materials", font=("TkDefaultFont", 18, "bold"), bd=10,
                                        labelanchor="n", bg='grey93')

        # add a canvas for the scroll bar
        self.mat_canvas = tk.Canvas(self.mats_frame, bg='grey93')
        self.mat_canvas.grid(row=0, column=0, sticky="nsew")

        # add a scroll bar to the canvas
        scroll_bar = tk.Scrollbar(self.mats_frame, orient="vertical", command=self.mat_canvas.yview)
        scroll_bar.grid(row=0, column=2, sticky='ns')
        self.mat_canvas.configure(yscrollcommand=scroll_bar.set)

        # create another frame to contain the buttons
        self.button_frame = tk.Frame(self.mat_canvas, bg='grey93')
        self.mat_canvas.create_window((0, 0), window=self.button_frame, anchor='nw')
        self.button_frame.bind("<Configure>", self.reset_scroll_region)  # allow new materials to appear

        # bind scroll wheel, but only when mouse is hovering in the button frame
        self.button_frame.bind('<Enter>', self._bound_to_mousewheel)
        self.button_frame.bind('<Leave>', self._unbound_to_mousewheel)

        # Add materials to the group of materials that can be selected for plotting and set up their respective buttons
        keys = self.item.keys()
        for n, material in enumerate(self.index):
            # default has just the materials with melting temperatures. In default file, just the pure metals
            self.add_material(material, n, keys)

        self.mats_frame.rowconfigure(0, weight=1)
        self.mats_frame.columnconfigure(0, weight=1)

        # update buttons to calc size
        self.button_frame.update_idletasks()
        # set scrolling region
        self.reset_scroll_region()

        f2 = tk.Frame(self.mats_frame)
        f2.grid(row=1)

        # add button to add or edit (or remove) materials
        new_mat_button = ttk.Button(f2, text="Add/Edit Material", command=self.add_or_edit_materials)
        new_mat_button.grid(row=0)

        # add button to save new spreadsheet
        excel_update_button = ttk.Button(f2, text="Save File", command=self.save_materials_in_excel)
        excel_update_button.grid(row=0, column=1)

        self.mats_frame.grid(row=0, column=0, sticky="nsew")

        self.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

    def plot_permeation(self):
        """ plot/remove the experimental permeation data in the overview plots """

        # Store the old dataframe for when we need to remove old data after new data is selected in permeation_plots
        old_data = self.exp_data
        # Get the new dataframe from the latest selected folder
        self.exp_data = self.storage.PermeationParameters
        if self.exp_data.empty and old_data.empty:
            # If no data has been loaded into HyPAT
            return

        if self.b.config('text')[-1] == 'Add Experimental Data':
            if self.exp_data.empty:
                return

            # Plot data from permeation_plots tab
            for file in self.exp_data.index:
                self.experimental_data2[file] = self.ax[2].plot(
                    1000/self.exp_data.loc[file, "Sample Temperature [K]"],
                    self.exp_data.loc[file, "Permeability [mol m^-1 s^-1 Pa^-0.5]"], '.')
            for file in self.exp_data.index:
                self.experimental_data0[file] = self.ax[0].plot(
                    1000/self.exp_data.loc[file, "Sample Temperature [K]"],
                    self.exp_data.loc[file, "Diffusivity [m^2 s^-1]"], '.')
            for file in self.exp_data.index:
                self.experimental_data1[file] = self.ax[1].plot(
                    1000/self.exp_data.loc[file, "Sample Temperature [K]"],
                    self.exp_data.loc[file, "Solubility [mol m^-3 Pa^-0.5]"], '.')
            self.b.config(text='Remove Experimental Data')  # Toggle the text so next time, it goes to the "else" part
        else:
            if old_data.empty:
                return
            # Remove all the plotted data points from the plots
            for file in old_data.index:
                self.experimental_data0[file].pop(0).remove()
                self.experimental_data1[file].pop(0).remove()
                self.experimental_data2[file].pop(0).remove()
            self.b.config(text='Add Experimental Data')  # Toggle the text so next time, it goes to the "if" part
        self.canvas.draw()
        self.toolbar.update()

    def toggle_multiple_labels(self):
        """ Toggles between being able to see labels while hovering over lines vs being able to see many labels
            by clicking on the lines """
        if self.b2.config('text')[-1] == 'Enable Multiple Labels':
            mplcursors.Cursor.remove(self.cursor)  # Removes the possibility of double labels
            self.cursor = mplcursors.cursor(self.artists2, multiple=True)  # Enable multiple labels
            # self.cursor.connect("add", lambda sel: sel.annotation.draggable(True))  # Doesn't work for some reason
            self.cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))  # Simplify text
            self.b2.config(text='Disable Multiple Labels')  # Toggle the text so next time, it goes to the "else" part
        else:
            mplcursors.Cursor.remove(self.cursor)  # Removes the possibility of double labels
            self.cursor = mplcursors.cursor(self.artists2, hover=2, multiple=False, highlight=True)  # Enable hovering
            self.cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))  # Simplify text
            self.b2.config(text='Enable Multiple Labels')  # Toggle the text so next time, it goes to the "if" part
        self.canvas.draw()
        self.toolbar.update()

    def fit_to_Arrhenius_eqn(self, new_data=False):
        """ Fit the diffusivity, solubility, and permeability derived at a variety of sample temperatures to an
            Arrhenius fit to determine the pre-exponential factors and activation energies """
        from scipy.optimize import curve_fit

        # define data for analysis according to whether it comes from the permeation plots tab or a user-selected file
        if self.exp_data.empty or self.b.config('text')[-1] == 'Add Experimental Data' or new_data:
            # If there is nothing in the dataframe that holds data or if the button for displaying experimental data is
            # toggled so that no data is currently displayed or if the button was clicked to select a new file...
            fit_filename = askopenfilename(filetypes=[('Excel File', '*.xlsx')])
            if not fit_filename:  # if no file is selected...
                return
            data = read_excel(fit_filename, header=0, engine="openpyxl")
        else:
            data = self.exp_data

        # Read in variables
        try:
            T = data["Sample Temperature [K]"]
            P_var = data["Permeability [mol m^-1 s^-1 Pa^-0.5]"]
            P_unc = data["Permeability Uncertainty [mol m^-1 s^-1 Pa^-0.5]"]
            D_var = data["Diffusivity [m^2 s^-1]"]
            D_unc = data["Diffusivity Uncertainty [m^2 s^-1]"]
            S_var = data["Solubility [mol m^-3 Pa^-0.5]"]
            S_unc = data["Solubility Uncertainty [mol m^-3 Pa^-0.5]"]
            R = 0.0083145  # kJ / (mol K)
        except KeyError:  # This is expected if the file doesn't have the correct headers
            tk.messagebox.showerror(title="Arrhenius Fit Error",
                                    message='File must be formatted like the files generated'
                                            ' by the "Export to Excel" button in the "Permeation Plots" tab."')
            return

        # The Arrhenius function
        def f(temp, energy, variable):
            return variable * np.exp(-energy / (R * temp))

        # Get weighted fits
        try:
            P_fit, P_error = curve_fit(f, T, P_var, p0=[50, P_var[0]], sigma=P_unc)
            P_uncert = np.sqrt(np.diag(P_error))
            D_fit, D_error = curve_fit(f, T, D_var, p0=[50, D_var[0]], sigma=D_unc)
            D_uncert = np.sqrt(np.diag(D_error))
            S_fit, S_error = curve_fit(f, T, S_var, p0=[50, S_var[0]], sigma=S_unc)
            S_uncert = np.sqrt(np.diag(S_error))
        except ValueError as e:  # Mostly to handle blank cells, but also handles other errors
            tk.messagebox.showerror(title="Arrhenius Fit Error", message="Value Error: " + str(e))
            return
        except Exception as e:  # in case of any other exception
            # If you want to see the traceback as part of the error text, change the below value to true
            want_traceback = False
            if want_traceback:
                import traceback
                full_traceback = "\n" + traceback.format_exc()
            else:
                full_traceback = ""
            tk.messagebox.showerror(title="Arrhenius Fit Error", message='Unknown Error with file. The following' +
                                          ' exception was raised: "' + str(e) + '".' + full_traceback + '\n\n')
            return

        # Create popup window to display results
        popup = tk.Tk()
        popup.wm_title("Hydrogen Transport Properties Arrhenius Fit")
        # Text for textbox
        fit_texts = "Weighted Arrhenius fit of selected data.\n\n" + \
                    "Pre-exponential factor for diffusivity: {:.2e} +/- {:.2e}".format(D_fit[1], D_uncert[1]) + \
                    " [m\u00b2 s\u207b\u00b9]\n" + \
                    "Activation energy for diffusivity: {:.2e} +/- {:.2e}".format(D_fit[0], D_uncert[0]) + \
                    " [kJ mol\u207b\u00b9]\n\n" + \
                    "Pre-exponential factor of solubility: {:.2e} +/- {:.2e}".format(S_fit[1], S_uncert[1]) + \
                    " [mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]\n" + \
                    "Activation energy for solubility: {:.2e} +/- {:.2e}".format(S_fit[0], S_uncert[0]) + \
                    " [kJ mol\u207b\u00b9]\n\n" + \
                    "Pre-exponential factor for permeability: {:.2e} +/- {:.2e}".format(P_fit[1], P_uncert[1]) + \
                    " [mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]\n" + \
                    "Activation energy for permeability: {:.2e} +/- {:.2e}".format(P_fit[0], P_uncert[0]) + \
                    " [kJ mol\u207b\u00b9]"

        # Create a textbox with the text in it
        from tkinter import font
        if platform.system() == 'Darwin':
            twidth = 60
        else:
            twidth = 75
        fittext = tk.Text(popup, font=font.nametofont("TkDefaultFont"), background=self.mats_frame.cget("background"),
                          padx=10, pady=10, spacing2=10, width=twidth, height=11)
        fittext.insert("insert", fit_texts)
        fittext.configure(state="disabled")  # Turn off editing the text
        fittext.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        # Add the buttons
        Arrhenius_button_frame = tk.Frame(popup)
        Arrhenius_button_frame.grid(row=1, column=0)
        b1 = tk.Button(master=Arrhenius_button_frame, text="Analyze Arrhenius for New Data",
                       command=lambda loading=True: self.fit_to_Arrhenius_eqn(loading))
        b1.grid(row=0, column=0)
        b2 = tk.Button(master=Arrhenius_button_frame, text="Save to Excel",
                       command=lambda pf=P_fit, pu=P_uncert, df=D_fit, du=D_uncert, sf=S_fit, su=S_uncert:
                       self.save_Arrhenius_to_excel(pf, pu, df, du, sf, su))
        b2.grid(row=0, column=1)

    def save_Arrhenius_to_excel(self, P_fit, P_uncert, D_fit, D_uncert, S_fit, S_uncert):
        """ Save the pre-exponential factors and activation energies determined via Arrhenius fit to Excel """
        filename = asksaveasfilename(initialdir=os.path.dirname(__file__), defaultextension=".xlsx")
        if filename != '':
            data = {'Value': [D_fit[1], D_fit[0], S_fit[1], S_fit[0], P_fit[1], P_fit[0]],
                    'Uncertainty': [D_uncert[1], D_uncert[0], S_uncert[1], S_uncert[0], P_uncert[1], P_uncert[0]],
                    'Units': ['[m\u00b2 s\u207b\u00b9]', '[kJ mol\u207b\u00b9]',
                              '[mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]', '[kJ mol\u207b\u00b9]',
                              '[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]', '[kJ mol\u207b\u00b9]']}

            # Creates pandas DataFrame.
            df = DataFrame(data, index=['Pre-exponential factor for diffusivity:', 'Activation energy for diffusivity:',
                                        'Pre-exponential factor for solubility:', 'Activation energy for solubility:',
                                        'Pre-exponential factor for permeability:', 'Activation energy for permeability:'])
            df.to_excel(filename)

    def reset_scroll_region(self, event=None):
        """ Derives initial scroll region and adapts the scroll region to the addition or removal of materials"""
        if self.added_material or not self.plotted_once or self.subtracted_material:
            bbox = self.mat_canvas.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.
            ROWS, COLS = len(list(self.storage.data.index)), 2  # Size of grid.
            ROWS_DISP = 25  # Number of rows to display.
            COLS_DISP = 2  # Number of columns to display.
            # Define the scrollable region as entire canvas with only the desired
            # number of rows and columns displayed.
            w, h = bbox[2] - bbox[1], bbox[3] - bbox[1]
            dw, dh = int((w / COLS) * COLS_DISP), int((h / ROWS) * ROWS_DISP)
            self.mat_canvas.configure(scrollregion=bbox, width=dw, height=dh)

            # don't reshape again until we have another change in materials
            self.added_material = False
            self.subtracted_material = False

    def add_material(self, material, n, keys):
        """ Handle the initial set up and further edits of the Materials in the overview plots tab
            material: string variable containing the material's name that is being added
            n: integer containing the material's row number
            keys: index of the list of materials; list of strings containing all the names of materials added so far """
        # Add the material to the list of materials if it isn't on it
        if material not in keys:
            # Determine whether the material is selected by default by checking if the melting point is included
            if not np.isnan(self.storage.data.loc[material, self.columns[-1]]):
                self.item[material] = tk.IntVar(value=1)  # Selected by default
            else:
                # NaN present for material without the melting points.
                self.item[material] = tk.IntVar(value=0)  # Not selected by default
        # Set up the check button of the material
        self.check_buttons.append(tk.Checkbutton(self.button_frame, text=material, variable=self.item[material],
                                                 bg='grey93', command=self.plot))
        self.check_buttons[n].bind("<Button-1>", self.selectstart)
        self.check_buttons[n].bind("<Shift-Button-1>", self.selectrange)
        self.check_buttons[n].grid(row=n, column=1, sticky="nsw")

        # Set up the color button of the material
        if material not in self.colors.keys():
            # assign colors the first time we create a plot
            color = self.cycle_base_colors()
            self.colors[material] = tk.StringVar(value=color)
        else:
            color = self.colors[material].get()  # In case the user changed the color at some point

        self.color_buttons[material] = MacButton(self.button_frame,
                                                 command=lambda data=material: self.change_color(data),
                                                 width=BUTTON_WIDTH, bg=color, fg=color,
                                                 highlightbackground=color,
                                                 activebackground=color)
        self.color_buttons[material].grid(row=n, column=0, sticky="nsew")
        self.button_frame.rowconfigure(n, weight=1)

    def _on_mousewheel(self, event):
        """ set scroll speed of mouse wheel """
        # need to divide event.delta by 120 on Windows - controls scroll speed
        # have to use integer division, as regular / returns a float
        self.mat_canvas.yview_scroll(-1 * event.delta//SCROLL_SCALE, "units")

    def _bound_to_mousewheel(self, event):
        """ enable scrolling via mouse wheel """
        self.mat_canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _unbound_to_mousewheel(self, event):
        """ disable scrolling via mouse wheel """
        self.mat_canvas.unbind_all("<MouseWheel>")

    def add_or_edit_materials(self):
        """ Calls the window that allows for editing the list of materials.
            Updates things as necessary if changes were made. """
        self.update()
        # get location of gui
        x = self.winfo_rootx()
        y = self.winfo_rooty()
        # When the Add or Edit Materials window is closed, "changed" is True if something was changed, and
        # edited_material is a string variable with the edited material's name in it
        changed, edited_material = EditMaterials(storage=self.storage, pos=(x, y)).show_aem_window()
        if changed:
            # check if a material has been added (or, for the elif, removed)
            if len(self.storage.data.index) == len(self.index) + 1:
                # update list of materials
                self.index = list(self.storage.data.index)

                # Update other relevant lists and dictionaries
                row = len(self.index) - 1  # use zero based index for row number
                keys = self.item.keys()
                self.add_material(edited_material, row, keys)

                self.added_material = True
            elif len(self.storage.data.index) == len(self.index) - 1:
                # update various dictionaries and lists
                element_pos = self.index.index(edited_material)  # Get the position of the removed material
                self.index = list(self.storage.data.index)
                self.check_buttons[element_pos].destroy()  # Destroys the buttons
                # Note that both destroying the button and removing its spot from the list/dict seem necessary
                self.check_buttons.pop(element_pos)
                self.color_buttons[edited_material].destroy()
                self.color_buttons.pop(edited_material)
                self.colors.pop(edited_material)
                self.item.pop(edited_material)

                # Update the buttons' positions in the grid to account for the removed material
                for n, material in enumerate(self.index):
                    self.check_buttons[n].grid(row=n)
                    self.color_buttons[material].grid(row=n)

                self.subtracted_material = True
            self.plot()

    def save_materials_in_excel(self):
        """ Save the dataframe containing the diffusivity, solubility, and permeability of each material to
            an Excel file """
        filename = asksaveasfilename(initialdir=os.path.dirname(__file__), defaultextension=".xlsx")
        if filename != '':
            self.storage.data.to_excel(filename)

    def create_plot_frame(self, row=0, column=0):
        """ entirely separate frame containing the diffusivity, solubility, and permeability plots """
        # make the plot in its own frame to avoid mixing pack and grid
        plot_frame = tk.Frame(self)
        self.fig = Figure()
        self.ax = [self.fig.add_subplot(1, 3, i) for i in range(1, 4)]
        # The below line helps keep the graphs looking nice when the window gets made small. Note that tight_layout is
        # an experimental feature as of 10/16/2021 and may be changed unpredictably.
        self.fig.set_tight_layout("True")

        # Format axes
        self.axC = [self.ax[i].twiny() for i in range(3)]
        for i in range(3):
            self.ax[i].tick_params(direction='in')
            # todo Supposedly this is really slow. Check to see if that is still the case
            self.axC[i].format_coord = lambda x, y: ''
            self.axC[i].tick_params(direction='in')

        # spread apart inner and outer plots
        box1 = self.ax[0].get_position()
        self.ax[0].set_position([box1.x0-0.025, box1.y0, box1.width, box1.height])

        box3 = self.ax[2].get_position()
        self.ax[2].set_position([box3.x0+0.025, box3.y0, box3.width, box3.height])

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()

        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        # add the buttons to the toolbar
        self.b = tk.Button(master=self.toolbar, text="Add Experimental Data", command=self.plot_permeation)
        self.b.pack(side="left", padx=1)
        self.b2 = tk.Button(master=self.toolbar, text="Enable Multiple Labels", command=self.toggle_multiple_labels)
        self.b2.pack(side="left", padx=1)

        # Add citation
        b3 = tk.Button(master=self.toolbar, text="Reference", command=self.reference)
        b3.pack(side="left", padx=1)

        # Facilitate an Arrhenius fit to experimental data
        b4 = tk.Button(master=self.toolbar, text="Arrhenius Fit", command=self.fit_to_Arrhenius_eqn)
        b4.pack(side="left", padx=1)

        self.toolbar.update()

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        plot_frame.grid(row=row, column=column, columnspan=2, sticky="nsew")

    def reference(self):
        """ Creates a message box with the overview plots reference inside """
        tk.messagebox.showinfo(title=None, message='M. Shimada, “Tritium Transport in Fusion Reactor Materials” ' +
                               'in Comprehensive Nuclear Materials 6, Second Edition, Elsevier (2020)')

    def plot(self, *args):
        """ Plot the three graphs in overview plots """
        # clear the plot
        for artist in [*self.artists1, *self.artists2]:  # combine into one list
            artist.remove()
        self.artists1, self.artists2 = [], []

        R = self.storage.R  # Gas constant
        # Minimum and maximum temperatures for extrapolating
        Tmin_ext = 373.15
        Tmax_ext = 1473.15
        T2 = np.linspace(Tmin_ext, Tmax_ext)
        for material in self.index:
            if self.item[material].get():  # if the material's check box is checked
                for i in range(3):  # For diffusivity, solubility, and permeability
                    val1 = self.storage.data.loc[material, self.columns[4*i+0]]  # D0, K0, P0
                    val2 = self.storage.data.loc[material, self.columns[4*i+1]]  # E_D, E_K, E_P
                    Tmin = self.storage.data.loc[material, self.columns[4*i+2]]
                    Tmax = self.storage.data.loc[material, self.columns[4*i+3]]

                    T = np.linspace(Tmin, Tmax)
                    res = val1*np.exp(-val2/(R * T))  # data
                    res2 = val1 * np.exp(-val2 / (R * T2))  # extrapolated data

                    ax = self.ax[i]

                    # dashed line for extrapolation, colored for data
                    a, = ax.semilogy(1000 / T2, res2, ':', color="lightgrey", label=material)
                    b, = ax.semilogy(1000 / T, res, color=self.colors[material].get(), label=material)
                    self.artists1.append(a)
                    self.artists2.append(b)
        # add a label to each material
        self.cursor = mplcursors.cursor(self.artists2, hover=mplcursors.HoverMode.Transient, highlight=True)
        self.cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))

        self.plotted_once = True

        def celsius_to_positions(labels):
            # given a temperature in Celsius, return the position of the tick mark to correlate with the K^-1 axis
            positions = 1000/(labels + 273.15)
            return positions

        self.xmin, self.max = 0.0005*1000, 0.0030*1000  # minimum and maximum values for the bottom x-axis

        # Format the x-axes for the diffusivity, solubility, and permeability plots
        for i in range(3):
            self.ax[i].minorticks_off()

            # bottom axis
            self.ax[i].set_xlabel(r"Temperature, 1000/T (K$^{-1}$)")
            self.ax[i].set_xlim(self.xmin, self.max)

            # top axis
            self.axC[i].set_xlabel(u"Temperature, T (\u2103)")
            self.axC[i].set_xlim(self.xmin, self.max)

            # redo the top ticks so that they line up with proper positions for Celsius
            new_ticks_labels = np.array([100, 200, 400, 800, 1600])
            self.axC[i].set_xticks(celsius_to_positions(new_ticks_labels))
            self.axC[i].set_xticklabels(new_ticks_labels)

        # Format the y-axes for diffusivity, solubility, and permeability
        self.ax[0].set_ylabel(r"Diffusivity, $D_H\, ($m$^2\,$s$^{-1})$")
        self.ax[0].set_ylim(1e-12, 1e-7)

        self.ax[1].set_ylabel(r"Solubility, $K_s\, ($mol m$^{-3}\,$Pa$^{-0.5})$")
        self.ax[1].set_ylim(1e-8, 1e+5)

        self.ax[2].set_ylabel(r"Permeability, $P_H\, ($mol m$^{-1}\, $s$^{-1}\, $Pa$^{-0.5})$")
        self.ax[2].set_ylim(1e-17, 1e-5)
        self.canvas.draw()
        self.toolbar.update()

    def cycle_base_colors(self):
        """ keep track of where in the cycle of base colors the cycle is. Reset to 0 when reach the list's end """
        desired_color = self.base_colors[self.current_base_color_idx]
        self.current_base_color_idx += 1

        if self.current_base_color_idx >= len(self.base_colors):
            # reset
            self.current_base_color_idx = 0
        return desired_color

    def change_color(self, material):
        """ opens a color wheel to allow user to specify dataset color """
        color_code = colorchooser.askcolor(title="Choose color")
        if color_code[1] is not None:
            color = color_code[1]
            # control the color of all parts of the button here.
            self.color_buttons[material].configure(bg=color, fg=color,
                                                   highlightbackground=color,
                                                   activebackground=color)
            self.colors[material].set(color_code[1])
            # redraw the plot
            self.plot()

    def selectstart(self, event):
        """ Keep track of the location of the last check button that was toggled in self.start """
        self.start = self.check_buttons.index(event.widget)

    def selectrange(self, event):
        """ toggle all check buttons between self.start and the location of the button just toggled.
            Reset self.start to be the location of the button just toggled. """
        start = self.start
        end = self.check_buttons.index(event.widget)
        sl = slice(min(start, end) + 1, max(start, end))
        for cb in self.check_buttons[sl]:
            cb.toggle()
        self.start = end


class EditMaterials(tk.Toplevel):
    """ popup box used to add, edit, and remove materials in the source spreadsheet. """

    def __init__(self, storage, pos, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # import custom widgets
        widgets = Widgets()
        self.add_entry = widgets.add_entry

        self.title("Add or Edit Materials")
        self.resizable(width=False, height=False)
        self.minsize(400, 170)

        # gui_x/y values determined by running self.updateidletasks() at the end of self.__init__ and then printing size
        # (Note, the window's default height changed from 742 to 750 due to some later edits, but was left as 742 here)
        gui_x = 1429
        gui_y = 742
        if platform.system() == 'Darwin':
            width = 1150
            height = 212
        else:
            width = 870
            height = 170
        pos_right = int(pos[0] + (gui_x-width)/2)
        pos_down = int(pos[1] + (gui_y-height)/3)
        self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        self.storage = storage
        self.changed = False  # store whether any material was added or removed or edited

        # store information
        self.material = tk.StringVar()
        self.tmelt = tk.DoubleVar(value=np.NaN)  # melting temperature

        self.D0 = tk.DoubleVar(value=np.NaN)  # Diffusivity
        self.E_D = tk.DoubleVar(value=np.NaN)  # activation energy for diffusion
        self.tmin_d = tk.DoubleVar(value=np.NaN)  # minimum temperature at which Diffusivity was measured
        self.tmax_d = tk.DoubleVar(value=np.NaN)  # maximum temperature at which Diffusivity was measured
        
        self.K0 = tk.DoubleVar(value=np.NaN)  # Solubility
        self.E_K = tk.DoubleVar(value=np.NaN)
        self.tmin_s = tk.DoubleVar(value=np.NaN)
        self.tmax_s = tk.DoubleVar(value=np.NaN)
        
        self.P0 = tk.DoubleVar(value=np.NaN)  # Permeability
        self.E_P = tk.DoubleVar(value=np.NaN)
        self.tmin_p = tk.DoubleVar(value=np.NaN)
        self.tmax_p = tk.DoubleVar(value=np.NaN)

        # store entry values and Combobox widget
        self.inputs = {}

        # create frames
        self.create_labels()
        self.create_inputs()

        button_frame = tk.Frame(self)
        button_frame.grid(row=0, column=2, sticky="e")

        submit_button = ttk.Button(button_frame, text="Submit", command=self.submit)
        submit_button.pack(side="right")

        quit_button = ttk.Button(button_frame, text="Cancel", command=self.destroy)
        quit_button.pack(side="right")

        remove_button = ttk.Button(button_frame, text="Remove", command=self.remove_material)
        remove_button.pack(side="right")

    def create_labels(self):
        """ Create label/Combobox/text that precedes the entries of the Add or Edit Materials window """
        parent = tk.Frame(self)
        parent.grid(row=0, column=0, columnspan=2, sticky="nsew")

        tk.Label(parent, text="Please enter the material's information.").grid(sticky="nsew", columnspan=2)
        tk.Label(parent, text="Material:").grid(row=1, column=0, sticky="nsw")

        # allow user to add new material or edit an existing one.
        self.inputs['material'] = ttk.Combobox(parent, textvariable=self.material)
        self.inputs['material']['values'] = [""]+list(self.storage.data.index)
        self.inputs['material'].bind('<<ComboboxSelected>>', self.populate_labels)
        self.inputs['material'].grid(row=1, column=1, sticky="nsew")

        self.add_entry(self, parent, variable=self.inputs, key='tmelt', text="Melting Temperature:", subscript='',
                       tvar=self.tmelt, units='[K]', row=1, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.storage.check_for_number(tvar, variable, key,
                                                                                             parent_frame=pf))
        self.inputs["tmelt"].grid(row=1, column=4)

    def create_inputs(self):
        """ Create the entries for editing the properties of the materials """
        parent = tk.Frame(self)
        parent.grid(row=1, column=0, columnspan=3, sticky="nsew")
        tk.Label(parent, text="Diffusivity").grid(row=0, column=1)
        self.add_entry(self, parent, variable=self.inputs, key="D0", text="D", innersubscript="0", innertext=":",
                       subscript="", tvar=self.D0, units="[m\u00b2 s\u207b\u00b9]", row=1,  # [m^2 s^-1]
                       in_window=True, command=lambda tvar, variable, key, pf: self.calc_DKP(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="E_D", text="E", innersubscript="D", innertext=":",
                       subscript="", tvar=self.E_D, units="[kJ mol\u207b\u00b9]", row=2,  # [kJ mol^-1]
                       in_window=True, command=lambda tvar, variable, key, pf: self.calc_E(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmin_d", text="min temp:",
                       subscript="", tvar=self.tmin_d, units="[K]", row=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmin(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmax_d", text="max temp:",
                       subscript="", tvar=self.tmax_d, units="[K]", row=4, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmax(tvar, variable, key, pf))

        tk.Label(parent, text="Solubility").grid(row=0, column=4)
        self.add_entry(self, parent, variable=self.inputs, key="K0", text="K", innersubscript="0", innertext=":",
                       subscript="", tvar=self.K0,
                       units=u"[mol m\u207b\u00B3 Pa\u207b\u2070\u1427\u2075]",  # mol m^-3 Pa^-0.5
                       row=1, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_DKP(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="E_K", text="E", innersubscript="K", innertext=":",
                       subscript="", tvar=self.E_K, units="[kJ mol\u207b\u00b9]",  # [kJ mol^-1]
                       row=2, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_E(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmin_s", text="min temp:",
                       subscript="", tvar=self.tmin_s, units="[K]", row=3, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmin(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmax_s", text="max temp:",
                       subscript="", tvar=self.tmax_s, units="[K]", row=4, column=3, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmax(tvar, variable, key, pf))

        tk.Label(parent, text="Permeability").grid(row=0, column=7)
        self.add_entry(self, parent, variable=self.inputs, key="P0", text="P", innersubscript="0", innertext=":",
                       subscript="", tvar=self.P0,  # units=[mol m^-1 s^-1 Pa^-0.5]
                       units=u"[mol m\u207b\u00b9 s\u207b\u00b9 Pa\u207b\u2070\u1427\u2075]",
                       row=1, column=6, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_DKP(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="E_P", text="E", innersubscript="P", innertext=":",
                       subscript="", tvar=self.E_P, units="[kJ mol\u207b\u00b9]",  # [kJ mol^-1]
                       row=2, column=6, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_E(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmin_p", text="min temp:",
                       subscript="", tvar=self.tmin_p, units="[K]", row=3, column=6, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmin(tvar, variable, key, pf))
        self.add_entry(self, parent, variable=self.inputs, key="tmax_p", text="max temp:",
                       subscript="", tvar=self.tmax_p, units="[K]", row=4, column=6, in_window=True,
                       command=lambda tvar, variable, key, pf: self.calc_tmax(tvar, variable, key, pf))

    def populate_labels(self, *args):
        """ load the data based on user choice of material """
        material = self.material.get()
        if material != "":  # If an actual material was selected...
            # disable editing in combobox
            self.inputs['material']['state'] = 'readonly'

            # fill in entry boxes
            '''  Deleting the current contents of of the entry box and updating it also updates 
                 self.inputs['tmelt' + 'var'], the variable that contains the entry box's string. However, it does this
                 without calling the verification function, which is why self.tmelt needs to get updated separately.'''
            self.inputs['tmelt'].delete(0, "end")
            self.inputs['tmelt'].insert(0, "{}".format(self.storage.data.loc[
                                                           material, ("melting temp. [K]", "")]))
            self.tmelt.set(self.storage.data.loc[material, ("melting temp. [K]", "")])

            # Diffusivity
            self.inputs['D0'].delete(0, "end")
            self.inputs['D0'].insert(0, "{}".format(self.storage.data.loc[
                                                        material, ("Diffusivity", "D0 [m^2 s^-1]")]))
            self.D0.set(self.storage.data.loc[material, ("Diffusivity", "D0 [m^2 s^-1]")])

            self.inputs['E_D'].delete(0, "end")
            self.inputs['E_D'].insert(0, "{}".format(self.storage.data.loc[
                                                         material, ("Diffusivity", "E_D [kJ mol^-1]")]))
            self.E_D.set(self.storage.data.loc[material, ("Diffusivity", "E_D [kJ mol^-1]")])

            self.inputs['tmin_d'].delete(0, "end")
            self.inputs['tmin_d'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Diffusivity", "min. temp. [K]")]))
            self.tmin_d.set(self.storage.data.loc[material, ("Diffusivity", "min. temp. [K]")])

            self.inputs['tmax_d'].delete(0, "end")
            self.inputs['tmax_d'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Diffusivity", "max. temp. [K]")]))
            self.tmax_d.set(self.storage.data.loc[material, ("Diffusivity", "max. temp. [K]")])

            # Solubility
            self.inputs['K0'].delete(0, "end")
            self.inputs['K0'].insert(0, "{}".format(self.storage.data.loc[
                                                        material, ("Solubility", "K0 [mol m^-3 Pa^-0.5]")]))
            self.K0.set(self.storage.data.loc[material, ("Solubility", "K0 [mol m^-3 Pa^-0.5]")])

            self.inputs['E_K'].delete(0, "end")
            self.inputs['E_K'].insert(0, "{}".format(self.storage.data.loc[
                                                         material, ("Solubility", "E_K [kJ mol^-1]")]))
            self.E_K.set(self.storage.data.loc[material, ("Solubility", "E_K [kJ mol^-1]")])

            self.inputs['tmin_s'].delete(0, "end")
            self.inputs['tmin_s'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Solubility", "min. temp. [K]")]))
            self.tmin_s.set(self.storage.data.loc[material, ("Solubility", "min. temp. [K]")])

            self.inputs['tmax_s'].delete(0, "end")
            self.inputs['tmax_s'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Solubility", "max. temp. [K]")]))
            self.tmax_s.set(self.storage.data.loc[material, ("Solubility", "max. temp. [K]")])

            # Permeability
            self.inputs['P0'].delete(0, "end")
            self.inputs['P0'].insert(0, "{}".format(self.storage.data.loc[
                                                        material, ("Permeability", "P0 [mol m^-1 s^-1 Pa^-0.5]")]))
            self.P0.set(self.storage.data.loc[material, ("Permeability", "P0 [mol m^-1 s^-1 Pa^-0.5]")])

            self.inputs['E_P'].delete(0, "end")
            self.inputs['E_P'].insert(0, "{}".format(self.storage.data.loc[
                                                         material, ("Permeability", "E_P [kJ mol^-1]")]))
            self.E_P.set(self.storage.data.loc[material, ("Permeability", "E_P [kJ mol^-1]")])

            self.inputs['tmin_p'].delete(0, "end")
            self.inputs['tmin_p'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Permeability", "min. temp. [K]")]))
            self.tmin_p.set(self.storage.data.loc[material, ("Permeability", "min. temp. [K]")])

            self.inputs['tmax_p'].delete(0, "end")
            self.inputs['tmax_p'].insert(0, "{}".format(self.storage.data.loc[
                                                            material, ("Permeability", "max. temp. [K]")]))
            self.tmax_p.set(self.storage.data.loc[material, ("Permeability", "max. temp. [K]")])

        else:
            # enable editing
            self.inputs['material']['state'] = 'normal'
            # clear entries
            ''' Deleting the current contents of of the entry box and updating it also updates 
                self.inputs['tmelt' + 'var'], the variable that contains the entry box's string. However, it does this
                without calling the verification function, which is why self.tmelt needs to get updated separately.'''
            self.inputs['tmelt'].delete(0, "end")
            self.inputs['tmelt'].insert(0, "{}".format(np.NaN))
            self.tmelt.set(np.NaN)

            # Diffusivity
            self.inputs['D0'].delete(0, "end")
            self.inputs['D0'].insert(0, "{}".format(np.NaN))
            self.D0.set(np.NaN)

            self.inputs['E_D'].delete(0, "end")
            self.inputs['E_D'].insert(0, "{}".format(np.NaN))
            self.E_D.set(np.NaN)

            self.inputs['tmin_d'].delete(0, "end")
            self.inputs['tmin_d'].insert(0, "{}".format(np.NaN))
            self.tmin_d.set(np.NaN)

            self.inputs['tmax_d'].delete(0, "end")
            self.inputs['tmax_d'].insert(0, "{}".format(np.NaN))
            self.tmax_d.set(np.NaN)

            # Solubility
            self.inputs['K0'].delete(0, "end")
            self.inputs['K0'].insert(0, "{}".format(np.NaN))
            self.K0.set(np.NaN)

            self.inputs['E_K'].delete(0, "end")
            self.inputs['E_K'].insert(0, "{}".format(np.NaN))
            self.E_K.set(np.NaN)

            self.inputs['tmin_s'].delete(0, "end")
            self.inputs['tmin_s'].insert(0, "{}".format(np.NaN))
            self.tmin_s.set(np.NaN)

            self.inputs['tmax_s'].delete(0, "end")
            self.inputs['tmax_s'].insert(0, "{}".format(np.NaN))
            self.tmax_s.set(np.NaN)

            # Permeability
            self.inputs['P0'].delete(0, "end")
            self.inputs['P0'].insert(0, "{}".format(np.NaN))
            self.P0.set(np.NaN)

            self.inputs['E_P'].delete(0, "end")
            self.inputs['E_P'].insert(0, "{}".format(np.NaN))
            self.E_P.set(np.NaN)

            self.inputs['tmin_p'].delete(0, "end")
            self.inputs['tmin_p'].insert(0, "{}".format(np.NaN))
            self.tmin_p.set(np.NaN)

            self.inputs['tmax_p'].delete(0, "end")
            self.inputs['tmax_p'].insert(0, "{}".format(np.NaN))
            self.tmax_p.set(np.NaN)

    def calc_DKP(self, tvar, variable, key, pf=None):
        """ Checks if user entry is a number and, if a number, checks to see if exactly two out of the three
            DKP values are entered in. If so, calculates and fills the third value. If the user entry is not a
            number, reverts entry to previous value """
        # Ensure the input is a number
        self.storage.check_for_number(tvar, variable, key, parent_frame=pf)

        # If two out of the three variables are submitted, derive the third
        D0, K0, P0 = np.isnan(self.D0.get()), np.isnan(self.K0.get()), np.isnan(self.P0.get())
        if (int(D0) + int(K0) + int(P0)) == 1:
            if D0:  # If D0 is nan...
                # Update all three parts related to the entry because various attempts to just update two didn't work
                variable['D0' + 'var'].set(value=round(self.P0.get()/self.K0.get(), 13))
                self.inputs['D0'].delete(0, "end")
                self.inputs['D0'].insert(0, "{}".format(round(self.P0.get()/self.K0.get(), 13)))
                self.D0.set(round(self.P0.get() / self.K0.get(), 13))
            elif K0:
                variable['K0' + 'var'].set(value=round(self.P0.get() / self.D0.get(), 13))
                self.inputs['K0'].delete(0, "end")
                self.inputs['K0'].insert(0, "{}".format(round(self.P0.get() / self.D0.get(), 13)))
                self.K0.set(round(self.P0.get() / self.D0.get(), 13))
            elif P0:
                variable['P0' + 'var'].set(value=round(self.D0.get() * self.K0.get(), 13))
                self.inputs['P0'].delete(0, "end")
                self.inputs['P0'].insert(0, "{}".format(round(self.D0.get() * self.K0.get(), 13)))
                self.P0.set(round(self.D0.get() * self.K0.get(), 13))
        return True

    def calc_E(self, tvar, variable, key, pf=None):
        """ Checks if user entry is a number and, if a number, checks to see if exactly two out of the three
            energy values are entered in. If so, calculates and fills the third value. If the user entry is not a
            number, reverts entry to previous value """
        # Ensure the input is a number
        self.storage.check_for_number(tvar, variable, key, parent_frame=pf)

        # If two out of the three variables are submitted, derive the third
        E_D, E_K, E_P = np.isnan(self.E_D.get()), np.isnan(self.E_K.get()), np.isnan(self.E_P.get())
        if int(E_D) + int(E_K) + int(E_P) == 1:
            if E_D:  # If E_D is nan...
                # Update all three parts related to the entry because various attempts to just update two didn't work
                variable['E_D' + 'var'].set(value=round(self.E_P.get() - self.E_K.get(), 13))
                self.inputs['E_D'].delete(0, "end")
                self.inputs['E_D'].insert(0, "{}".format(round(self.E_P.get() - self.E_K.get(), 13)))
                self.E_D.set(round(self.E_P.get() - self.E_K.get(), 13))
            elif E_K:
                variable['E_K' + 'var'].set(value=round(self.E_P.get() - self.E_D.get(), 13))
                self.inputs['E_K'].delete(0, "end")
                self.inputs['E_K'].insert(0, "{}".format(round(self.E_P.get() - self.E_D.get(), 13)))
                self.E_K.set(round(self.E_P.get() - self.E_D.get(), 13))
            elif E_P:
                variable['E_P' + 'var'].set(value=round(self.E_D.get() + self.E_K.get(), 13))
                self.inputs['E_P'].delete(0, "end")
                self.inputs['E_P'].insert(0, "{}".format(round(self.E_D.get() + self.E_K.get(), 13)))
                self.E_P.set(round(self.E_D.get() + self.E_K.get(), 13))
        return True

    def calc_tmin(self, tvar, variable, key, pf=None):
        """ Checks if user entry is a number and, if a number, checks to see if exactly two out of the three
            temperature values are entered in. If so, calculates and fills the third value. If the user entry is not a
            number, reverts entry to previous value """
        # Ensure the input is a number
        self.storage.check_for_number(tvar, variable, key, parent_frame=pf)

        # If two out of the three variables are submitted, derive the third
        tmin_d, tmin_s, tmin_p = np.isnan(self.tmin_d.get()), np.isnan(self.tmin_s.get()), np.isnan(self.tmin_p.get())
        if int(tmin_d) + int(tmin_s) + int(tmin_p) == 1:
            if tmin_d:  # If tmin_d is nan...
                # Update all three parts related to the entry because various attempts to just update two didn't work
                variable['tmin_d' + 'var'].set(value=max((self.tmin_p.get(), self.tmin_s.get())))
                self.inputs['tmin_d'].delete(0, "end")
                self.inputs['tmin_d'].insert(0, "{}".format(max((self.tmin_p.get(), self.tmin_s.get()))))
                self.tmin_d.set(max((self.tmin_p.get(), self.tmin_s.get())))
            elif tmin_s:
                variable['tmin_s' + 'var'].set(value=max((self.tmin_p.get(), self.tmin_d.get())))
                self.inputs['tmin_s'].delete(0, "end")
                self.inputs['tmin_s'].insert(0, "{}".format(max((self.tmin_p.get(), self.tmin_d.get()))))
                self.tmin_s.set(max((self.tmin_p.get(), self.tmin_d.get())))
            elif tmin_p:
                variable['tmin_p' + 'var'].set(value=max((self.tmin_d.get(), self.tmin_s.get())))
                self.inputs['tmin_p'].delete(0, "end")
                self.inputs['tmin_p'].insert(0, "{}".format(max((self.tmin_d.get(), self.tmin_s.get()))))
                self.tmin_p.set(max((self.tmin_d.get(), self.tmin_s.get())))
        return True

    def calc_tmax(self, tvar, variable, key, pf=None):
        """ Checks if user entry is a number and, if a number, checks to see if exactly two out of the three
            temperature values are entered in. If so, calculates and fills the third value. If the user entry is not a
            number, reverts entry to previous value """
        # Ensure the input is a number
        self.storage.check_for_number(tvar, variable, key, parent_frame=pf)

        # If two out of the three variables are submitted, derive the third
        tmax_d, tmax_s, tmax_p = np.isnan(self.tmax_d.get()), np.isnan(self.tmax_s.get()), np.isnan(self.tmax_p.get())
        if int(tmax_d) + int(tmax_s) + int(tmax_p) == 1:
            if tmax_d:  # If tmax_d is nan...
                # Update all three parts related to the entry because various attempts to just update two didn't work
                variable['tmax_d' + 'var'].set(value=min((self.tmax_p.get(), self.tmax_s.get())))
                self.inputs['tmax_d'].delete(0, "end")
                self.inputs['tmax_d'].insert(0, "{}".format(min((self.tmax_p.get(), self.tmax_s.get()))))
                self.tmax_d.set(min((self.tmax_p.get(), self.tmax_s.get())))
            elif tmax_s:
                variable['tmax_s' + 'var'].set(value=min((self.tmax_p.get(), self.tmax_d.get())))
                self.inputs['tmax_s'].delete(0, "end")
                self.inputs['tmax_s'].insert(0, "{}".format(min((self.tmax_p.get(), self.tmax_d.get()))))
                self.tmax_s.set(min((self.tmax_p.get(), self.tmax_d.get())))
            elif tmax_p:
                variable['tmax_p' + 'var'].set(value=min((self.tmax_d.get(), self.tmax_s.get())))
                self.inputs['tmax_p'].delete(0, "end")
                self.inputs['tmax_p'].insert(0, "{}".format(min((self.tmax_d.get(), self.tmax_s.get()))))
                self.tmax_p.set(min((self.tmax_d.get(), self.tmax_s.get())))
        return True

    def submit(self):
        """ Handles potential new materials or edits to old materials"""
        # check that all entry boxes are filled (allow ing an exception for the melting temp)
        if sum(map(np.isnan, [self.D0.get(), self.E_D.get(), self.tmin_d.get(), self.tmax_d.get(),
                              self.K0.get(), self.E_K.get(), self.tmin_s.get(), self.tmax_s.get(),
                              self.P0.get(), self.E_P.get(), self.tmin_p.get(), self.tmax_p.get()])) > 0:
            tk.messagebox.showwarning(title="Missing Entry",
                                      message="Please fill in remaining boxes before continuing.")
            self.deiconify()
        else:  # if all entries have been filled in...
            # popup message dependent on if it's an edit or an addition
            if self.material.get() in self.storage.data.index:
                c_message = "Edit the recorded properties of {} in the source file for the materials?".format(
                    self.material.get()) + \
                            '\n \u2022 Selecting "yes" will change them here and in the original Excel sheet.' + \
                            '\n \u2022 Selecting "no" will change them here but not in the original spreadsheet.' + \
                            '\n \u2022 Selecting "cancel" will return you to the Add or Edit Materials window.'
            else:
                c_message = "Add {} to the source file for the materials?".format(self.material.get()) + \
                            '\n \u2022 Selecting "yes" will add it here and in the original Excel sheet.' + \
                            '\n \u2022 Selecting "no" will add it here but not in the original spreadsheet.' + \
                            '\n \u2022 Selecting "cancel" will return you to the Add or Edit Materials window.'
            ans = tk.messagebox.askyesnocancel(
                title="Confirmation",
                message=c_message)
            if ans:
                # Change the selected material to have the new properties (or create that material if it didn't exist)
                self.storage.data.loc[
                    self.material.get()] = [self.D0.get(), self.E_D.get(), self.tmin_d.get(), self.tmax_d.get(),
                                            self.K0.get(), self.E_K.get(), self.tmin_s.get(), self.tmax_s.get(),
                                            self.P0.get(), self.E_P.get(), self.tmin_p.get(), self.tmax_p.get(),
                                            self.tmelt.get()]

                # Change the material dataframe reserved for editing the source file
                self.storage.material_data.loc[
                    self.material.get()] = [self.D0.get(), self.E_D.get(), self.tmin_d.get(), self.tmax_d.get(),
                                            self.K0.get(), self.E_K.get(), self.tmin_s.get(), self.tmax_s.get()]

                # overwrite the old source file with the updated dataframe
                material_filename = os.path.join('datafiles', 'material_data.xlsx')
                self.storage.material_data.to_excel(material_filename)
                # This next bit is needed because to_excel formats named indexes poorly
                import openpyxl
                material_data_file = openpyxl.load_workbook(material_filename)
                material_data_file["Sheet1"].delete_rows(3)
                material_data_file["Sheet1"]['A2'] = "Material"
                material_data_file.save(material_filename)

                # Check if a melting temperature was provided. If so, add it to the dataframe reserved for updating
                # the melting temperatures source file, then update that source file.
                if not np.isnan(self.tmelt.get()):
                    self.storage.melting_tempK.loc[self.material.get()] = [self.tmelt.get()]
                    self.storage.melting_tempK.to_excel(os.path.join('datafiles', 'melting_tempK.xlsx'))

                self.changed = True
                self.destroy()
            elif ans is False:
                # Change the selected material to have the new properties (or create that material if it didn't exist)
                self.storage.data.loc[
                    self.material.get()] = [self.D0.get(), self.E_D.get(), self.tmin_d.get(), self.tmax_d.get(),
                                            self.K0.get(), self.E_K.get(), self.tmin_s.get(), self.tmax_s.get(),
                                            self.P0.get(), self.E_P.get(), self.tmin_p.get(), self.tmax_p.get(),
                                            self.tmelt.get()]
                self.changed = True
                self.destroy()
            else:
                self.deiconify()  # bring back the pop-up window

    def remove_material(self):
        """ Handles removal of old materials """
        ans = tk.messagebox.askyesnocancel(
            title="Confirmation",
            message='Do you want to also remove {} from the source file for materials?'.format(self.material.get()) +
                    '\n \u2022 Selecting "yes" will remove it from here and the original Excel sheet.' +
                    '\n \u2022 Selecting "no" will remove it from here but not the original spreadsheet.' +
                    '\n \u2022 Selecting "cancel" will return you to the Add or Edit Materials window.')

        if ans:
            try:
                # Remove the material from the dataframe that influences the program (storage.data) and from the
                # dataframe that will be used to replace the source file (which only gets edited when if specifically
                # selected to be edited)
                self.storage.data.drop(index=self.material.get(), inplace=True)
                if self.material.get() in self.storage.material_data.index:  # Check the material's in this dataframe
                    self.storage.material_data.drop(index=self.material.get(), inplace=True)

                    # overwrite the old source file with the updated dataframe
                    material_filename = os.path.join('datafiles', 'material_data.xlsx')
                    self.storage.material_data.to_excel(material_filename)
                    # This next bit is needed because to_excel formats named indexes poorly
                    import openpyxl
                    material_data_file = openpyxl.load_workbook(material_filename)
                    material_data_file["Sheet1"].delete_rows(3)
                    material_data_file["Sheet1"]['A2'] = "Material"
                    material_data_file.save(material_filename)

                    # Check if a melting temperature was provided. If so, remove it from the dataframe reserved for
                    # updating the melting temperatures source file, then update that source file.
                    if self.material.get() in self.storage.melting_tempK.index:
                        self.storage.melting_tempK.drop(index=self.material.get(), inplace=True)
                        self.storage.melting_tempK.to_excel(os.path.join('datafiles', 'melting_tempK.xlsx'))

                self.changed = True
                self.destroy()  # Closes Add/Edit Materials window
            except KeyError:
                # If no material has been selected
                tk.messagebox.showwarning(title="Missing Selection",
                                          message="Please select a material to remove.")
                self.deiconify()
        elif ans is False:
            try:
                # Remove the material from the dataframe that influences the program (storage.data)
                self.storage.data.drop(index=self.material.get(), inplace=True)
                self.changed = True
                self.destroy()  # Closes Add/Edit Materials window
            except KeyError:
                # If no material has been selected
                tk.messagebox.showwarning(title="Missing Selection",
                                          message="Please select a material to remove.")
                self.deiconify()
        else:
            self.deiconify()  # bring back the pop-up window

    def show_aem_window(self):
        """ This allows the window to open and then, when closed, trigger the main gui to update the materials if
            any changes were made """
        self.deiconify()
        self.wait_window()

        return self.changed, self.material.get()
