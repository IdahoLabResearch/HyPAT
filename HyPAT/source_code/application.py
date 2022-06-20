""" Top Level containing the appearance of the GUI
This module creates the base container for the program and loads the notebook + tabs.
"""
import tkinter as tk
from tkinter import ttk, font
from . import data_storage
from . import permeation_estimates
from . import overview_plots
from . import permeation_plots


class Application(tk.Tk):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Change the default font for text boxes to TkDefaultFont. This mostly affects
        # the help window of the permeation_plot.py's settings, but may affect other places
        default_font = font.nametofont("TkFixedFont")
        default_font.configure(family=font.nametofont("TkDefaultFont").cget("family"),
                               size=font.nametofont("TkDefaultFont").cget("size"))

        self.title("Hydrogen Permeation Analysis Tool")
        self.resizable(width=None, height=None)
        self.minsize(900, 600)

        # size of gui based of size of widgets
        # this could be calculated later but would result in the gui appearing then moving to the center.
        width = 1429
        height = 750
        # position of top left corner of gui based on monitor resolution
        pos_right = int(self.winfo_screenwidth() / 2 - width / 2)
        pos_down = int(self.winfo_screenheight() / 2 - height / 2)
        # set position
        self.geometry("{}x{}+{}+{}".format(width, height, pos_right, pos_down))

        # Ensure the program stops when the window is closed (counteracts a bug brought in from using
        # constrained_layouts = True in the permeation_plots.py code)
        self.protocol("WM_DELETE_WINDOW", self.quit_me)

        self.storage = data_storage.Storage()
        self.notebook = ttk.Notebook()
        self.load_pages()
        self.notebook.grid(row=0, sticky='nsew')

        top = self.winfo_toplevel()
        top.rowconfigure(0, weight=1)  # make row zero of top level window stretchable
        top.columnconfigure(0, weight=1)  # make column zero of top level window stretchable

    def quit_me(self):
        """ Ensure the program ends when window is closed """
        # This counteracts a bug brought in by constrained_layout=True in the permeation_plots code
        # Function found at https://stackoverflow.com/a/55206851
        self.quit()
        self.destroy()

    def load_pages(self):
        """ Load the GUI with the pages containing the calculations """
        front_page = permeation_estimates.InputForm(self.notebook, self.storage)  # Permeation estimates
        experiment_plots_page = permeation_plots.PermeationPlots(self.notebook, self.storage)  # Permeation plots
        overview_plots_page = overview_plots.Plots(self.notebook, self.storage)  # Theoretical plots (Overview)
        self.notebook.add(front_page, text="Permeation Estimates")
        self.notebook.add(experiment_plots_page, text="Permeation Plots")
        self.notebook.add(overview_plots_page, text="Overview Plots")
