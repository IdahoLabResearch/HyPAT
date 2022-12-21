# Hydrogen Permeation and Absorption Tool (HyPAT v2.0.0)
## About 
The Hydrogen Permeation and Absorption Tool (HyPAT v2.0.0)  is an application that streamlines data analysis for hydrogen 
permeation and absorption experiments to measure the following hydrogen transport properties: permeability, diffusivity, 
and solubility. Specifically, HyPAT analyzes data from gas-driven permeation through metals using the build-up in closed 
volume method and data from absorption experiments using a Sieverts'-type apparatus. A built-in literature database 
provides direct comparison for experimental results as well as timescale and permeation rate estimates for experiment 
preparation based on the user-input parameters. 

Find more details about HyPAT v1.0.0 in the SoftwareX publication ["HyPAT: A GUI for high-throughput gas-driven hydrogen permeation data analysis."](https://doi.org/10.1016/j.softx.2022.101284) A publication on HyPAT v2.0.0 is in progress.

## Description
The Hydrogen Permeation and Absorption Tool (HyPAT) provides a convenient application that allows the user to 
prepare permeation experiment estimates, quickly analyze multiple data files, and compare results to literature values.

Features:

* Calculates the hydrogen transport properties of permeability, diffusivity, and solubility from experimental pressure-rise permeation data.
* Calculates the hydrogen transport properties of permeability, diffusivity, and solubility from experimental absorption data obtained from a Sieverts'-type apparatus.
* Batchwise data analysis of multiple files.
* Direct comparison of measured hydrogen transport properties to known literature values.
* Predicts the permeation rates and diffusion times of permeation experiments.
* Easily customizable to adapt to user needs and experimental setups.

## Examples of usage of HyPAT 

###  Usage: Data Analysis with "Permeation Plots" Tab

This tab calculates the hydrogen transport properties permeability, diffusivity, and solubility from experimental 
permeation data using pressure build-up in a secondary side closed volume while allowing for batchwise data analysis of 
multiple files. The required user inputs and associated errors are highlighted in yellow. Thickness (mm) is the only 
required sample property. The experimental apparatus parameters are the secondary side volume (m<sup>3</sup>) and 
permeation surface area (m<sup>2</sup>). 

Clicking the "Steady State Variables" button beneath the top-right plot allows the user to edit the 
tolerance for determining steady state, the delay until steady state, and the number of points used to find the leak 
rate and steady state rate. The tolerance is the value at which the second derivative of the secondary side pressure 
with respect to time (∂<sup>2</sup>P/∂t<sup>2</sup>) is approximated to be zero. The delay until steady state is the 
minimum number of seconds to wait after the isolation valve is opened until initiating analysis to detect steady state. 
The leak range is the number of data points prior to opening the isolation valve used to determine the leak rate and 
initial values. The steady state range is the number of data points used to determine final values and when the 
permeation rate has reached a steady state. 

Click the "Choose folder" button to select the folder containing experimental data. The graphs will then self-populate. 
Note: Only XLS and XLSX files are processed by the program; all other files are ignored. HyPAT assumes data is organized 
such that the data of each instrument (or set of instruments) has its own column in the Excel sheet without a header.

The bottom left graph displays the permeability calculated from each data set plotted against temperature. Use the 
"Current Measurement" drop-down menu to change to the diffusivity, solubility, or flux calculated from each file plotted 
against temperature (or against pressure, in the case of flux).

Also visible are a pressure versus time plot, a permeability versus time plot, and a comparison of different diffusivity 
optimizations plot, all using data from the first file in the selected folder. Use the "Current File" drop-down menu to 
use data from a different file for these plots.

By clicking the "Settings" button, you gain the ability to edit some of how the program processes the data. 
Specifically, the following parameters are editable:

* The columns of the Excel sheet from which the program reads data for each instrument,
* The row at which the program starts reading data from the Excel sheet,
* The number of rows at the end of the Excel sheet to ignore when reading data,
* The number of thermocouples used to measure the gas temperature,
* Conversion factors for each instrument to facilitate appropriate unit conversions,
* Constant and proportional errors of each instrument.

###  Usage: Data Analysis with "Absorption Plots" Tab

The "Absorption Plots" tab calculates the hydrogen transport properties permeability, 
diffusivity, and solubility from experimental absorption data obtained using a Sieverts'-type apparatus. HyPAT does this 
while allowing for batchwise data analysis of multiple files. The required user inputs and associated errors are 
highlighted in yellow. Thickness (mm), mass (g), volume (cm<sup>3</sup>), and molar mass (g mol<sup>-1</sup>) are the 
required sample properties. The required experimental apparatus parameters are the initial container volume 
(cm<sup>3</sup>) and the sample container volume (cm<sup>3</sup>). HyPAT also requires the user to submit whether the 
experiment type is "Single," meaning each data file has the sample measured at only one pressure, or "Isotherm," 
meaning at least one data file has the sample measured at multiple pressures along an isotherm.

Clicking the "Equilibrium Variables" button immediately beneath the top-right plot allows the user to edit the 
tolerance for determining equilibrium, the delay until equilibrium, the number of points used to find the equilibrium, 
and the number of data points used to find initial and equilibrium values. The tolerance is the value at which the 
second derivative of the pressure with respect to time (∂<sup>2</sup>P/∂t<sup>2</sup>) is approximated to be zero. 
The delay until equilibrium is the minimum number of seconds to wait after the isolation valve is opened until 
initiating analysis to detect equilibrium. The initial values range is the number of data points prior to opening the 
isolation valve used to determine the initial values. The equilibrium range is the number of data points used to 
determine final values and when the absorption rate has reached an equilibrium. 

Click the "Choose folder" button to select the folder containing experimental data. The graphs will then self-populate. 
Note: Only XLS and XLSX files are processed by the program; all other files are ignored. HyPAT assumes data is organized 
such that the data of each instrument (or set of instruments) has its own column in the Excel sheet without a header.

The bottom left graph displays the solubility calculated from each data set plotted against temperature. Use the 
"Current Measurement" drop-down menu to change to the diffusivity or permeability calculated from each file plotted 
against temperature. The bottom middle graph displays the final pressure of each data set plotted against the final 
composition. The points are color coded according to final sample temperature, yielding a 
pressure-composition-temperature (PCT) plot.

Also visible are a raw data plot and a comparison of different diffusivity optimizations plot, each using data from the 
first file in the selected folder. Use the "Current File" drop-down menu to use data from a different file for the raw 
data and diffusivity comparison plots.

By clicking the "Settings" button, you gain the ability to edit some of how the program processes the data. 
Specifically, the following parameters are editable:

* The columns of the Excel sheet from which the program reads data for each instrument,
* The row at which the program starts reading data from the Excel sheet,
* The number of rows at the end of the Excel sheet to ignore when reading data,
* The number of thermocouples used to measure the gas temperature initially and after the isolation valve is opened,
* Conversion factors for each instrument to facilitate appropriate unit conversions,
* Constant and proportional errors of each instrument.

### Usage: Comparison with Literature Values with "Overview Plots" Tab

The "Overview Plots" tab allows quick comparison to literature values of diffusivity, solubility, and permeability 
taken from M. Shimada (2020). The user can select and remove materials displayed in these overview plots using the 
materials column. Materials can be added to, modified, or removed from the database with the "Add/Edit Material" button.

The button labelled "Arrhenius Fit" calculates the pre-exponential factors and activation energies of diffusivity, 
solubility, and permeability for a material.

Reference:

M. Shimada, "Tritium Transport in Fusion Reactor Materials" in Comprehensive Nuclear Materials 6, Second Edition, 
Elsevier (2020); https://doi.org/10.1016/B978-0-12-803581-8.11754-0

### Usage: Predict Outcomes with "Permeation Estimates" Tab

The "Permeation Estimates" tab provides predictions for permeation experiments based on specified user inputs. The required user inputs are 
highlighted in yellow. The sample properties are thickness (mm) and material. The experimental apparatus parameters are 
O-ring/sealing method, calibrated leak rate (mol s<sup>-1</sup> Torr<sup>-1</sup>), secondary side volume (cc), and 
accumulation time (hr). The specific test parameters are sample temperature (°C) and primary calculated partial pressure 
(Torr).

The results reported in the "FINAL OUTPUT" section include estimated time-lag (s), permeation rate 
(mol s<sup>-1</sup>), and final secondary pressure (Torr).

By clicking the "Settings" button, O-rings can be added, edited, or removed. The default values for calibrated leak rate 
and secondary side volume can also be edited there. The options for sample material can be edited using the "Add/Edit 
Material" button on the "Overview Plots" tab.

## Installation Instructions

1.	Select the button "Code" near the top of this directory.
2.	Select "Download ZIP" from the dropdown menu that appears.
3.	Extract the downloaded file to a folder of your choice.
4.	Open the command line, then move to the directory "HyPAT" within the downloaded folder.
5.	Run the following code to open the application: python main.py
 
* Note:
  * Your system may require the following code instead: python3 main.py
  * These instructions assume you already have Python and the required Python packages installed.

Python Interpreter: 3.8, 3.9, or 3.10

Python Packages: Matplotlib, Pandas, NumPy, tkmacosx, mplcursors, SciPy, and openpyxl

Compatible: macOS and Windows 10

## Plans

There are no current updates planned for HyPAT. 

## Contribution Guidelines
 
If you notice a feature that you want or a bug that needs to be fixed, please contact Chase Taylor or Thomas Fuerst.
 
## Authors and Acknowledgement 

George Evans

Joseph Watkins

Thomas Fuerst

Chase Taylor 

Joey Watkins (https://github.com/joeymwatkins) built most of the original functionality. George Evans 
(https://github.com/GeorgeEvans0) added the "Absorption Plots" tab, cleaned up errors, made visual changes, and 
added functionality. Thomas Fuerst (https://github.com/FuerstT) was the primary supervisor over this project, with Chase 
Taylor (https://github.com/taylchas) as secondary supervisor.

This work was prepared for the U.S. Department of Energy, Office of Fusion Energy Sciences, under the DOE Idaho Field 
Office contract number DE-AC07–05ID14517.

## License

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an 
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the file "LICENSE" for the 
specific language governing permissions and limitations under the License.
