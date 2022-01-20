# Hydrogen Permeation Analysis Tool (HyPAT) 
## About 
The Hydrogen Permeation Analysis Tool (HyPAT) is an application that streamlines data analysis for gas-driven permeation through metals using the buildup in closed volume method to measure the following hydrogen transport properties: permeability, diffusivity, and solubility. A built-in literature database provides direct comparison for experimental results as well as timescale and permeation rate estimates for experiment preparation based on the user-input parameters. 
## Description
The Hydrogen Permeation Analysis Tool (HyPAT) provides a convenient application that allows the user to prepare experiment estimates, quickly analyze multiple data files, and compare results to literature values based on gas-driven permeation through metals using the buildup method.

Features:

* Calculates the hydrogen transport properties of permeability, diffusivity, and solubility from experimental pressure-rise permeation data.
* Batchwise data analysis of multiple files.
* Direct comparison of measured hydrogen transport properties to known literature values.
* Predicts the permeation rates and diffusion times of experiments.
* Easily customizable to adapt to user needs and experimental setups.

## Examples of usage of HyPAT 

###  Usage: Data Analysis with “Permeation Plots" Tab

This tab calculates the hydrogen transport properties permeability, diffusivity, and solubility from experimental permeation data using pressure buildup in a secondary side closed volume while allowing for batchwise data analysis of multiple files. The required user inputs and associated errors are highlighted in yellow. Thickness (mm) is the only required sample property. The experimental apparatus parameters are the secondary side volume (m<sup>3</sup>) and permeation surface area (m<sup>2</sup>). The tolerance for determining a steady state is also adjustable. The tolerance is the value at which the second derivative of the secondary pressure with respect to time (∂<sup>2</sup>P/∂t<sup>2</sup>) is approximated to be zero.

Click the “Choose folder” button to select the folder containing experimental data. The graphs will then self-populate. Note: Only .xls and .xlsx files are processed by the program; all other files are ignored. HyPAT assumes data is organized such that the data of each instrument (or set of instruments) has its own column in the excel sheet without a header.

The bottom left graph displays the permeability calculated from each data set plotted against temperature. Use the “Current Measurement” drop-down menu to change to the diffusivity, solubility, or flux calculated from each file plotted against temperature (or pressure, in the case of flux).

Also visible are a pressure versus time plot, a permeability versus time plot, and a comparison of different diffusivity optimizations plot, all using data from the first file in the selected folder. Use the “Current File” drop-down menu to use data from a different file for these plots.

By clicking the “Settings” button, you gain the ability to edit some of how the program processes the data. Specifically, the following parameters are editable:

* The column of the excel sheet from which data is read for each instrument,
* The row at which the program starts reading data from the excel sheet,
* The number of thermocouples used to measure the gas temperature,
* Conversion factors for each instrument to allow the user to update for appropriate unit conversions,
* Constant and proportional errors of each instrument.

### Usage: Comparison with Literature Values with "Overview Plots" Tab

The “Overview Plots” tab allows quick comparison to literature values of diffusivity, solubility, and permeability taken from M. Shimada (2020). The user can select and remove materials displayed in these overview plots using the materials column. Materials can be added to, modified, or removed from the database with the “Add/Edit Material” button.

Reference:

M. Shimada, “Tritium Transport in Fusion Reactor Materials” in Comprehensive Nuclear Materials 6, Second Edition, Elsevier (2020); https://doi.org/10.1016/B978-0-12-803581-8.11754-0

### Usage: Predict Outcomes with "Permeation Estimates” Tab

This tab provides predictions for permeation experiments based on specified user inputs. The required user inputs are highlighted in yellow. The sample properties are thickness (mm) and material. The experimental apparatus parameters are O-ring/sealing method, calibrated leak rate (mol s<sup>-1</sup> Torr<sup>-1</sup>), secondary side volume (cc), and accumulation time (hr). The specific test parameters are sample temperature (°C) and primary side pressure (Torr).

The results reported in the "FINAL OUTPUT” section include estimated time lag, permeation rate, and final secondary pressure.

By clicking the “Settings” button, O-rings can be added, edited, or removed. The default values for calibrated leak rate and secondary side volume can also be edited there. The options for sample material can be edited using the “Add/Edit Material” button on the Overview Plots tab.

## Installation Instructions

1.	Select the button “Code” near the top of this directory
2.	Select “Download ZIP” from the dropdown menu that appears
3.	Extract the downloaded file to a folder of your choice
4.	Open the command line, then move to the directory “HyPAT” within the downloaded folder
5.	Run the following code to open the application: python main.py

Required Python Packages: matplotlib, pandas, numpy, tkmacosx, mplcursors, scipy, and openpyxl

Compatible: macOS and Windows 10

## Plans

 A new tab is currently under development to analyze hydrogen absorption data in a Sieverts’ type apparatus. 
 
 ## Contribution Guidelines
 
 If you notice a feature that you want or a bug that needs to be fixed, please contact Chase Taylor or Thomas Fuerst.
 
 ## Authors and Acknowledgement 

Joseph Watkins

George Evans

Thomas Fuerst

Chase Taylor 

Joey Wakins (https://github.com/joeymwatkins) built most of the original functionality. George Evans (https://github.com/GeorgeEvans0) cleaned up errors, made various visual fixes, and added functionality. Thomas Fuerst (https://github.com/FuerstT) was the primary supervisor over this project, with Chase Taylor (https://github.com/taylchas) as secondary supervisor.

Work supported through the INL Laboratory Directed Research & Development (LDRD) Program under DOE Idaho Operations Office Contract DE-AC07-05ID14517

## License

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the file “LICENSE” for the specific language governing permissions and limitations under the License.
