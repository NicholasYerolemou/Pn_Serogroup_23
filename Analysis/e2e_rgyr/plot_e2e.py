# Author: Nicholas Yerolemou - Sep 2024

# Usage
# 1. Load your data into VMD
# 2. Create Molecule objects in main() for each molecule's e2e data you want to plot
# 3. Call appropriate methods to plot the data and pass the Molecule object/objects
# 4. Run the script with "python3 plot_e2e.py"

from dataclasses import dataclass
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec


#Simulation Variables
dcd_freq = 250
time_step = 1000000  # ns, 1fs = 1 million ns
extraction_stride = 100  # Stride used when extracting frames from original dcd file
vmd_stride = 1  # stride used when loading data into vmd
stride = extraction_stride * vmd_stride


#Define a dataclass for molecules
@dataclass
class Molecule:
    Name: str
    PATH: str
    FILENAME: str
    Data: tuple[list[int], list[int]] = ([],[])

    #default values
    line_x_limit: tuple[int, int] = (0,1050)
    line_y_limit: tuple[int, int] = (0,80)
    hist_x_limit: tuple[int, int] = (5, 70)
    hist_y_limit: tuple[int, int] = (0, 0.08)

    colour: str = 'k'
    annotation_pos: tuple[int, int] = (25, 0.07)

    fontScale:float = 1.0 #scale all the fonts on a figure


    #Constructor, reads data in when object initialized
    def __post_init__(self):
        self.Data = self.read_in_data(self.PATH + self.FILENAME)

    #reads in data
    def read_in_data(self,File: str) -> tuple[list[int], list[int]]:
        # Read in the data
        X, Y = [], []
        for line in open(File, 'r'):
            values = [float(s) for s in line.split()]
            X.append(((dcd_freq / time_step) * stride * values[0]))
            Y.append(values[1])
        return X, Y


def Main():

    #Define data for each molecule you wish to plot
    Pn23bb = Molecule(Name = "Pn23bb",PATH = "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/e2e/",FILENAME="Pn23bb_6RU_0_to_1000ns_e2e.txt",colour='dimgrey',fontScale=0.0)
    Pn23bb_Rha = Molecule(Name = "Pn23bb_Rha",PATH = "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/e2e/",FILENAME="Pn23bb_Rha_6RU_0_to_1000ns_e2e.txt",colour='darkgreen',fontScale=0.0)
    Pn23B = Molecule(Name = "Pn23B",PATH = "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23B_6RU/Analysis/e2e/",FILENAME="Pn23B_6RU_0_to_1000ns_e2e.txt",colour='darkorange',fontScale=0.0)
    Pn23F = Molecule(Name = "Pn23F_V2",PATH = "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/e2e/",FILENAME="Pn23F_6RU_V2_0_to_1000ns_e2e.txt",colour='darkblue',fontScale=0.0)
    Pn23A = Molecule(Name = "Pn23A",PATH = "/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/9RU/Pn23A_9RU/Analysis/e2e/",FILENAME="Pn23A_9RU_0_to_1000ns_e2e.txt",colour = "darkred")


    #Molecule Specific Formatting
    Pn23bb.annotation_pos = (10,0.07)

    Pn23A.hist_y_limit = (0,0.2)
    Pn23A.annotation_pos = (40,0.170)

    #Plot data
    # line_graph(Pn23bb,"Pn23bb 6RU E2E 0 to 1000ns")
    # histogram(Pn23bb,"Pn23bb 6RU E2E 0 to 1000ns")
    # combined(Pn23bb,"Pn23bb 6RU E2E 0 to 1000ns")

    # Mols = [Pn23bb,Pn23bb_Rha,Pn23B,Pn23F,Pn23A]
    Mols = [Pn23bb,Pn23bb_Rha,Pn23B,Pn23F,Pn23A]
    plot_multiple_mols(Mols,"E2E All")
    

def line_graph(mol:Molecule, Title: str, ax: plt.Axes = None) -> None:

    save = False
    if ax is None:
        save = True
        fig, ax = plt.subplots(figsize=[10, 6], dpi=160)


    X,Y = mol.Data

    equlibration_run = int(200 * time_step / (dcd_freq * stride))  # The first 200ns

    # plot the first 200ns as lighter
    ax.plot(X[0:equlibration_run], Y[0:equlibration_run], color=mol.colour, linewidth=1, alpha=0.55, label="Equilibration")
    # plot the remaining 800ns
    ax.plot(X[equlibration_run:], Y[equlibration_run:], color=mol.colour, linewidth=1, label="Production run")

    # Plot mean line
    mean_y = np.mean(Y)
    ax.axhline(mean_y, color='k', linestyle='--', linewidth=2, label=f'Mean: {mean_y:.2f}')
    ax.annotate(str(int(round(mean_y,0)))+ "\u212B",(1010, mean_y+2),fontsize=24, fontweight='bold')

    # Change the x and y limits here
    ax.set_xlim(mol.line_x_limit)
    ax.set_ylim(mol.line_y_limit)
    ax.tick_params(axis='both', labelsize=18)
    ax.set_title(Title, fontsize=24)
    ax.set_xlabel("Time (ns)", fontsize=22)
    ax.set_ylabel('Length (\u212B)', fontsize=22)
    ax.grid(visible=True, which='major', axis='both', linestyle='--', linewidth=0.5)

    if save:
        # plt.savefig(mol.PATH+Title + ".png")
        plt.show()


def histogram(mol:Molecule, Title: str, ax: plt.Axes = None) -> None:
    save = False
    if ax is None:
        save = True
        fig, ax = plt.subplots(figsize=(9, 7), dpi=120)

    X,Y = mol.Data

    equlibration_run = int(200 * time_step / (dcd_freq * stride))  # The first 200ns

    # First 200ns
    # ax.hist(Y[0:equlibration_run], color=mol.colour, alpha=0.55, label='200ns Equilibration', density=True, edgecolor='white',linewidth=0.3, bins=35, histtype='bar')
    # Remaining data (200ns-end)
    n, bins, patches = ax.hist(Y[equlibration_run:], color=mol.colour, label='800ns production run', density=True, edgecolor='white', linewidth=0.3, bins=35, histtype='bar')

    # Calculate standard deviation (sd)
    sd = str(int(round(np.std(Y), 0)))
    ax.annotate("\u03C3 = " + sd, mol.annotation_pos, fontsize=26, fontweight='bold')


    offset = 0.11*mol.hist_y_limit[1]
    

    # Plot mode 
    # Find the bin with the highest count
    # Calculate the mode as the midpoint of this bin
    max_bin_index =  np.argmax(n)
    mode = (bins[max_bin_index] + bins[max_bin_index + 1]) / 2
    ax.axvline(mode, color='k', linestyle='-.', linewidth=2, label=f'Mode: {mode:.2f}')
    ax.annotate(f"{int(round(mode,0))} \u212B", (mol.annotation_pos[0],mol.annotation_pos[1]-offset), fontsize=26, fontweight='bold')


    # Add labels
    ax.set_title(Title, fontsize=24)
    ax.set_xlabel('Length (\u212B)', fontsize=22)
    ax.set_ylabel('Probability', fontsize=22)
    ax.tick_params(axis='both', labelsize=18)
    ax.grid(True, linestyle="--", linewidth=0.5)
    ax.set_xlim(mol.hist_x_limit)
    ax.set_ylim(mol.hist_y_limit)

    if save:
        plt.savefig(mol.PATH+Title+ " histogram.png")
        plt.show()

#linegraaph and Histogram combined
def combined(mol:Molecule, Title: str) -> None:

    fig = plt.figure(figsize=(15, 5), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 2])

    # Plot the line graph
    ax1 = plt.subplot(gs[0])
    line_graph(mol,Title, ax1)
    ax1.tick_params(axis='both', labelsize=4)
    ax1.set_xlabel("Time (ns)", fontsize=10, fontweight='bold')
    ax1.set_ylabel('r (\u212B)', fontsize=10, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.grid(visible=False)

    # Plot the histogram
    ax2 = plt.subplot(gs[1])
    histogram(mol, Title, ax2)
    ax2.set_xlabel('r (\u212B)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Probability', fontsize=10, fontweight='bold')
    ax2.tick_params(axis='both', labelsize=8)
    ax2.grid(False)



    plt.tight_layout()
    # plt.savefig(mol.PATH+Title + " combined.png")
    plt.show()

#line graph and histogram for each molecule given on one figure
def plot_multiple_mols(Mols: list[Molecule], Title: str) -> None:
    num_mols = len(Mols)

    fig = plt.figure(figsize=(10, 5 * num_mols))  # Scales size of figure based on number of sublots
    gs = fig.add_gridspec(num_mols, 2, width_ratios=[5, 1], height_ratios=[1] * num_mols)  # creates N rows with 2 columns each and sets their width/height
    # plt.title(Title)

    # loop through each molecule and plot its data
    for i, mol in enumerate(Mols):
        ax1 = plt.subplot(gs[i, 0])
        ax2 = plt.subplot(gs[i, 1])

        # Plot the line graph
        line_graph(mol, "", ax1)
        ax1.tick_params(axis='both', labelsize=28)
        ax1.set_xlabel("Time (ns)", fontsize=36*mol.fontScale, fontweight='bold')
        ax1.set_ylabel('r (\u212B)', fontsize=36, fontweight='bold')
        ax1.grid(visible=False)

        # Plot the histogram
        histogram(mol, "", ax2)
        ax2.set_xlabel('r (\u212B)', fontsize=33*mol.fontScale, fontweight='bold')
        ax2.set_ylabel('Probability', fontsize=33, fontweight='bold')
        ax2.tick_params(axis='both', labelsize=26)
        ax2.grid(False)


        # mol_labels = ["Pn23bb","Pn23bb+Rha","Pn23B","Pn23F", "Pn23A"]
        # RU_labels = ["6RU", "6RU","6RU","6RU","9RU"]
        # labels = ["A","B","C","D","E"]
        # Add label to the left of the line graph
        # ax1.annotate(labels[i], xy=(-0.08, 0.45), xycoords='axes fraction', fontsize=32, fontweight='bold', ha='center', va='center')
        # ax1.annotate(mol_labels[i], xy=(-0.13, 0.5), xycoords='axes fraction', fontsize=18, ha='center', va='center')
        # ax1.annotate(RU_labels[i], xy=(-0.13, 0.35), xycoords='axes fraction', fontsize=18, ha='center', va='center')


    plt.tight_layout()
    plt.subplots_adjust(
    left=0.035,  # Space between the left edge of the figure and the plot
    right=0.993, # Space between the right edge of the figure and the plot
    bottom=0.057, # Space between the bottom edge of the figure and the plot
    top=0.988,   # Space between the top edge of the figure and the plot
    wspace=0.112, # Width space between subplots (if multiple)
    hspace=0.21  # Height space between subplots (if multiple)
)
    # plt.savefig("/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/"+Title)
    plt.show()


if __name__ == "__main__":
    Main()