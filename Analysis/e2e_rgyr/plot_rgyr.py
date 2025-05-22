# Author: Nicholas Yerolemou - Sep 2024

# Usage
# 1. Load your data into VMD
# 2. Update Simulation variables
# 3. Create Molecule objects in main() for each molecule's data
# 4. Call appropriate methods to plot the data and pass the Molecule object/objects
# 5. Run the script with "python3 plot_rgyr.py"

from dataclasses import dataclass
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec

# Simulation Variables
dcd_freq = 250
time_step = 1000000  # ns, 1fs = 1 million ns
extraction_stride = 100  # Stride used when extracting frames from original dcd file
vmd_stride = 1  # stride used when loading data into vmd
stride = extraction_stride * vmd_stride

# Define a dataclass for molecules
@dataclass
class Molecule:
    Name: str
    PATH: str
    FILENAME: str
    Data: tuple[list[int], list[int]] = ([], [])
    
    # default values
    line_x_limit: tuple[int, int] = (0, 1050)
    line_y_limit: tuple[int, int] = (6, 24)
    hist_x_limit: tuple[int, int] = (6, 22)
    hist_y_limit: tuple[int, int] = (0, 0.3)
    
    colour: str = 'k'
    annotation_pos: tuple[int, int] = (8, 0.2)
    
    # Constructor, reads data when object is initialized
    def __post_init__(self):
        self.Data = self.read_in_data(self.PATH + self.FILENAME)
    
    def read_in_data(self, File: str) -> tuple[list[int], list[int]]:
        X, Y = [], []
        for line in open(File, 'r'):
            values = [float(s) for s in line.split()]
            X.append(((dcd_freq / time_step) * stride * values[0]))
            Y.append(values[1])
        return X, Y


# Main function to load and plot data
def Main():
    # Define molecules
    Pn23bb = Molecule(Name="Pn23bb 6RU", PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/rgyr/", FILENAME="Pn23bb_6RU_0_to_1000ns_rgyr.txt", colour='darkblue')
    Pn23bb_Rha = Molecule(Name="Pn23bb+Rha 6RU", PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/rgyr/", FILENAME="Pn23bb_Rha_6RU_0_to_1000ns_rgyr.txt", colour='black')
    Pn23B = Molecule(Name="Pn23B 6RU", PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23B_6RU/Analysis/rgyr/", FILENAME="Pn23B_6RU_0_to_1000ns_rgyr.txt", colour='purple')
    Pn23F = Molecule(Name="Pn23F 6RU", PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/rgyr/", FILENAME="Pn23F_6RU_V2_0_to_1000ns_rgyr.txt", colour='green')
    Pn23A = Molecule(Name="Pn23A 9RU", PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/9RU/Pn23A_9RU/Analysis/rgyr/", FILENAME="Pn23A_9RU_0_to_1000ns_rgyr.txt", colour='darkred')
    
    
    Pn23A.hist_y_limit = (0,0.8)
    # Plot data
    # combined(Pn23bb, "Pn23bb 1000ns")
    # combined(Pn23bb_Rha, "Pn23bb_Rha 1000ns")

    # Plot multiple molecules
    mols = [Pn23bb, Pn23bb_Rha, Pn23B, Pn23F, Pn23A]
    plot_multiple_mols(mols, "All Molecules")

# Plotting functions
def line_graph(mol: Molecule, Title: str, ax: plt.Axes = None) -> None:
    """Plot line graph of single molecules data"""
    save = False
    if ax is None:
        save = True
        fig, ax = plt.subplots(figsize=[10, 6], dpi=160)
        ax.set_title(Title, fontsize=24)
    
    X, Y = mol.Data
    equlibration_run = int(200 * time_step / (dcd_freq * stride))
    
    ax.plot(X[:equlibration_run], Y[:equlibration_run], color=mol.colour, linewidth=0.5, alpha=0.55, label="Equilibration")
    ax.plot(X[equlibration_run:], Y[equlibration_run:], color=mol.colour, linewidth=0.5, label="Production run")
    
    ax.set_xlim(mol.line_x_limit)
    ax.set_ylim(mol.line_y_limit)
    
    ax.set_xlabel("Time (ns)", fontsize=22)
    ax.set_ylabel('Length (\u212B)', fontsize=22)
    ax.grid(True, which='major', linestyle='--', linewidth=0.5)

    if save:
        plt.savefig(f"{mol.PATH}{mol.Name}_line.png")
        plt.show()

def histogram(mol: Molecule, Title: str, ax: plt.Axes = None) -> None:
    """Plot histogram of single molecules data"""
    save = False
    if ax is None:
        save = True
        fig, ax = plt.subplots(figsize=(9, 7), dpi=120)
        ax.set_title(Title, fontsize=24)
    
    X, Y = mol.Data
    equlibration_run = int(200 * time_step / (dcd_freq * stride))
    
    # ax.hist(Y[:equlibration_run], color=mol.colour, alpha=0.55, label='Equilibration', density=True, edgecolor='white', bins=35, histtype='bar')
    n, bins, patches = ax.hist(Y[equlibration_run:], color=mol.colour, label='Production run', density=True, linewidth=0.3, edgecolor='white', bins=35, histtype='bar')
    
    sd = round(np.std(Y), 2)
    ax.annotate(f"\u03C3 = {sd}",mol.annotation_pos, fontsize=15, fontweight='bold')
    

    offset = 0.11*mol.hist_y_limit[1]#vertical distance between sd and mode annotations

    # Plot mode 
    # Find the bin with the highest count
    # Calculate the mode as the midpoint of this bin
    max_bin_index =  np.argmax(n)
    mode = (bins[max_bin_index] + bins[max_bin_index + 1]) / 2
    ax.axvline(mode, color='k', linestyle='-.', linewidth=1, label=f'Mode: {mode:.2f}')
    ax.annotate(f"{int(round(mode,0))} \u212B", (mol.annotation_pos[0],mol.annotation_pos[1]-offset), fontsize=15, fontweight='bold')


    #Add labels
    ax.set_title(Title, fontsize=24)
    ax.set_xlabel('Length (\u212B)', fontsize=22)
    ax.set_ylabel('Probability', fontsize=22)
    ax.tick_params(axis='both', labelsize=18)
    ax.grid(True, linestyle="--", linewidth=0.5)
    ax.set_xlim(mol.hist_x_limit)
    ax.set_ylim(mol.hist_y_limit)

    if save:
        plt.savefig(f"{mol.PATH}{Title}_histogram.png")
        plt.show()

def combined(mol: Molecule, Title: str) -> None:
    """Plot one molecules line graph and histogram data on same axes"""
    fig = plt.figure(figsize=(15, 5), dpi=160)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 2])

    ax1 = plt.subplot(gs[0])
    line_graph(mol, Title, ax1)

    ax2 = plt.subplot(gs[1])
    histogram(mol, Title, ax2)


    plt.tight_layout()
    plt.savefig(f"{mol.PATH}{mol.Name}_combined.png")
    plt.show()

def plot_multiple_mols(mols: list[Molecule], Title: str) -> None:
    """Plot multiple molecules line graph and histogram data on same axes"""
    num_mols = len(mols)
    fig = plt.figure(figsize=(12, 4 * num_mols))
    gs = fig.add_gridspec(num_mols, 2, width_ratios=[5, 1], height_ratios=[1] * num_mols)

    for i, mol in enumerate(mols):
        ax1 = plt.subplot(gs[i, 0])
        ax2 = plt.subplot(gs[i, 1])

        line_graph(mol, "", ax1)
        ax1.tick_params(axis='both', labelsize=18)
        ax1.set_xlabel("Time (ns)", fontsize=20, fontweight='bold')
        ax1.set_ylabel('Length (\u212B)', fontsize=20, fontweight='bold')
        ax1.grid(visible=False)



        histogram(mol, "", ax2)
        ax2.set_xlabel('Length (\u212B)', fontsize=18, fontweight='bold')
        ax2.set_ylabel('Probability', fontsize=18, fontweight='bold')
        ax2.tick_params(axis='both', labelsize=16)
        ax2.grid(False)

        # ax1.annotate(chr(65 + i), xy=(-0.25, 0.42), xycoords='axes fraction', fontsize=32, fontweight='bold')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    Main()
