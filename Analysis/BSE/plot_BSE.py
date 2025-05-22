# Author: N. Yerolemou - Sep 2024
# Calulate and plot Block Standard Averaging based on E2E or RGYR data

# Usage
# 1. Extract data for each molecules E2E distance or Radius of Gyration or both before running this script
# 2. Create BSE objects for each molecule in main(), setting the file paths and names and other molecule specific values, such as sim length as needed
# 3. Plot the data in main(). Plot the e2e BSE data with plot_e2e_BSE, rgyr BSE data with plot_rgyr_BSE or both with plot_both_BSE
# 4. Run the script with python3 plot_BSE.py

from dataclasses import dataclass, field
import os
import matplotlib.pyplot as plt
import numpy as np

@dataclass
class BSE:
    #This class stores all info related to a single molecule's BSE data
    Name: str

    E2E_PATH: str #Path to e2e distance data
    E2E_FILENAME: str #Name of e2e distance data file

    RGYR_PATH: str #Path to rgyr data
    RGYR_FILENAME: str #Name of rgyr data file

    BSE_output_PATH: str #Path to where you want BSE output file saved

    simLength: int = 1000  # Length of simulation in ns
    maxBlockSize: int = 0.1 * simLength  # Maximum block size in ns, recommended 5-10% of sim length
    samplingFactor: int = 1  # Frame stride - speeds up calculation

    force_recalculate: bool = False  # Force the script to recalculate the BSE file rather than reading in existing file

    # Variables for storing BSE data, Nind, and Tcorr for both E2E and RGYR
    Nind_E2E: float = 0.0
    Tcorr_E2E: float = 0.0
    BSE_data_E2E: list = field(default_factory=list)

    Nind_RGYR: float = 0.0
    Tcorr_RGYR: float = 0.0
    BSE_data_RGYR: list = field(default_factory=list)

    #Plotting variables
    colour:str = 'r'
    linestyle:str = "-"

    #Default constructor, checks if BSE data needs to be calculated or can be read in
    def __post_init__(self):
        self.bse_e2e_file = os.path.join(self.BSE_output_PATH, f"{self.Name}_E2E_BSE.txt")
        self.bse_rgyr_file = os.path.join(self.BSE_output_PATH, f"{self.Name}_RGYR_BSE.txt")

        # Check if the BSE data has already been calculated or needs recalculating
        #Calculate E2E BSE data
        if self.force_recalculate or not os.path.exists(self.bse_e2e_file):
            print(f"Calculating BSE for {self.E2E_FILENAME}")
            e2e_values = self.read_time_series(self.E2E_PATH+self.E2E_FILENAME)
            if e2e_values: self.write_BSE(e2e_values, 'E2E')
        else:
            print("Using existing BSE E2E file...")
            self.readBSE(self.bse_e2e_file,"E2E")

        #Calculate RGYR BSE data
        if self.force_recalculate or not os.path.exists(self.bse_rgyr_file):
            print(f"Calculating BSE for {self.RGYR_FILENAME}")
            rgyr_values = self.read_time_series(self.RGYR_PATH+self.RGYR_FILENAME)
            if rgyr_values: self.write_BSE(rgyr_values, 'RGYR')
        else:
            print("Using existing BSE RGYR file...")
            self.readBSE(self.bse_rgyr_file,"RGYR")
        
    def read_time_series(self, file) -> list:
        """Reads in time series data from the given file."""
        values = []
        try:
            with open(file, "r") as f:
                for line in f:
                    frame, value = line.split()
                    values.append(float(value))
        except FileNotFoundError:
            print(f"Error: Could not open file '{file}'. File not found.")
            return None
        except Exception as e:
            print(f"Error: Could not read file '{file}'. Reason: {e}")
            return None
        
        return values

    def readBSE(self, filepath: str, dataType: str) -> None:
        print("Reading in data for "+self.Name)
        """Reads BSE data from the file and stores it in the appropriate variable for E2E or RGYR."""
        with open(filepath, 'r') as file:
            lines = file.readlines()

            # First line contains Nind and Tcorr
            header = lines[0].strip().split(',')
            for item in header:
                if "Nind" in item:
                    Nind_value = float(item.split('=')[1])
                if "Tcorr" in item:
                    Tcorr_value = float(item.split('=')[1])

            # Store the values in the appropriate variables
            if dataType == 'E2E':
                self.Nind_E2E = Nind_value
                self.Tcorr_E2E = Tcorr_value
                self.BSE_data_E2E = [tuple(map(float, line.split())) for line in lines[2:]]  # Skip second line
            elif dataType == 'RGYR':
                self.Nind_RGYR = Nind_value
                self.Tcorr_RGYR = Tcorr_value
                self.BSE_data_RGYR = [tuple(map(float, line.split())) for line in lines[2:]]

    def write_BSE(self, values: list, dataType: str) -> None:
        """Calculates and writes BSE data, with correlation coefficients as the header."""
        BSEvalues = []
        X = []
        timeFactor = self.simLength / len(values) * 1000  # Picoseconds per frame
        timeFactor = round(timeFactor, 1) / 1000

        for blockSize in range(1, int(self.maxBlockSize / timeFactor), self.samplingFactor):
            count = 0
            valueSum = 0
            averageArr = []
            for v in values:
                if count >= blockSize:
                    averageArr.append(valueSum / count)
                    valueSum = 0
                    count = 0
                valueSum += v
                count += 1

            stdDev = np.std(averageArr)
            numBlocks = len(values) // blockSize
            X.append(round(blockSize * timeFactor, 2))
            BSEvalues.append(stdDev / np.sqrt(numBlocks))

        # Calculate final BSE and correlation values
        final_BSE = np.mean(BSEvalues[-10:])
        N_independent = pow(np.std(values) / final_BSE, 2)
        correlation_time = self.simLength / N_independent

        # Write out the data, with correlation coefficients as header
        Output_File = self.BSE_output_PATH + self.Name + "_" + dataType + "_BSE.txt"
        with open(Output_File, "w") as data_output:
            # Write correlation values as header to the BSE file
            data_output.write(f"#Correlation Values: Nind={N_independent:.3f}, Tcorr={correlation_time:.3f}\n")
            data_output.write(f"#Simlength:{self.simLength}ns, MaxBlockSize:{self.maxBlockSize}ns\n")
            for x, bse in zip(X, BSEvalues):
                data_output.write(f"{x} {bse}\n")  # Write actual data

        # # Store values in the appropriate variables
        if dataType == 'E2E':
            self.Nind_E2E = N_independent
            self.Tcorr_E2E = correlation_time
            self.BSE_data_E2E = list(zip(X, BSEvalues))
        elif dataType == 'RGYR':
            self.Nind_RGYR = N_independent
            self.Tcorr_RGYR = correlation_time
            self.BSE_data_RGYR = list(zip(X, BSEvalues))

def main():
    # Initialize the BSE object with your file paths and settings
    Pn23bb = BSE(
        Name = "Pn23bb",
        E2E_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/e2e/",
        E2E_FILENAME="Pn23bb_6RU_0_to_1000ns_e2e.txt",
        RGYR_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/rgyr/",
        RGYR_FILENAME = "Pn23bb_6RU_0_to_1000ns_rgyr.txt",
        BSE_output_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/BSE/",
        force_recalculate=False, #Force recaclulation of BSE values
        colour = 'dimgrey'
    )

    Pn23bb_Rha = BSE(
        Name = "Pn23bb+Rha",
        E2E_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/e2e/",
        E2E_FILENAME="Pn23bb_Rha_6RU_0_to_1000ns_e2e.txt",
        RGYR_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/rgyr/",
        RGYR_FILENAME = "Pn23bb_Rha_6RU_0_to_1000ns_rgyr.txt",
        BSE_output_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/BSE/",
        force_recalculate=False, #Force recaclulation of BSE values
        colour = 'darkgreen'
    )

    Pn23B = BSE(
        Name = "Pn23B",
        E2E_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23B_6RU/Analysis/e2e/",
        E2E_FILENAME="Pn23B_6RU_0_to_1000ns_e2e.txt",
        RGYR_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23B_6RU/Analysis/rgyr/",
        RGYR_FILENAME = "Pn23B_6RU_0_to_1000ns_rgyr.txt",
        BSE_output_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23B_6RU/Analysis/BSE/",
        force_recalculate=False, #Force recaclulation of BSE values
        colour = 'darkorange'
    )

    Pn23F = BSE(
        Name = "Pn23F",
        E2E_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/e2e/",
        E2E_FILENAME="Pn23F_6RU_V2_0_to_1000ns_e2e.txt",
        RGYR_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/rgyr/",
        RGYR_FILENAME = "Pn23F_6RU_V2_0_to_1000ns_rgyr.txt",
        BSE_output_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/BSE/",
        force_recalculate=False, #Force recaclulation of BSE values
        colour = 'darkblue'

    )

    Pn23A = BSE(
        Name = "Pn23A",
        E2E_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/9RU/Pn23A_9RU/Analysis/e2e/",
        E2E_FILENAME="Pn23A_9RU_0_to_1000ns_e2e.txt",
        RGYR_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/9RU/Pn23A_9RU/Analysis/rgyr/",
        RGYR_FILENAME = "Pn23A_9RU_0_to_1000ns_rgyr.txt",
        BSE_output_PATH="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/9RU/Pn23A_9RU/Analysis/BSE/",
        force_recalculate=False, #Force recaclulation of BSE values
        colour = 'darkred'
    )

    mols = [Pn23bb,Pn23bb_Rha,Pn23B,Pn23F,Pn23A]
    plot_e2e_BSE(mols)
    plot_rgyr_BSE(mols)
    # plot_both_BSE(mols,mols)

def plot_e2e_BSE(BSEs: list[BSE]) -> None:
    """Plots the BSE data for the E2E distances for a list of BSE objects."""
    plt.figure(figsize=(10, 6),dpi=160)

    for bse in BSEs:
        if bse.BSE_data_E2E:
            # Extract block sizes and BSE values for E2E
            block_sizes, bse_values = zip(*bse.BSE_data_E2E)
            plt.plot(block_sizes, bse_values, label=f"{bse.Name}: Nind={bse.Nind_E2E:.2f}, Tcorr={bse.Tcorr_E2E:.2f}ns", color=bse.colour, linestyle=bse.linestyle)

    # Add titles and labels
    plt.title("BSE on end-to-end distance",fontsize=38)
    plt.xlabel("Block Size (ns)",fontsize=30)
    plt.ylabel("Block Standard Error",fontsize=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Add legend
    plt.legend(fontsize=20,loc='upper left')
    
    # Display plot
    plt.grid(True)
    plt.tight_layout()
    # plt.savefig("")
    plt.show()

def plot_rgyr_BSE(BSEs: list[BSE]) -> None:
    """Plots the BSE data for the radius of gyration for a list of BSE objects."""
    plt.figure(figsize=(10, 6),dpi=160)

    for bse in BSEs:
        if bse.BSE_data_RGYR: #Check the data is not None
            # Extract block sizes and BSE values for RGYR
            block_sizes, bse_values = zip(*bse.BSE_data_RGYR) #Check values are not none
            plt.plot(block_sizes, bse_values, label=f"{bse.Name}: Nind={bse.Nind_RGYR:.2f}, Tcorr={bse.Tcorr_RGYR:.2f}ns", color=bse.colour, linestyle=bse.linestyle)

    # Add titles and labels
    plt.title("BSE on radius of gyration",fontsize=38)
    plt.xlabel("Block Size (ns)",fontsize=30)
    plt.ylabel("Block Standard Error",fontsize=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    # Add legend
    plt.legend(fontsize=20, loc='upper left')

    # Display plot
    plt.grid(True)
    plt.tight_layout()
    # plt.savefig("")
    plt.show()

def plot_both_BSE(E2E_BSEs: list[BSE], RGYR_BSEs: list[BSE]) -> None:
    """Plots the BSE data for the E2E distances and RGYR for a list of BSE objects."""
    plt.figure(figsize=(10, 6))

    for bse in E2E_BSEs: 
        if bse.BSE_data_E2E: #Check the data is not None
            # Extract block sizes and BSE values for E2E
            block_sizes, bse_values = zip(*bse.BSE_data_E2E)
            plt.plot(block_sizes, bse_values, label=f"{bse.Name}: Nind={bse.Nind_E2E:.2f}, Tcorr={bse.Tcorr_E2E:.2f}ns", color=bse.colour)

    for bse in RGYR_BSEs:
        if bse.BSE_data_RGYR: #Check the data is not None
            # Extract block sizes and BSE values for RGYR
            block_sizes, bse_values = zip(*bse.BSE_data_RGYR) #Check values are not none
            plt.plot(block_sizes, bse_values, label=f"{bse.Name}: Nind={bse.Nind_RGYR:.2f}, Tcorr={bse.Tcorr_RGYR:.2f} ns", color=bse.colour)

    plt.title("")
    plt.xlabel("Block Size (ns)",fontsize=26)
    plt.ylabel("BSE",fontsize=26)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    # Add legend
    plt.legend(fontsize=16)

    # Display plot
    plt.grid(True)
    plt.tight_layout()
    # plt.savefig("")
    plt.show()

if __name__ == "__main__":
    main()
