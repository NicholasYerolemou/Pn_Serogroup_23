
#Author: N.Yerolemou
#Plot Potential Mean free energy
#Usage
#Create a PMF dataclass for each glycosidic linkage you have in main()
#Read in the pmf data using read_PMF_data
#Plot the figure using plot_contourmap

from dataclasses import dataclass,field
import scipy.ndimage
from matplotlib import cm
import matplotlib.pyplot as plt 
import numpy as np
np.seterr(divide='ignore', invalid='ignore')


#Dataclass for PMF data, can be PHI,PSI,OMEGA,EPSILON on X or Y axis
@dataclass
class PMF:
    LinkageName: str = ""
    PATH: str = ""
    Filename:str = ""
    
    X: list[float] = None
    Y: list[float] = None
    Energy: list[float] = None

    x_axis_label: str = "phi"
    y_axis_label: str = "psi"
    
    colours: cm = field(default_factory=lambda: cm.gray) #Contour map colours

    def __post_init__(self):#Default constructor loads data
        self.X,self.Y,self.Energy = self.read_PMF_data(self.PATH + self.Filename)

    def read_PMF_data(self,File):#reads in data from given file
        x,y,energy = [],[],[]
        with open(File) as f:
            lines = f.readlines()
            xline=[]
            yline=[]
            zline=[]
            for line in lines:
                if line[0]!='#' and len(line)>=3:
                    nextX=float(line.split()[0])
                    if xline and (xline[-1]!=nextX): #end of row
                        x.append(xline)
                        y.append(yline)
                        energy.append(zline)
                        xline=[]
                        yline=[]
                        zline=[]
                    xline.append(nextX)
                    yline.append(float(line.split()[1]))
                    zline.append(float(line.split()[2]))
            x.append(xline) #do last append
            y.append(yline)
            energy.append(zline)
        return (x,y,energy)


@dataclass
class singleLinkage:
    # Data for a single glycosidic linkage
    PHI_atoms: list[int] = None  # atoms that form PHI
    PSI_atoms: list[int] = None
    EPSILON_atoms: list[int] = None
    OMEGA_atoms: list[int] = None

    occurance: str = None  # A letter indicating which linkage in the molecule this is

    # Data
    PHI: list[float] = None
    PSI: list[float] = None
    EPSILON: list[float] = None
    OMEGA: list[float] = None


@dataclass
class Dihedrals:
    LinkageName: str = ""
    PATH: str = ""
    Filename: str = ""

    Frames: list[int] = None
    per_linkage_data: list[singleLinkage] = None
    
    colours: cm = field(default_factory=lambda: cm.coolwarm)  # Contour map colours

    def __post_init__(self):
        # Initialize an empty list to hold singleLinkage data
        self.per_linkage_data = []
        self.read_Dihedral_data(self.PATH + self.Filename)

    def read_Dihedral_data(self, File: str):
        frames = []
        current_linkage = None

        with open(File, "r") as file:
            for line in file:
                if not line.strip():
                    continue

                # Process linkage headers (e.g., # Linkage Occurance B)
                if line.startswith("# Linkage Occurance") or not(current_linkage):
                    occurance = line.split()[-1].strip()
                    current_linkage = singleLinkage(
                        occurance=occurance, PHI_atoms=[], PSI_atoms=[], EPSILON_atoms=[], OMEGA_atoms=[], PHI=[], PSI=[], EPSILON=[], OMEGA=[])
                    self.per_linkage_data.append(current_linkage)

                # Process PHI atom list
                elif line.startswith("#PHI Atoms:"):
                    atoms = list(map(int, line.split(":")[1].strip().split()))
                    current_linkage.PHI_atoms = atoms

                # Process PSI atom list
                elif line.startswith("#PSI Atoms:"):
                    atoms = list(map(int, line.split(":")[1].strip().split()))
                    current_linkage.PSI_atoms = atoms

                # Process EPSILON atom list
                elif line.startswith("#EPSILON Atoms:"):
                    atoms = list(map(int, line.split(":")[1].strip().split()))
                    current_linkage.EPSILON_atoms = atoms

                # Process OMEGA atom list
                elif line.startswith("#OMEGA Atoms:"):
                    atoms = list(map(int, line.split(":")[1].strip().split()))
                    current_linkage.OMEGA_atoms = atoms

                # Process dihedral data
                elif not line.startswith("#"):
                    parts = line.strip().split(",")
                    if len(parts) >= 3:
                        try:
                            frame = int(parts[0])
                            phi_val = float(parts[1])
                            psi_val = float(parts[2])

                            epsilon_val = float(parts[3]) if len(parts) > 3 and parts[3].strip() else None
                            omega_val = float(parts[4]) if len(parts) > 4 and parts[4].strip() else None

                            frames.append(frame)
                            current_linkage.PHI.append(phi_val)
                            current_linkage.PSI.append(psi_val)
                            current_linkage.EPSILON.append(epsilon_val)
                            current_linkage.OMEGA.append(omega_val)

                        except ValueError as e:
                            print(f"Skipping line due to conversion error: {line.strip()}. Error: {e}")
                    else:
                        print(f"Skipping line with insufficient data: {line.strip()}")

        self.Frames = frames



def Main():

    #Create PMF and Dihedral objects for each linkage type

    #Below is a PMF object that holds x and y values with an associated energy
    Gal_14_Rha_PMF =PMF(LinkageName="Gal_14_Rha", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/PMF/Single Linkages/bDGal14bLRha/",Filename="bDGal14bLRha_PMF.pmf")
    #Below is Dihedral object that holds PHI,PSI values and OMEGA,EPSILON values if they exist.
    Gal_14_Rha_Dihed = Dihedrals(LinkageName="Gal_14_Rha", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/Dihedrals/200_to_1000ns/",Filename="Pn23bb_bDGal_14_bLRha_Dihedrals.txt")


    Glc_14_Gal_PMF = PMF(LinkageName="Glc_14_Gal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/PMF/Single Linkages/bDGlc_14_bDGal/",Filename="bDGlc_14_bDGal_PMF.pmf")
    Glc_14_Gal_Dihed = Dihedrals(LinkageName="Glc_14_Gal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/Dihedrals/200_to_1000ns/",Filename="Pn23bb_bDGlc_14_bDGal_Dihedrals.txt")

    Rha_14_Glc_PMF = PMF(LinkageName="Rha_14_Glc", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/PMF/Single Linkages/bLRha14bDGlc/",Filename="bLRha14bDGlc_PMF.pmf")
    Rha_14_Glc_Dihed = Dihedrals(LinkageName="Rha_14_Glc", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_6RU/Analysis/Dihedrals/200_to_1000ns/",Filename="Pn23bb_bLRha_14_bDGlc_Dihedrals.txt")

    Rha_12_Gal_PMF = PMF(LinkageName="aLRha_12_bDGal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/PMF/Single Linkages/aLRha12bDGal/",Filename="aLRha12bDGal_PMF.pmf")
    Rha_12_Gal_Dihed = Dihedrals(LinkageName="aLRha12bDGal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23bb_Rha_6RU/Analysis/Dihedrals/200_to_1000ns/",Filename="Pn23bb+Rha_aLRha_12_bDGal_Dihedrals.txt")


    Gro_2P3_Gal_Dihed = Dihedrals(LinkageName="Gro_2P3_Gal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/Pn23/6RU/Pn23F_6RU_V2/Analysis/Dihedrals/200_to_1000ns/",Filename="Pn23F_6RU_G2P_3_Gal_Dihedrals.txt")
    
    
    #Plot a combination of PMF and Dihderal data on the same plot
    plot_both(pmf=Gal_14_Rha_PMF,dihedrals=Gal_14_Rha_Dihed,title="bDGal(1->4)bLRha",x_axis="phi",y_axis="psi")#This plots PHI vs PSI of the Gal_14_Rha linkage
    plot_both(pmf=Glc_14_Gal_PMF,dihedrals=Glc_14_Gal_Dihed,title="bDGlc(1->4)bDGal",x_axis="phi",y_axis="psi")
    plot_both(pmf=Rha_14_Glc_PMF,dihedrals=Rha_14_Glc_Dihed,title="bLRha(1->4)bDGlc",x_axis="phi",y_axis="psi")

    
    # plot_both(pmf=Rha_12_Gal_PMF,dihedrals=Rha_12_Gal_Dihed,title="aLRha(1->2)bDGal",x_axis="phi",y_axis="psi")
    # plot_both(pmf=None,dihedrals=Gro_2P3_Gal_Dihed,title="Gro_2P3_Gal",x_axis="phi",y_axis="psi")#This plots just the dihedrals PHI vs PSI of the Gro_2P3_Gal linkage
    # plot_both(pmf=None,dihedrals=Gro_2P3_Gal_Dihed,title="Gro_2P3_Gal",x_axis="phi",y_axis="Omega") #This plots just the dihedrals PHI vs OMEGA of the Gro_2P3_Gal linkage
    # plot_both(pmf=None,dihedrals=Gro_2P3_Gal_Dihed,title="Gro_2P3_Gal",x_axis="phi",y_axis="Epsilon")
    


class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

#Plot PMF contour map
def plot_contourmap(pmf:PMF,title:str="",ax:plt.Axes = None):

    save = False
    if ax is None:
        fig, ax = plt.subplots(figsize=[6, 5], dpi=160)
        save = True


    #smoothing data
    energy = scipy.ndimage.gaussian_filter(pmf.Energy, sigma=1) #bigger sigma = more smoothing; can go <1

    levs=[1,2,3,4,5,6,7,8,9,10] #levels to draw contours at
    CS = ax.contour(pmf.X, pmf.Y, energy, levs, cmap=pmf.colours,zorder=1)

    
    #add the level labels
    CS.levels = [nf(val) for val in CS.levels]
    ax.clabel(CS, CS.levels, inline=True, fmt='%r ', fontsize=8,zorder=1)
    # CB = plt.colorbar(CS, shrink=0.8, extend='both') #If you want the colour scale on the plot


    plt.title(title)
    plt.xlabel(rf'$\{pmf.x_axis_label.lower()}$', fontsize=12)
    plt.ylabel(rf'$\{pmf.y_axis_label.lower()}$', fontsize=12)


    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_aspect('equal')

    if(save):
        plt.show()

#Plot dihedral data
def plot_Dihedral_data(dihed: Dihedrals, title: str, ax: plt.Axes = None,x_axis: str="phi",y_axis: str="psi"):
    save = False
    if ax is None:
        fig, ax = plt.subplots(figsize=[6, 5], dpi=160)
        save = True

    for linkage in dihed.per_linkage_data:
        x_data = getattr(linkage, x_axis.upper(), None)#get x-axis data
        y_data = getattr(linkage, y_axis.upper(), None)#get y-axis data

        # Only plot if both x_data and y_data exist
        if x_data is not None and y_data is not None:
            ax.hist2d(
                x_data, y_data, 
                bins=(180, 180), 
                density=True, 
                cmin=0.00001, 
                range=[[-180, 180], [-180, 180]], 
                cmap=dihed.colours, 
                norm=cm.colors.Normalize(vmin=0, vmax=0.0007, clip=False),
                zorder=2
            )

    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_aspect('equal')


    plt.title(title)
    plt.xlabel(rf'$\{x_axis.lower()}$', fontsize=12)  # plot x axis symbol based on x_axis_label
    plt.ylabel(rf'$\{y_axis.lower()}$', fontsize=12)  # plot y axis symbol based on y_axis_label

    if(save):
        plt.show()

#Plot PMF and Dihedral data on same figure
def plot_both(pmf: PMF=None, dihedrals: Dihedrals=None, title: str = "PMF and Dihedral Data",x_axis: str="phi",y_axis: str="psi"):
    fig, ax = plt.subplots(figsize=[6, 5], dpi=160)
    
    # Plot PMF data (contour plot)
    if pmf:
        plot_contourmap(pmf, title=title, ax=ax)
    
    # Overlay Dihedral data (2D histogram)
    if dihedrals:
        plot_Dihedral_data(dihedrals, title=title, ax=ax,x_axis=x_axis,y_axis=y_axis)
    
    # Show the combined plot
    plt.savefig(dihedrals.PATH+title)
    plt.title(title)
    plt.show()

if __name__ == "__main__":
    Main()