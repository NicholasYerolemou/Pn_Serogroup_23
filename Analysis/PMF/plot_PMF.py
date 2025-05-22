
#Author: N.Yerolemou
#Plot Potential Mean free energy
#Usage
#Create a PMF dataclass for each glycosidic linkage you have in main()
#Plot the figure using plot_contourmap

from dataclasses import dataclass, field
import scipy.ndimage
from matplotlib import cm
import matplotlib.pyplot as plt 

#Dataclass for PMF data, handles all data manipulation and processing
@dataclass
class PMF:
    LinkageName: str = ""
    PATH: str = ""
    Filename: str = ""
    PHI: list = field(default_factory=list)
    PSI: list = field(default_factory=list)
    Energy: list = field(default_factory=list)
    colours: cm = field(default_factory=lambda: cm.coolwarm)

    def __post_init__(self):
        self.read_PMF_data()

    def read_PMF_data(self):
        # Read in data
        with open(self.PATH + self.Filename) as f:
            lines = f.readlines()
            xline = []
            yline = []
            zline = []
            for line in lines:
                if line[0] != '#' and len(line) >= 3:
                    nextX = float(line.split()[0])
                    if xline and (xline[-1] != nextX):  # end of row
                        self.PHI.append(xline)
                        self.PSI.append(yline)
                        self.Energy.append(zline)
                        xline = []
                        yline = []
                        zline = []
                    xline.append(nextX)
                    yline.append(float(line.split()[1]))
                    zline.append(float(line.split()[2]))
            self.PHI.append(xline)
            self.PSI.append(yline)
            self.Energy.append(zline)



def Main():
    Gro1P_O3_bDGal = PMF(LinkageName="Gro1PO3bDGal", PATH ="/home/nicholas-yerolemou/Documents/UCT/PhD/Simulation/Pn18/PMF/Gro1PO3bDGal/",Filename="Gro1PO3bDGal.pmf")
    # Glc14Gal_old = PMF(LinkageName="bDGlc14bDGal", PATH ="/home/nicholas-yerolemou/Documents/UCT/Masters/Simulation/PMF/Single Linkages/bDGlc14bDGal/",Filename="bDGlc14bDGal_PMF.pmf")


    plot_contourmap(Gro1P_O3_bDGal,"Gro1P_O3_bDGal")
    


class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

def plot_contourmap(pmf:PMF,title:str):

    fig, ax = plt.subplots(figsize=[6, 5], dpi=160)
    #smoothing data
    energy = scipy.ndimage.gaussian_filter(pmf.Energy, sigma=1) #bigger sigma = more smoothing; can go <1

    levs=[1,2,3,4,5,6,7,8,9,10] #levels to draw contours at
    CS = ax.contour(pmf.PHI, pmf.PSI, energy, levs, cmap=pmf.colours,zorder=1)

    
    #add the level labels
    CS.levels = [nf(val) for val in CS.levels]
    ax.clabel(CS, CS.levels, inline=True, fmt='%r ', fontsize=8)
    # CB = plt.colorbar(CS, shrink=0.8, extend='both') #If you want the colour scale on the plot


    plt.title(title)
    plt.xlabel(r'$\phi$', fontsize=12)
    plt.ylabel(r'$\psi$', fontsize=12)

    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_aspect('equal')

    plt.savefig(pmf.PATH+pmf.LinkageName)
    plt.show()


if __name__ == "__main__":
    Main()