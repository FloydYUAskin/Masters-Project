from wave_functions import *
from function_fitting_final import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

class EnergyDiagram():
    def __init__(self, desiredMass,accuracy):
        self.Mass = masses[desiredMass]
        self.desiredMass = desiredMass
        self.target = targetValues[desiredMass]
        self.boxVariables, self.oscVariables = FitFunction(0.3,0.6,1,4, accuracy, masses[desiredMass],targetValues[desiredMass])
        self.accuracy = accuracy

    def createEnergyMesh(self):
        nRangeBox = np.arange(1,4,1)
        nRangeOsc = np.arange(0,3,1)
        nMeshBox = np.meshgrid(nRangeBox,nRangeBox,nRangeBox)
        nMeshOsc = np.meshgrid(nRangeOsc,nRangeOsc,nRangeOsc)

        self.particleInABox = particleInABox3DEnergy(self.boxVariables, nMeshBox, self.Mass) - particleInABox3DEnergy(self.boxVariables, [1,1,1], self.Mass)
        self.simpleHarmonicOscillator = simpleHarmonicOscillator3DEnergy(self.oscVariables, nMeshOsc, self.Mass) - simpleHarmonicOscillator3DEnergy(self.oscVariables, [0,0,0], self.Mass)

        self.particleInABox[self.target==0] = 0
        self.simpleHarmonicOscillator[self.target==0] = 0

    def plotEnergyLevels(self):
        matplotlib.style.use('ggplot')
        import matplotlib.patches as mpatches
        self.fig, ((self.ax1, self.ax2), (self.ax3,self.ax4)) = plt.subplots(2,2)

        self.fig.suptitle(f'Energy Diagrams for {self.desiredMass}')
        simulatedEnergy = mpatches.Patch(color='red', label='Simulated Energy')
        actualEnergy = mpatches.Patch(color='blue', label='Actual Energy')
        plt.legend(handles=[simulatedEnergy, actualEnergy],loc = (0.,-.25))

        self.ax1.hlines(self.particleInABox[:,0,0],0,1)
        self.ax1.hlines(self.target[:,0,0],0,1,color='b')
        self.ax1.set_title('Particle In a Box nx Excitations',fontsize=10)
        self.ax1.set_ylabel('Energy (cm^-1)')
        self.ax1.set_xticks([])

        self.ax2.hlines(self.simpleHarmonicOscillator[:,0,0],0,1)
        self.ax2.hlines(self.target[:,0,0],0,1,color='b')
        self.ax2.set_title('Simple Harmonic Oscillator nx Excitations',fontsize=10)
        self.ax2.set_ylabel('Energy (cm^-1)')
        self.ax2.set_xticks([])

        self.ax3.hlines(self.particleInABox[:,:,:],0,1)
        self.ax3.hlines(self.target[:,:,:],0,1,color='b')
        self.ax3.set_title('Particle In a Box all Excitations',fontsize=10)
        self.ax3.set_ylabel('Energy (cm^-1)')
        self.ax3.set_xticks([])

        self.ax4.hlines(self.simpleHarmonicOscillator[:,:,:],0,1)
        self.ax4.hlines(self.target[:,:,:],0,1,color='b')
        self.ax4.set_title('Simple Harmonic Oscillator all Excitations',fontsize=10)
        self.ax4.set_ylabel('Energy (cm^-1)')
        self.ax4.set_xticks([])

        plt.savefig(f'Graphs-Diagrams/{self.desiredMass}-EnergyDiagram')
        

if __name__ == '__main__':
    Deuterium = EnergyDiagram('D', int(sys.argv[1]))
    Deuterium.createEnergyMesh()
    Deuterium.plotEnergyLevels()
    Hydrogen = EnergyDiagram('H2', int(sys.argv[1]))
    Hydrogen.createEnergyMesh()
    Hydrogen.plotEnergyLevels()
