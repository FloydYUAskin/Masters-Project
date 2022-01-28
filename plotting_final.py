import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from wave_functions import *
from scipy import optimize
from function_fitting_final import targetValuesH, targetValuesD, targetValues, masses, EnergyDifferenceBox, EnergyDifferenceSHM
import os

methods = [
'Nelder-Mead',
'Powell',
'CG',
'BFGS',
'Newton-CG',
'L-BFGS-B',
'TNC',
'COBYLA',
'SLSQP',
'trust-constr',
'dogleg',
'trust-ncg',
'trust-exact',
'trust-krylov'
]

class EnergyDiagram:
    def __init__(self, particle, method):
        self.mass = masses[particle]
        self.particle = particle
        self.method = method
        self.targetValues = targetValues[particle]

        self.BoxValues = optimize.minimize(EnergyDifferenceBox, [1,1,1], args=(self.particle),method = self.method).x
        self.SHMValues = optimize.minimize(EnergyDifferenceSHM, [1,1,1], args=(self.particle),method = self.method).x

        self.createEnergyArray()
        self.createPlots()

    def createEnergyArray(self):
        self.BoxEnergyLevels = np.arange(1,3,1)
        self.SHMEnergyLevels = np.arange(0,2,1)

        self.BoxMesh = np.meshgrid(self.BoxEnergyLevels, self.BoxEnergyLevels, self.BoxEnergyLevels)
        self.SHMMesh = np.meshgrid(self.SHMEnergyLevels, self.SHMEnergyLevels, self.SHMEnergyLevels)

        self.BoxEnergies = particleInABox3DEnergy(self.BoxValues, self.BoxMesh, self.mass)
        self.SHMEnergies = simpleHarmonicOscillator3DEnergy(self.SHMValues, self.SHMMesh, self.mass)

        # self.BoxEnergies[self.targetValues == 0] = 0
        # self.SHMEnergies[self.targetValues == 0] = 0

    def createPlots(self):
        matplotlib.style.use('ggplot')
        self.fig, ((self.ax1,self.ax2),(self.ax3,self.ax4)) = plt.subplots(2,2)

        self.fig.suptitle(f'Energy Diagrams for {self.particle}')
        simulatedEnergy = mpatches.Patch(color='red', label='Simulated Energy')
        actualEnergy = mpatches.Patch(color='blue', label='Actual Energy')
        plt.legend(handles=[simulatedEnergy, actualEnergy],loc = (0.,-.25))

        self.ax1.hlines(self.BoxEnergies[:,0,0],0,1)
        self.ax1.hlines(self.targetValues[:,0,0],0,1,color='b',linestyles='dashed')
        self.ax1.set_title('Particle In a Box nx Excitations',fontsize=10)
        self.ax1.set_ylabel('Energy (cm^-1)')
        self.ax1.set_xticks([])
        self.ax1.set_ylim([0,150])

        self.ax2.hlines(self.SHMEnergies[:,0,0],0,1)
        self.ax2.hlines(self.targetValues[:,0,0],0,1,color='b',linestyles='dashed')
        self.ax2.set_title('Simple Harmonic Oscillator nx Excitations',fontsize=10)
        self.ax2.set_ylabel('Energy (cm^-1)')
        self.ax2.set_xticks([])
        self.ax2.set_ylim([0,150])

        self.ax3.hlines(self.BoxEnergies[:,:,:],0,1)
        self.ax3.hlines(self.targetValues[:,:,:],0,1,color='b',linestyles='dashed')
        self.ax3.set_title('Particle In a Box all Excitations',fontsize=10)
        self.ax3.set_ylabel('Energy (cm^-1)')
        self.ax3.set_xlabel(f'Box Dimensions: {self.BoxValues} \n SHM Constants: {self.SHMValues}',fontsize=5)
        self.ax3.set_xticks([])

        self.ax4.hlines(self.SHMEnergies[:,:,:],0,1)
        self.ax4.hlines(self.targetValues[:,:,:],0,1,color='b',linestyles='dashed')
        self.ax4.set_title('Simple Harmonic Oscillator all Excitations',fontsize=10)
        self.ax4.set_ylabel('Energy (cm^-1)')
        self.ax4.set_xticks([])

        plt.savefig(f'Graphs-Diagrams/{self.particle}-EnergyDiagram-{self.method}')


if __name__ == '__main__':
    for method in methods:
        try:
            EnergyDiagram('D',method)
            EnergyDiagram('H2',method)
            os.system('cls')
            print(method)
        except ValueError:
            pass
