import numpy as np
import matplotlib.pyplot as plt
import os
from wave_functions import *
import sys
from scipy import optimize

targetValuesH = np.zeros((3,3,3))
targetValuesH[1,0,0] = 73
targetValuesH[0,1,0] = 80.2
targetValuesH[0,0,1] = 101.1
targetValuesH[2,0,0] = 159.3
targetValuesH[1,1,0] = 160.2
targetValuesH[1,0,1] = 173.3
targetValuesH[0,2,0] = 181.9
targetValuesH[0,1,1] = 184.1
targetValuesH[0,0,2] = 219.4

targetValuesD = np.zeros((3,3,3))
targetValuesD[1,0,0] = 46.1
targetValuesD[0,1,0] = 46.5
targetValuesD[0,0,1] = 61.8
targetValuesD[2,0,0] = 99.0
targetValuesD[1,1,0] = 99.1
targetValuesD[1,0,1] = 106.5
targetValuesD[0,2,0] = 114.2
targetValuesD[0,1,1] = 114.5
targetValuesD[0,0,2] = 137.1

masses = {'H2':h2Mass, 'D':deuteriumMass}
targetValues = {'H2':targetValuesH, 'D':targetValuesD}

def EnergyDifferenceBox(variables, particle):
    mass = masses[particle]
    target = targetValues[particle]

    transitionList = [[Index, Value] for Index, Value in np.ndenumerate(target) if Value != 0]
    difference = 0
    for [transition,value] in transitionList:
        transitionBox = particleInABox3DEnergy(variables, [transition[0]+1, transition[1]+1,transition[2]+1], mass) - particleInABox3DEnergy(variables, [1,1,1],mass)
        difference += np.abs(transitionBox - value)

    return difference

def EnergyDifferenceSHM(variables, particle):
    mass = masses[particle]
    target = targetValues[particle]

    transitionList = [[Index, Value] for Index, Value in np.ndenumerate(target) if Value != 0]
    difference = 0
    for [transition,value] in transitionList:
        transitionBox = simpleHarmonicOscillator3DEnergy(variables, transition, mass) - simpleHarmonicOscillator3DEnergy(variables, [0,0,0],mass)
        difference += np.abs(transitionBox - value)

    return difference

if __name__ == '__main__':
    minimisedBoxH2 = optimize.minimize(EnergyDifferenceBox, [1,1,1], args=('H2'))
    minimisedBoxD = optimize.minimize(EnergyDifferenceBox, [1,1,1], args=('D'))
    minimisedSHMH2 = optimize.minimize(EnergyDifferenceSHM, [1,1,1],args=('H2'))
    minimisedSHMD = optimize.minimize(EnergyDifferenceSHM, [1,1,1],args=('D'))

    print(f'Box H2: {minimisedBoxH2.x}\n Box D: {minimisedBoxD.x} \n SHM H2: {minimisedSHMH2.x} \n SHM D: {minimisedSHMD.x}')
