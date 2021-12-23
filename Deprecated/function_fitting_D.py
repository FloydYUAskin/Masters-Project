import numpy as np
import matplotlib.pyplot as plt
import os
from wave_functions import *

boxVarMin, boxVarMax = 0.4, 0.6
oscVarMin, oscVarMax = 1, 3
accuracy = 50

targetValues = np.zeros((3,3,3))
targetValues[1,0,0] = 46.1
targetValues[0,1,0] = 46.5
targetValues[0,0,1] = 61.8
targetValues[2,0,0] = 99.0
targetValues[1,1,0] = 99.1
targetValues[1,0,1] = 106.5
targetValues[0,2,0] = 114.2
targetValues[0,1,1] = 114.5
targetValues[0,0,2] = 137.1

transitionList = [[Index, Value] for Index, Value in np.ndenumerate(targetValues) if Value != 0]

boxVariableArray = np.linspace(boxVarMin, boxVarMax, accuracy)
oscVariableArray = np.linspace(oscVarMin, oscVarMax, accuracy)

differenceListBox = []
variableListBox = []

for xIndex, xVariable in enumerate(boxVariableArray):
    for yIndex, yVariable in enumerate(boxVariableArray):
        for zIndex, zVariable in enumerate(boxVariableArray):
            currentStep = (xIndex * accuracy ** 2) + (yIndex * accuracy) + zIndex
            print(f'{currentStep * 100 / accuracy**3}%')
            difference = 0
            for [transition, value] in transitionList:
                transitionBox = particleInABox3DEnergy([xVariable, yVariable, zVariable], [transition[0]+1, transition[1]+1, transition[2]+1], deuteriumMass) - particleInABox3DEnergy([xVariable, yVariable, zVariable], [1,1,1], deuteriumMass)
                difference += np.abs(transitionBox - value)
            variableListBox.append([xVariable, yVariable, zVariable])
            differenceListBox.append(difference)
            os.system('clear')

differenceListBox = np.array(differenceListBox)
variableListBox = np.array(variableListBox)

differencesInd = np.argsort(differenceListBox)

orderedDifferenceBox = differenceListBox[differencesInd]
orderedVariablesBox = variableListBox[differencesInd]

differenceList = []
variableList = []

for xIndex, xVariable in enumerate(oscVariableArray):
    for yIndex, yVariable in enumerate(oscVariableArray):
        for zIndex, zVariable in enumerate(oscVariableArray):
            currentStep = (xIndex * accuracy ** 2) + (yIndex * accuracy) + zIndex
            print(f'{currentStep * 100 / accuracy**3}%')
            difference = 0
            for [transition, value] in transitionList:
                transitionBox = simpleHarmonicOscillator3DEnergy([xVariable, yVariable, zVariable], [transition[0], transition[1], transition[2]], deuteriumMass) - simpleHarmonicOscillator3DEnergy([xVariable, yVariable, zVariable], [0,0,0], deuteriumMass)
                difference += np.abs(transitionBox - value)
            variableList.append([xVariable, yVariable, zVariable])
            differenceList.append(difference)
            os.system('clear')


differenceList = np.array(differenceList)
variableList = np.array(variableList)

differencesInd = np.argsort(differenceList)

orderedDifference = differenceList[differencesInd]
orderedVariables = variableList[differencesInd]

closestOsc = [simpleHarmonicOscillator3DEnergy(orderedVariables[0], [n,0,0], deuteriumMass) - simpleHarmonicOscillator3DEnergy(orderedVariables[0], [0,0,0], deuteriumMass) for n in range(0,3)]
print(f'Deviance in Simple  Harmonic OScillator: {closestOsc - targetValues[:,0,0]}')
print(f'With Spring constants of: {orderedVariables[0]}')
closestBox = [particleInABox3DEnergy(orderedVariablesBox[0], [n,1,1],deuteriumMass) - particleInABox3DEnergy(orderedVariablesBox[0], [1,1,1],deuteriumMass) for n in range(1,4)]
print(f'Deviance on Particle in a Box: {closestBox - targetValues[:,0,0]}')
print(f'Box Size: {orderedVariablesBox[0]}')


plt.hlines(closestBox, 0, 1)
plt.hlines(closestOsc, 1, 2)
plt.hlines(targetValues[:,0,0], 0, 2, color='r')
plt.show()
