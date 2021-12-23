import numpy as np
import matplotlib.pyplot as plt
import os
from wave_functions import *

boxVarMin, boxVarMax = 0.4, 0.7
oscVarMin, oscVarMax = 1, 4
accuracy = 50

targetValues = np.zeros((3,3,3))
targetValues[1,0,0] = 73
targetValues[0,1,0] = 80.2
targetValues[0,0,1] = 101.1
targetValues[2,0,0] = 159.3
targetValues[1,1,0] = 160.2
targetValues[1,0,1] = 173.3
targetValues[0,2,0] = 181.9
targetValues[0,1,1] = 184.1
targetValues[0,0,2] = 219.4

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
                transitionBox = particleInABox3DEnergy([xVariable, yVariable, zVariable], [transition[0]+1, transition[1]+1, transition[2]+1], h2Mass) - particleInABox3DEnergy([xVariable, yVariable, zVariable], [1,1,1], h2Mass)
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
                transitionBox = simpleHarmonicOscillator3DEnergy([xVariable, yVariable, zVariable], [transition[0], transition[1], transition[2]], h2Mass) - simpleHarmonicOscillator3DEnergy([xVariable, yVariable, zVariable], [0,0,0], h2Mass)
                difference += np.abs(transitionBox - value)
            variableList.append([xVariable, yVariable, zVariable])
            differenceList.append(difference)
            os.system('clear')


differenceList = np.array(differenceList)
variableList = np.array(variableList)

differencesInd = np.argsort(differenceList)

orderedDifference = differenceList[differencesInd]
orderedVariables = variableList[differencesInd]

closestOsc = [simpleHarmonicOscillator3DEnergy(orderedVariables[0], [n,0,0], h2Mass) - simpleHarmonicOscillator3DEnergy(orderedVariables[0], [0,0,0], h2Mass) for n in range(0,3)]
print(f'Deviance in Simple  Harmonic OScillator: {closestOsc - targetValues[:,0,0]}')
print(f'With Spring constants of: {orderedVariables[0]}')
closestBox = [particleInABox3DEnergy(orderedVariablesBox[0], [n,1,1],h2Mass) - particleInABox3DEnergy(orderedVariablesBox[0], [1,1,1],h2Mass) for n in range(1,4)]
print(f'Deviance on Particle in a Box: {closestBox - targetValues[:,0,0]}')
print(f'Box Size: {orderedVariablesBox[0]}')


plt.hlines(closestBox, 0, 1)
plt.hlines(closestOsc, 1, 2)
plt.hlines(targetValues[:,0,0], 0, 2, color='r')
plt.show()
