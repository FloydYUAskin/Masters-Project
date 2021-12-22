import numpy as np
import matplotlib.pyplot as plt
import os
from wave_functions import *

def createResults(function, variableSets, nSet, accuracy, m, target):
    xVariableArray = np.linspace(variableSets[0][0], variableSets[0][1], accuracy)
    yVariableArray = np.linspace(variableSets[1][0], variableSets[1][1], accuracy)
    zVariableArray = np.linspace(variableSets[2][0], variableSets[2][1], accuracy)

    resultsList = []
    differenceList = []
    variablesList = []

    numberOfSteps = accuracy ** 3

    for xIndex, xVariable in enumerate(xVariableArray):
        for yIndex, yVariable in enumerate(yVariableArray):
            for zIndex, zVariable in enumerate(zVariableArray):
                stepsCompleted = (xIndex * accuracy**2) + (yIndex * accuracy) + (zIndex)
                print(f'{stepsCompleted*100/numberOfSteps}%')
                os.system('clear')
                testResults = function([xVariable, yVariable, zVariable], nSet, m) - function([xVariable, yVariable, zVariable], [nSet[0][0,0,0],nSet[1][0,0,0],nSet[2][0,0,0]], m)
                differenceArray = target - testResults
                totalDifference = np.linalg.norm(differenceArray)
                print(totalDifference)

                resultsList.append(testResults)
                differenceList.append(totalDifference)
                variablesList.append([xVariable,yVariable,zVariable])

    return resultsList, differenceList, variablesList

def nearestVariables(resultsList, differenceList, variablesList):
    resultsList = np.array(resultsList)
    differenceList = np.array(differenceList)
    variablesList = np.array(variablesList)

    differenceListInd = differenceList.argsort()

    differenceListSorted = differenceList[differenceListInd]
    resultsListSorted = resultsList[differenceListInd]
    variablesListSorted = variablesList[differenceListInd]

    return variablesListSorted[0]


if __name__ == '__main__':

    nRangeBox = np.arange(0,3,1)
    nSetBox = np.meshgrid(nRangeBox, nRangeBox, nRangeBox)

    nRangeOsc = np.arange(1,4,1)
    nSetOscillator = np.meshgrid(nRangeOsc, nRangeOsc, nRangeOsc)

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

    nSetBox[0][targetValues==0] = 0
    nSetBox[1][targetValues==0] = 0
    nSetBox[2][targetValues==0] = 0

    nSetOscillator[0][targetValues==0] = 0
    nSetOscillator[1][targetValues==0] = 0
    nSetOscillator[2][targetValues==0] = 0

    # results, differences, variables = createResults(particleInABox3DEnergy, [[1e-1,10], [1e-1,10], [1e-1,10]], nSetBox, 10, h2Mass, targetValues)
    # closestBoxSize = nearestVariables(results, differences, variables)
    #
    # results, differences, variables = createResults(simpleHarmonicOscillator3DEnergy, [[1e-5,1e-3], [1e-5,1e-3], [1e-5,1e-3]], nSetOscillator, 10, h2Mass, targetValues)
    # closestOscillatorStrength = nearestVariables(results, differences, variables)
    #
    # closestBox = particleInABox3DEnergy(closestBoxSize, nSetBox, h2Mass) - particleInABox3DEnergy(closestBoxSize, [0,0,0], h2Mass)
    # closestOscillator = simpleHarmonicOscillator3DEnergy(closestOscillatorStrength, nSetOscillator, h2Mass) - simpleHarmonicOscillator3DEnergy(closestOscillatorStrength, [1,1,1], h2Mass)

    print(targetValues, nSetOscillator)
