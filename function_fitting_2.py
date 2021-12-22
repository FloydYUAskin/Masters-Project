import numpy as np
import matplotlib.pyplot as plt
from wave_functions import *
import warnings

def createResultsMesh(function, variableSet, nSet, accuracy, target, m):
    xVariableArray = np.linspace(variableSet[0][0], variableSet[0][1], accuracy)
    yVariableArray = np.linspace(variableSet[1][0], variableSet[1][1], accuracy)
    zVariableArray = np.linspace(variableSet[2][0], variableSet[2][1], accuracy)


    resultsMesh = np.zeros((accuracy, accuracy, accuracy))

    for xIndex, xVariable in enumerate(xVariableArray):
        for yIndex, yVariable in enumerate(yVariableArray):
            for zIndex, zVariable in enumerate(zVariableArray):

                simulatedResults = function([xVariable, yVariable, zVariable], nSet, m)

                differenceArray = np.abs(simulatedResults - target)
                totalDifference = np.nansum(differenceArray)

                resultsMesh[xIndex, yIndex, zIndex] = totalDifference

    return resultsMesh, (xVariableArray, yVariableArray, zVariableArray)

def createVariableSet(resultsMesh, variableLists):
    minimumResultIndex = np.argmin(resultsMesh)
    print(minimumResultIndex)



if __name__ == '__main__':
    nRangeOscillator = np.linspace(0,3,1)
    nRangeBox = np.linspace(1,4,1)

    nSetOscillator = np.meshgrid(nRangeOscillator, nRangeOscillator, nRangeOscillator)
    nSetBox = np.array(np.meshgrid(nRangeBox, nRangeBox, nRangeBox))
    print(nSetBox.shape)

    targetValues = particleInABox3DEnergy([5,5,5], nSetBox, h2Mass)

    results, variables = createResultsMesh(particleInABox3DEnergy, [[0,10], [0,10], [0,10]], nSetBox, 10, targetValues, h2Mass)
    createVariableSet(results, variables)
