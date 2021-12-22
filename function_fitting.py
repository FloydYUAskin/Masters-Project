import numpy as np
import matplotlib.pyplot as plt
from wave_functions import *
import warnings

def closestVariableEnergy(function, variableSet, nSet, accuracy, target):
    xVariableSet = np.linspace(variableSet[0][0], variableSet[0][1], accuracy)
    yVariableSet = np.linspace(variableSet[1][0], variableSet[1][1], accuracy)
    zVariableSet = np.linspace(variableSet[2][0], variableSet[2][1], accuracy)

    variableMesh = np.array(np.meshgrid(xVariableSet, yVariableSet, zVariableSet))
    resultsMesh = np.ones_like(variableMesh[0])

    for xIndex, xValue in enumerate(xVariableSet):
        for yIndex, yValue in enumerate(yVariableSet):
            for zIndex, zValue in enumerate(zVariableSet):
                testValues = function([xValue, yValue, zValue], nSet, h2Mass)
                difference = np.abs(testValues - target)
                totalDifference = np.nansum(difference)

                resultsMesh[xIndex, yIndex, zIndex] = totalDifference


    minimumIndex = np.where(resultsMesh == np.nanmin(resultsMesh))


    xClosestVariable = variableMesh[0][minimumIndex]
    yClosestVariable = variableMesh[1][minimumIndex]
    zClosestVariable = variableMesh[2][minimumIndex]

    return resultsMesh, [xClosestVariable, yClosestVariable, zClosestVariable], minimumIndex, [xVariableSet, yVariableSet, zVariableSet]


def variableSearchEnergy(function, variableSet, nSet, accuracy, depth, target):
    resultsMesh, variableList, minimumIndex, variableSets = closestVariableEnergy(function, variableSet, nSet, accuracy, target)
    print(minimumIndex)
    for iteration in range(depth):
        xMinimumIndex = minimumIndex[0][0] - 1
        xMaximumIndex = minimumIndex[0][0] + 1
        yMinimumIndex = minimumIndex[1][0] - 1
        yMaximumIndex = minimumIndex[1][0] + 1
        zMinimumIndex = minimumIndex[2][0] - 1
        zMaximumIndex = minimumIndex[2][0] + 1

        if xMaximumIndex >= len(variableSets[0]):
            xMaximumIndex -= 1
        if yMaximumIndex >= len(variableSets[1]):
            yMaximumIndex -= 1
        if zMaximumIndex >= len(variableSets[2]):
            zMaximumIndex -= 1

        xMinimumVariable = variableSets[0][xMinimumIndex]
        xMaximumVariable = variableSets[0][xMaximumIndex]
        yMinimumVariable = variableSets[1][yMinimumIndex]
        yMaximumVariable = variableSets[1][yMaximumIndex]
        zMinimumVariable = variableSets[2][zMinimumIndex]
        zMaximumVariable = variableSets[2][zMaximumIndex]

        variableSet = [xMinimumVariable, xMaximumVariable], [yMinimumVariable, yMaximumVariable], [zMinimumVariable, zMaximumVariable]

        resultsMesh, variableList, minimumIndex, variableSets = closestVariableEnergy(function, variableSet, nSet, accuracy, target)



    return variableList



if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nRange = np.arange(0,3,1)
        nSetOscillator = np.array(np.meshgrid(nRange,nRange,nRange)).astype('float')
        nRange = np.arange(1,4,1)
        nSetBox = np.array(np.meshgrid(nRange,nRange,nRange)).astype('float')

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

        targetValues[targetValues==0] = np.nan
        nSetOscillator[targetValues==0] = np.nan
        nSetBox[targetValues==0] = np.nan



        closestVariablesParticleInABox = variableSearchEnergy(particleInABox3DEnergy, [[0,10],[0,10],[0,10]], nSet, 10, 10, targetValues)
        closestVariablesSimpleHarmonicOscillator = variableSearchEnergy(simpleHarmonicOscillator3DEnergy, [[0,10],[0,10],[0,10]], nSet, 10, 10, targetValues)
        print(closestVariablesParticleInABox, closestVariablesSimpleHarmonicOscillator)

        approximateBox = particleInABox3DEnergy(closestVariablesParticleInABox, nSet, h2Mass)
        approximateOscillator = simpleHarmonicOscillator3DEnergy(closestVariablesSimpleHarmonicOscillator, nSet, h2Mass)

        plt.hlines(approximateBox[:,0,0],0,1)
        plt.hlines(approximateOscillator[:,0,0],1,2)
        plt.hlines(targetValues[:,0,0],0,2,color='r')
        plt.show()

        print(targetValues[:,0,0], approximateOscillator[:,0,0])
