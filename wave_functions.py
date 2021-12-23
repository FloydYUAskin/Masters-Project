import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite, factorial

hbar = 1
e = 1
electronMass = 1
bohrRadius = 1

hartree_to_cm = 2625.5
h2Mass = 1837.15267 * 2
deuteriumMass = 3673.230534 * 2

def particleInABox1D(L, n, Array):
    return np.sin(n * np.pi * Array / L) * np.sqrt(2 / L)

def particleInABox3D(L, n, Position):
    xWave = particleInABox1D(L[0], n[0], Position[0])
    yWave = particleInABox1D(L[1], n[1], Position[1])
    zWave = particleInABox1D(L[2], n[2], Position[2])

    return xWave * yWave * zWave

def particleInABox3DEnergy(L, n, m):
    constant = (2 * np.pi)**2 / (8 * m)
    bracket = (n[0]**2 / L[0]**2) + (n[1]**2 / L[1]**2) + (n[2]**2 / L[2]**2)

    return constant * bracket * hartree_to_cm


def simpleHarmonicOscillator1D(k, n, m, Array):
    omega = np.sqrt(k / m)
    constant = np.sqrt(np.sqrt(m * omega / (np.pi*hbar))) * 1./np.sqrt((2 ** n) * factorial(n))
    hermitePolynomial = hermite(n)
    wave = constant * hermitePolynomial(Array * np.sqrt(m * omega/hbar)) * np.exp(-Array**2 * m * omega / (hbar*2))

    return wave

def simpleHarmonicOscillator3D(k, n, m, Position):
    xWave = simpleHarmonicOscillator1D(k[0], n[0], m, Position[0])
    yWave = simpleHarmonicOscillator1D(k[1], n[1], m, Position[1])
    zWave = simpleHarmonicOscillator1D(k[2], n[2], m, Position[2])

    return xWave * yWave * zWave

def simpleHarmonicOscillator3DEnergy(k, n, m):
    xOmega = np.sqrt(k[0] / m)
    yOmega = np.sqrt(k[1] / m)
    zOmega = np.sqrt(k[2] / m)

    E = ((n[0] + 0.5) * xOmega) + ((n[1] + 0.5) * yOmega) + ((n[2] + 0.5) * zOmega)
    return E * hartree_to_cm

def calculateTransitions(energies):
    energyLevelTransitions = np.zeros(energies[:,0,0].shape)
    for energyLevel in energies[:,0,0]:
        energyTransitions = [energyLevel - energyLevelOther for energyLevelOther in energies[:,0,0]]
        energyLevelTransitions = np.vstack([energyLevelTransitions, energyTransitions])

    energyLevelTransitions = np.delete(energyLevelTransitions,0,axis=1)
    return energyLevelTransitions


if __name__ == '__main__':
    import matplotlib
    matplotlib.style.use('ggplot')

#----------------------------------------------------------------------------Particle In a Box---------------------------------------------------------------------------------

    positionArray = np.linspace(0,1,100)
    positionMesh = np.meshgrid(positionArray, positionArray, positionArray)
    energyLevelArray = np.arange(0,5,1)
    energyLevelMesh = np.meshgrid(energyLevelArray, energyLevelArray, energyLevelArray)
    lengthArray = np.linspace(0,5,5)
    lengthMesh = np.meshgrid(lengthArray, lengthArray, lengthArray)

    fig, axes = plt.subplots(2,2)

    particleBoxDistribution = particleInABox3D([1,1,1], [1,1,1], positionMesh) ** 2
    particleBoxEnergies = particleInABox3DEnergy([1,1,1], energyLevelMesh, h2Mass)
    particleBoxEnergiesLength = particleInABox3DEnergy(lengthMesh, [1,1,1], h2Mass)
    particleBoxDistribution1D = [particleInABox1D(1, n, positionArray) ** 2 for n in range(5)]

    for array in particleBoxDistribution1D:
        axes[0][1].plot(positionArray, array)

    axes[0][0].pcolor(particleBoxDistribution[:,:,50])
    axes[1][0].pcolor(particleBoxDistribution[:,50,:])

    axes[1][1].hlines(particleBoxEnergiesLength[:,1,1], -1, 1)

    plt.show()

#-------------------------------------------------------------------------------Simple Harmonic Oscillator--------------------------------------------------------------------

    positionArray = np.linspace(-3, 3, 100)
    positionMesh = np.meshgrid(positionArray, positionArray, positionArray)
    energyLevelArray = np.arange(0,5,1)
    energyLevelMesh = np.meshgrid(energyLevelArray, energyLevelArray, energyLevelArray)
    kArray = np.linspace(0,10,10)
    kMesh = np.meshgrid(kArray, kArray, kArray)

    fig, axes = plt.subplots(2,2)

    simpleHarmonics1D = np.array([simpleHarmonicOscillator1D(k,0,1,positionArray) ** 2 for k in range(5)])
    for array in simpleHarmonics1D:
        axes[0][0].plot(positionArray, array)

    simpleHarmonics3D = simpleHarmonicOscillator3D([0.5,0.5,0.5], [0,0,0], h2Mass, positionMesh) ** 2
    simpleHarmonics3DEnergy = simpleHarmonicOscillator3DEnergy(kMesh, [0,0,0], h2Mass)

    axes[1][0].pcolor(simpleHarmonics3D[:,:,50])
    axes[0][1].pcolor(simpleHarmonics3D[:,50,:])
    axes[1][1].hlines(simpleHarmonics3DEnergy[:,0,0], -1, 1, colors = 'r')

    plt.show()

    transitionEnergies = calculateTransitions(simpleHarmonics3DEnergy)
    plt.hlines(transitionEnergies[1],0,1)
