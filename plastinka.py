from plastinka_functions import recalculateStress 
from classes import thickness 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import sys


xRange = 20
yRange = 100
force = -1000
numberOfNodes = 3

np.set_printoptions(threshold=sys.maxsize) # Вывод всего np.array без пропуска элементов

fig, (ax, ax1, ax2) = plt.subplots(1,3)

cax = fig.add_axes([0.8, 0.25, 0.05, 0.65])
axSlider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
forceSlider = fig.add_axes([0.25, 0.05, 0.65, 0.03])
fig.subplots_adjust(bottom=0.25, right=0.75)
nodeAmountSlider = Slider(axSlider, 'number of nodes', 2, xRange, valinit=numberOfNodes, valstep=1)
forceSlider = Slider(forceSlider, 'force', 0, 5_000_000, valinit=force, valstep=10000)

recalculateStress(numberOfNodes, force, thickness, xRange, yRange,ax,cax, fig, ax1, ax2)

def updateNumberOfNodes(val):
    global numberOfNodes,ax, line, fig, cax
    numberOfNodes = val
    recalculateStress(numberOfNodes, force,thickness, xRange, yRange, ax,cax,fig, ax1, ax2)


def updateForce(val):
    global numberOfNodes, force, ax, line, fig, cax
    force = -val
    recalculateStress(numberOfNodes, force, thickness, xRange, yRange, ax,cax, fig, ax1, ax2)

nodeAmountSlider.on_changed(updateNumberOfNodes)
forceSlider.on_changed(updateForce)

plt.show()
