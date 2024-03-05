from functions import sortDOFS, fillGlobalMatrix, findPinnedDOFS, deletePinnedRowsAndColumns, printNodes, calculateQMatrix, fillFMatrix, fillElements, recalculateStress 
from classes import thickness 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import sys


xRange = 100
force = 1000
elements = []
points = []
numberOfNodes = 3

np.set_printoptions(threshold=sys.maxsize) # Вывод всего np.array без пропуска элементов

globalMatrix = np.zeros([numberOfNodes*numberOfNodes * 2,numberOfNodes*numberOfNodes * 2]) ## Заполнение глобальной матрицы нулями
basePoints = np.linspace(0,xRange,numberOfNodes)
fig, ax = plt.subplots()

cax = fig.add_axes([0.8, 0.25, 0.05, 0.65])
axSlider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
forceSlider = fig.add_axes([0.25, 0.05, 0.65, 0.03])
fig.subplots_adjust(bottom=0.25, right=0.75)
nodeAmountSlider = Slider(axSlider, 'number of nodes', 2, xRange, valinit=numberOfNodes, valstep=1)
forceSlider = Slider(forceSlider, 'force', 100, 100000, valinit=force, valstep=100)

recalculateStress(numberOfNodes, force, thickness, xRange, xRange,ax,cax, fig)
# for i in range(basePoints.size):
#     for j in range(basePoints.size):
#       points.append([basePoints[i] , basePoints[j]])
#
# points = np.array(points)
# # print(points)
# tri = Delaunay(points) # Триангуляция
# mesh = ax.triplot(points[:,0], points[:,1], tri.simplices)
#
#
# elements = fillElements(points[tri.simplices], tri.simplices, thickness, force, xRange, xRange)
# # printNodes(elements)
# globalMatrix = fillGlobalMatrix(elements, globalMatrix)
#
# DOFS = findPinnedDOFS(elements)
#
# # print('DOFS_TO_STAY', DOFS[1])
# matrixF = fillFMatrix(elements, DOFS[1])
# # print('F', matrixF) 
# globalMatrix = deletePinnedRowsAndColumns(DOFS[0], globalMatrix)
#
# matrixQ = calculateQMatrix(globalMatrix, matrixF)
# # print('q',matrixQ) 
# # print('global', globalMatrix)
# maxStrain = 0
# strains = []
# for i in range(len(elements)):
#     elements[i].calculateStrain(matrixQ, list(DOFS[1]))
#     strains.append(elements[i].strain[0])
#     if maxStrain < elements[i].strain[0]: 
#         maxStrain = elements[i].strain[0]
#
#
# tpc = ax.tripcolor(points[:, 0], points[:,1],strains, cmap="gist_rainbow_r")
# fig.colorbar(tpc)
# for e in elements:
#     x = [] 
#     for n in e.nodes: 
#         x.append(n.coordinates)
#     t = plt.Polygon(x, color='red')
#     ax.gca().add_patch(t)
# print('max', maxStrain)
def updateNumberOfNodes(val):
    global numberOfNodes,points, elements, globalMatrix,ax, line, fig
    numberOfNodes = val
    recalculateStress(numberOfNodes, force,thickness, xRange, xRange, ax,cax,fig)


def updateForce(val):
    global numberOfNodes, force, ax, line, fig
    force = val
    recalculateStress(numberOfNodes, force, thickness, xRange, xRange, ax,cax, fig)

nodeAmountSlider.on_changed(updateNumberOfNodes)
forceSlider.on_changed(updateForce)

plt.show()
