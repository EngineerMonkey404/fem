    
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.multiarray import dtype
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation
from matplotlib.widgets import Slider
from random import randrange


class Element: 
    def __init__(self, nodes, thickness):
        self.stress = 0
        self.nodes = nodes
        self.calculateDetJ()
        self.calculateB()
        self.calculateD()
        self.area = 0.5 * self.detJ
        self.thickness = thickness
        self.matrixK = self.thickness * self.area * np.matmul(np.matmul(np.transpose(self.matrixB), self.matrixD), self.matrixB)

    def calculateDetJ(self):
        self.detJ = self.nodeDif(0,1,3) * self.nodeDif(1,2,3) - self.nodeDif(0,2,3) * self.nodeDif(1,1,3)

    def calculateB(self ):
        self.matrixB = 1 / self.detJ * np.array([[self.nodeDif(1,2,3), 0, self.nodeDif(1,3,1), 0, self.nodeDif(1,1,2), 0],[0,self.nodeDif(0,3,2), 0, self.nodeDif(0,1,3), 0, self.nodeDif(0,2,1)],[self.nodeDif(0,3,2), self.nodeDif(1,2,3), self.nodeDif(0,1,3), self.nodeDif(1,3,1), self.nodeDif(0,2,1), self.nodeDif(1,1,2)]]) 

    def calculateD(self):
        self.matrixD = yunga / (1 - kPuas**2) * np.array([[1, kPuas, 0],[kPuas, 1, 0],[0, 0, (1-kPuas) / 2]])

    def nodeDif(self,axis, node1, node2):
        return self.nodes[node1 - 1].coordinates[axis] - self.nodes[node2 - 1].coordinates[axis]


class Node:
    def __init__(self, nodeNumber, coordinates, force, pinned):
        self.nodeNumber = nodeNumber
        self.coordinates = coordinates
        self.force = force
        self.pinned = pinned
        self.DOFNumber = [(nodeNumber + 1) * 2 - 1, (nodeNumber + 1) * 2] 

kPuas=0.3
yunga=2000*1e6
thickness = 0.5
xRange = 100
force = 1000
elements = []
points = []
numberOfNodes = 3


def fillElements(triCoordinates, triSimplices):
    global elements, thickness
    seenNodes = {}
    for i in range(len(triSimplices)):
        nodes = []
        for j in range(3):
            nodeCoordinates = triCoordinates[i][j]
            node = None
            if seenNodes.get(triSimplices[i,j]):
                nodes.append(seenNodes.get(triSimplices[i,j]))
                continue
            elif nodeCoordinates[0] == xRange and (nodeCoordinates[1] == 0 or nodeCoordinates[1] == xRange): # right top and bot nodes
                node = Node(triSimplices[i][j], nodeCoordinates, force, True)
            elif nodeCoordinates[0] == xRange and not (nodeCoordinates[1] == 0 or nodeCoordinates[1] == xRange): # right pillar with forces
                node = Node(triSimplices[i][j], nodeCoordinates, force, False)
            elif nodeCoordinates[0] == 0: # left pillar
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True)             
            elif nodeCoordinates[1] == 0 and nodeCoordinates[0] != 0 and nodeCoordinates[0] != xRange: # bottom w/0 first and last
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True) 
            elif nodeCoordinates[1] == xRange and nodeCoordinates[0] != 0 and nodeCoordinates[0] != xRange: # top w/0 first and last
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True)
            else: 
                node = Node(triSimplices[i][j], nodeCoordinates, 0, False)
            # print('node', node.nodeNumber)
            nodes.append(node)
            seenNodes[triSimplices[i][j]] = node
        e = Element(nodes, thickness)
        elements.append(e)

    print("seenNodes", seenNodes)
    return elements

def sortDOFS(nodes): 
    result = []
    for i in range(len(nodes)):
        print(nodes[i].DOFNumber)
        result.append(nodes[i].DOFNumber[0])
        result.append(nodes[i].DOFNumber[1])
    result.sort()
    return result

def fillGlobalMatrix(elements, matrix):
    print('le', len(elements))
    for i in range(len(elements)):
        element = elements[i]
        DOFS = sortDOFS(element.nodes)
        print('DOFS', DOFS)
        for j in range(len(element.matrixK)):
            for k in range(len(element.matrixK[j])):
                matrix[DOFS[j] - 1, DOFS[k] - 1] += element.matrixK[j][k]
        return matrix


def findPinnedDOFS(elements):
    DOFSToDelete = set()
    DOFSToStay = set()
    for element in elements:
        for node in element.nodes:
            if node.pinned: 
                DOFSToDelete.add(node.DOFNumber[0])
                DOFSToDelete.add(node.DOFNumber[1])
            else:
                DOFSToStay.add(node.DOFNumber[0])
                DOFSToStay.add(node.DOFNumber[1])
    return [DOFSToDelete, DOFSToStay]

def deletePinnedRowsAndColumns(numbersToDelete,globalMatrix):
    deleted = 0
    print('DOFSOTDELETE',numbersToDelete)
    for number in numbersToDelete:
        print(number)
        globalMatrix = np.delete(globalMatrix, number - 1 - deleted, 0)
        globalMatrix = np.delete(globalMatrix, number - 1 - deleted, 1)
        deleted += 1
    return globalMatrix


def printNodes(elements):
    for element in elements:
        for node in element.nodes:
            print('node', node.pinned, node.coordinates, node.DOFNumber)


globalMatrix = np.empty([numberOfNodes*numberOfNodes * 2,numberOfNodes*numberOfNodes * 2])
basePoints = np.linspace(0,xRange,numberOfNodes)
fig, ax = plt.subplots()
axSlider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
fig.subplots_adjust(left=0.25, bottom=0.25)
nodeAmountSlider = Slider(axSlider, 'number of nodes', 2, xRange, valinit=numberOfNodes, valstep=1)

for i in range(basePoints.size):
    for j in range(basePoints.size):
      points.append([basePoints[i] , basePoints[j]])

points = np.array(points)
# print(points)
tri = Delaunay(points)
mesh = ax.triplot(points[:,0], points[:,1], tri.simplices)
line, = ax.plot(points[:,0], points[:,1], 'o')


# print('points', points)
# print('trianglus', tri.simplices)
# print(points[tri.simplices])
elements = fillElements(points[tri.simplices], tri.simplices)
printNodes(elements)
for i in range(len(elements)): 

    print('element', i) 
    for j in range(3):
        print('nodes', elements[i].nodes[j].DOFNumber, elements[i].nodes[j].coordinates)
print('empty global', globalMatrix)
globalMatrix = fillGlobalMatrix(elements, globalMatrix)
print('global before', globalMatrix)

DOFS = findPinnedDOFS(elements)
globalMatrix = deletePinnedRowsAndColumns(DOFS[0], globalMatrix)
print('global', globalMatrix)
def updateSlider(val):
    global points, elements, globalMatrix
    newPoints = []
    print('updated', val)
    basePoints = np.linspace(0,100,val)
    for i in range(basePoints.size):
        for j in range(basePoints.size):
            points = newPoints.append([basePoints[i] , basePoints[j]])
    points =  np.array(newPoints)
    ax.cla()
    line.set_data(points[:, 0], points[:, 1])
    tri = Delaunay(points)
    ax.triplot(points[:,0], points[:,1], tri.simplices)
    ax.plot(points[:,0], points[:,1], 'o')
    # print('triangles', tri.simplices)
    # print('triangles coordinates', points[tri.simplices]) 
    elements = fillElements(points[tri.simplices], tri.simplices)
    globalMatrix = fillGlobalMatrix(elements, globalMatrix)

    DOFS = findPinnedDOFS(elements)
    globalMatrix = deletePinnedRowsAndColumns(DOFS[0], globalMatrix)
    print('global', globalMatrix)
    fig.canvas.draw_idle()
nodeAmountSlider.on_changed(updateSlider)

plt.show()
