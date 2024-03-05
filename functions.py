from scipy.spatial import Delaunay
import numpy as np
from classes import Element, Node
from matrepr import mdisplay, mprint

def fillElements(triCoordinates, triSimplices, thickness, force, borderX, borderY):
    seenNodes = {}
    elements = []
    for i in range(len(triSimplices)):
        nodes = []
        for j in range(3):
            nodeCoordinates = triCoordinates[i][j]
            node = None
            if seenNodes.get(triSimplices[i,j]):
                nodes.append(seenNodes.get(triSimplices[i,j]))
                continue
            elif nodeCoordinates[0] == borderX and (nodeCoordinates[1] == 0 or nodeCoordinates[1] == borderY): # right top and bot nodes
                node = Node(triSimplices[i][j], nodeCoordinates, force, True)
            elif nodeCoordinates[0] == borderX and not (nodeCoordinates[1] == 0 or nodeCoordinates[1] == borderY): # right pillar with forces
                node = Node(triSimplices[i][j], nodeCoordinates, force, False)
            elif nodeCoordinates[0] == 0: # left pillar
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True)             
            elif nodeCoordinates[1] == 0 and nodeCoordinates[0] != 0 and nodeCoordinates[0] != borderX: # bottom w/0 first and last
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True) 
            elif nodeCoordinates[1] == borderY and nodeCoordinates[0] != 0 and nodeCoordinates[0] != borderX: # top w/0 first and last 
                node = Node(triSimplices[i][j], nodeCoordinates, 0, True) 
            else: 
                node = Node(triSimplices[i][j], nodeCoordinates, 0, False)
            # print('node', node.nodeNumber)
            nodes.append(node)
            seenNodes[triSimplices[i][j]] = node
        e = Element(nodes, thickness)
        elements.append(e)

    # print("seenNodes", seenNodes)
    return elements

def calculateQMatrix(globalMatrix, matrixF):
    # return np.matmul(np.transpose(globalMatrix), matrixF)
    return np.linalg.solve(globalMatrix, matrixF)


def sortDOFS(nodes): 
    result = []
    for i in range(len(nodes)):
        # print(nodes[i].DOFNumber)
        result.append(nodes[i].DOFNumber[0])
        result.append(nodes[i].DOFNumber[1])
    result.sort()
    return result

def fillGlobalMatrix(elements, matrix):
    print(len(elements))
    for i in range(len(elements)):
        element = elements[i]
        DOFS = sortDOFS(element.nodes)
        # print('el', element)
        # print('DFS', DOFS)
        mprint(element.matrixK)
        for j in range(len(element.matrixK)):
            for k in range(len(element.matrixK[j])):
                matrix[DOFS[j] - 1, DOFS[k] - 1] += element.matrixK[j][k]
    print("MATRICA")
    mprint(matrix)
    return matrix


def fillFMatrix(elements, dofsToStay):
    DOFTS = list(dofsToStay)
    DOFTS.sort()
    seen = set()
    matrixF = []
    for e in elements: 
        for n in e.nodes:
            if n.DOFNumber[0] in DOFTS and not (n.DOFNumber[0] in seen):
                matrixF.append(n.force)
            if n.DOFNumber[1] in DOFTS and not (n.DOFNumber[1] in seen):
                matrixF.append(n.force)
            seen.add(n.DOFNumber[1])
            seen.add(n.DOFNumber[0])
    return matrixF



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
    # print('DOFSOTDELETE',numbersToDelete)
    for number in numbersToDelete:
        # print(number)
        globalMatrix = np.delete(globalMatrix, number - 1 - deleted, 0)
        globalMatrix = np.delete(globalMatrix, number - 1 - deleted, 1)
        deleted += 1
    return globalMatrix


def printNodes(elements):
    for element in elements:
        # print('K_MATRIX', element.matrixK)
        DOFS = sortDOFS(element.nodes)
        mprint(element.matrixK, title=None, row_labels = DOFS, col_labels=DOFS, fill_value="--", num_after_dots=100)
        for node in element.nodes:
            print('node', node.pinned, node.force,node.coordinates, node.DOFNumber)



def recalculateStress(numberOfNodes, force,thickness, borderX, borderY, ax,cax,fig,):
    print(force)
    newPoints = []
    basePoints = np.linspace(0,100,numberOfNodes)
    globalMatrix = np.zeros([numberOfNodes*numberOfNodes*2,numberOfNodes* numberOfNodes*2])
    for i in range(basePoints.size):
        for j in range(basePoints.size):
            points = newPoints.append([basePoints[i] , basePoints[j]])
    points =  np.array(newPoints)
    ax.cla()
    line, = ax.plot(points[:,0], points[:,1], 'o')
    line.set_data(points[:, 0], points[:, 1])
    tri = Delaunay(points)
    ax.triplot(points[:,0], points[:,1], tri.simplices)
    ax.plot(points[:,0], points[:,1], 'o')
    elements = fillElements(points[tri.simplices], tri.simplices,thickness, force, borderX, borderY)
    globalMatrix = fillGlobalMatrix(elements, globalMatrix)

    # mprint(globalMatrix, max_rows=18, max_cols=18)
    DOFS = findPinnedDOFS(elements)
    globalMatrix = deletePinnedRowsAndColumns(DOFS[0], globalMatrix)

    # print('points', points)
    # print('trianglus', tri.simplices)
    # print('node coordinates', points[tri.simplices])
    # printNodes(elements)

    # mprint(globalMatrix,row_labels=list(DOFS[1]), col_labels=list(DOFS[1]), max_rows=18, max_cols=18)
    matrixF = fillFMatrix(elements, DOFS[1])
    # mprint(matrixF, col_labels=list(DOFS[1]))
    matrixQ = calculateQMatrix(globalMatrix, matrixF)
    # mprint(matrixQ, col_labels=list(DOFS[1]))
    strains = []
    maxStrain = 0
    for i in range(len(elements)):
        elements[i].calculateStrain(matrixQ, list(DOFS[1]))
        strains.append(elements[i].strain)
    strains[0] = 0
    cax.clear()
    tpc = ax.tripcolor(points[:, 0], points[:,1],strains, cmap="jet")
    fig.colorbar(tpc,ax=ax,cax=cax )
    # print('global', globalMatrix)
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
