
import numpy as np


def calculateQMatrix(globalMatrix, matrixF):
    print('G', globalMatrix)
    return np.matmul(np.transpose(globalMatrix), matrixF)
    # return np.linalg.solve(globalMatrix, matrixF)


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


def fillFMatrix(elements, dofsToStay):
    dofsToStay = list(dofsToStay)
    seen = set()
    print('DFST', dofsToStay)
    matrixF = []
    for e in elements: 
        for n in e.nodes:
            if n.DOFNumber[0] in dofsToStay and not (n.DOFNumber[0] in seen):
                matrixF.append(n.force)
            if n.DOFNumber[1] in dofsToStay and not (n.DOFNumber[1] in seen):
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
