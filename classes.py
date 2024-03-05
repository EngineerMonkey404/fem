
import matrepr
import numpy as np
from matrepr import mprint
from math import sqrt
kPuas=0.3
yunga=2000*1e6
thickness = 0.5

class Element: 
    def __init__(self, nodes, thickness):
        self.stress = 0 # Sigma напряжения
        self.nodes = nodes # узлы 
        self.calculateDetJ() # detJ = x13 - 
        self.calculateB()
        self.calculateD()
        self.area = 0.5 * abs(self.detJ)
        self.thickness = thickness
        self.matrixK = self.thickness * self.area * np.matmul(np.matmul(np.transpose(self.matrixB), self.matrixD), self.matrixB) # t * A * transpose(B) * D * B -- должно быть верно
        self.strain = [] # Эпсилон - линейные деформации


    def calculateStrain(self, matrixQ, dofsToStay):
        # print('MQ', matrixQ)
        mq = {}
        q = []
        for i in range(len(matrixQ)):
            mq[dofsToStay[i]] = matrixQ[i]
        for n in self.nodes: 
            if n.DOFNumber[0] in mq.keys(): 
                q.append(mq[n.DOFNumber[0]])
            else:
                q.append(0)
            if n.DOFNumber[1] in mq.keys(): 
                q.append(mq[n.DOFNumber[1]])
            else: 
                q.append(0)

        # mprint(self.matrixB)
        # mprint(self.matrixD)
        DB = np.matmul(self.matrixD, self.matrixB) 
        self.strain = np.matmul(DB, np.array(q))
        # print('STRAIN', self.strain)
        self.strain = sqrt(self.strain[0]**2 - self.strain[0]*self.strain[1] + self.strain[1]**2 + 3 * self.strain[2]**2)
        # print(self.strain)
        


    def calculateDetJ(self): # Верно
        self.detJ = self.nodeDif(0,1,3) * self.nodeDif(1,2,3) - self.nodeDif(0,2,3) * self.nodeDif(1,1,3) # x13 * y23 - x23 * y13 - Jacobian

    def calculateB(self ): # Верно
        self.matrixB = 1 / self.detJ * np.array([[self.nodeDif(1,2,3), 0, self.nodeDif(1,3,1), 0, self.nodeDif(1,1,2), 0],
                                                 [0,self.nodeDif(0,3,2), 0, self.nodeDif(0,1,3), 0, self.nodeDif(0,2,1)],
                                                 [self.nodeDif(0,3,2), self.nodeDif(1,2,3), self.nodeDif(0,1,3), self.nodeDif(1,3,1), self.nodeDif(0,2,1), self.nodeDif(1,1,2)]])  

    def calculateD(self): # Верно
        self.matrixD = (yunga / (1 - kPuas**2)) * np.array([[1, kPuas, 0],[kPuas, 1, 0],[0, 0, (1-kPuas) / 2]])

    def nodeDif(self,axis, node1, node2): # axis - x или y, и два узла
        return self.nodes[node1 - 1].coordinates[axis] - self.nodes[node2 - 1].coordinates[axis]


class Node:
    def __init__(self, nodeNumber, coordinates, force, pinned):
        self.nodeNumber = nodeNumber # номер узла
        self.coordinates = coordinates # координаты
        self.force = force # силы примененная к узлу
        self.pinned = pinned # закреплен или нет
        self.DOFNumber = [(nodeNumber + 1) * 2 - 1, (nodeNumber + 1) * 2] # Номера степеней свободы, 0ой элемент - по оси x, 1-ый по оси y
