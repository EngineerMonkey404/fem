import numpy as np
import matplotlib.pyplot as plt
import sys
from sympy import *
from matplotlib.tri import Triangulation
I = 0.00013
# I = 0.003
E = 65e6
EI = E*I
F = 3043478

print(EI)
np.set_printoptions(threshold=sys.maxsize)
fig, ax = plt.subplots()
fig1, ax1 = plt.subplots()

xRange = 0.02
yRange = 0.1

yNodes = 100 
xNodes = 100

yPoints = np.linspace(0, yRange, yNodes)
xPoints = np.linspace(0, xRange, xNodes)

points = []

for xp in xPoints:
    for yp in yPoints: 
        points.append([xp,yp])


points = np.array(points)


y = Symbol('y')
a0 = Symbol('a0') 
a1 = Symbol('a1') 
a2 = Symbol('a2') 
a3 = Symbol('a3') 
a4 = Symbol('a4') 
a5 = Symbol('a5') 
a6 = Symbol('a6') 

xApprox = a2*y**2 + a3*y**3 + a4*y**4 + a5*y**5 + a6*y**6


print(xApprox)
# print(diff(xApprox, y))
# print(diff(xApprox,y,2))
# print(integrate(diff(xApprox, x, 2), (y,0, yRange)))
ISE = integrate((EI / 2) * diff(xApprox, y, 2) ** 2,(y, 0, yRange))

W = integrate(F * xApprox, (y,0, yRange))

PE = ISE - W

Eq1 = diff(PE, a2)
Eq2 = diff(PE, a3)
Eq3 = diff(PE, a4)
Eq4 = diff(PE, a5)
Eq5 = diff(PE, a6)

solution = solve([Eq1, Eq2, Eq3, Eq4, Eq5], [a2, a3,a4,a5,a6])

# print(ISE)
# [print(eq) for eq in [Eq1, Eq2, Eq3, Eq4, Eq5]]
# print(solution)

print(xApprox.subs(solution))
stressPoints = []
for p in points:
    if p[0] == 0 :
        solution['y'] = p[1]
        # print(xApprox.subs(solution))
        stressPoints.append([xApprox.subs(solution), p[1]])
    # if p[0] == xRange:
    #     solution['y'] = p[1]
    #     stressPoints.append([xApprox.subs(solution) + xRange, p[1]])

stressPoints = np.array(stressPoints,dtype='float64')

tri = Triangulation(points[:,0], points[:,1])
ax.triplot(points[:,0], points[:,1], tri.triangles)


# print('points', points)
# print('sPoints', stressPoints)

tri1 = Triangulation(stressPoints[:,0], stressPoints[:,1])
ax1.plot(stressPoints[:,0], stressPoints[:,1])

p1 = plot(xApprox.subs(solution))
p1.show()
plt.show()

