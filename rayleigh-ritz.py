import numpy as np
import matplotlib.pyplot as plt
import sys
from sympy import *
from matplotlib.tri import Triangulation

I = 0.003
E = 20_000_000
EI = E*I
F = 45

print(EI)
np.set_printoptions(threshold=sys.maxsize)
fig, ax = plt.subplots()
fig1, ax1 = plt.subplots()

xRange = 5
yRange = 100

yNodes = 9 
xNodes = 5


yPoints = np.linspace(0, yRange, yNodes)
xPoints = np.linspace(0, xRange, xNodes)

points = []

for xp in xPoints:
    for yp in yPoints: 
        points.append([xp,yp])


stressPoints = points
points = np.array(points)


x = Symbol('x')
a0 = Symbol('a0') 
a1 = Symbol('a1') 
a2 = Symbol('a2') 
a3 = Symbol('a3') 
a4 = Symbol('a4') 
a5 = Symbol('a5') 
a6 = Symbol('a6') 

yApprox = a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6


print(yApprox)
print(diff(yApprox, x))
print(diff(yApprox,x,2))
print(integrate(diff(yApprox, x, 2), (x,0, yRange)))
ISE = integrate((EI / 2) * diff(yApprox, x, 2) ** 2,(x, 0, yRange))

W = integrate(45000 * yApprox, (x,0, yRange))

PE = ISE - W

Eq1 = diff(PE, a2)
Eq2 = diff(PE, a3)
Eq3 = diff(PE, a4)
Eq4 = diff(PE, a5)
Eq5 = diff(PE, a6)

solution = solve([Eq1, Eq2, Eq3, Eq4, Eq5], [a2, a3,a4,a5,a6])

print(ISE)
[print(eq) for eq in [Eq1, Eq2, Eq3, Eq4, Eq5]]
print(solution)

print(yApprox.subs(solution))


for p in stressPoints:
    if p[1] <= 25:
        continue
    solution['x'] = p[0]
    p[0] = p[0] + yApprox.subs(solution)

stressPoints = np.array(stressPoints,dtype='float64')

tri = Triangulation(points[:,0], points[:,1])
ax.triplot(points[:,0], points[:,1], tri.triangles)


print('points', points)
print('sPoints', stressPoints)

tri1 = Triangulation(stressPoints[:,0], stressPoints[:,1])
ax1.triplot(stressPoints[:,0], stressPoints[:,1], tri1.triangles)

p1 = plot(yApprox.subs(solution))
p1.show()
plt.show()

