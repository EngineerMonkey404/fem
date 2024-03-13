import numpy as np
import matplotlib.pyplot as plt
import sys
from sympy import *
from matplotlib.tri import Triangulation

I = 0.00013
# I = 0.003
E = 2e11
EI = E*I
F = 3043478

zadelka = 0.02
yRange = 0.8
xRange = 0.02
yNodes = 100
xNodes = 100
yPoints = np.linspace(0, yRange, yNodes)
xPoints = np.linspace(0, xRange, xNodes)



points = []

for xp in xPoints:
    for yp in yPoints: 
        points.append([xp,yp])


points = np.array(points)

yPoints = np.linspace(0, yRange, yNodes)
xPoints = np.linspace(0, xRange, xNodes)
fig, ax = plt.subplots()
a0 = Symbol('a0') 
a1 = Symbol('a1') 
a2 = Symbol('a2') 
a3 = Symbol('a3') 
a4 = Symbol('a4') 
a5 = Symbol('a5')
a6 = Symbol('a6') 
y = Symbol('y')
W0 = 1 
W1 = y
W2 = y**2
W3 = y**3
W4 = y**4
W5 = y**5
W6 = y**6

xApprox = a0*W0 + a1*W1 + a2*W2 + a3*W3 + a4*W4+ a5*W5 + a6*W6 

Eq0 = integrate(diff(W0, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W0 * F*y / EI, (y, 0, yRange))
Eq1 = integrate(diff(W1, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W1 * F*y / EI, (y, 0, yRange))
Eq2 = integrate(diff(W2, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W2 * F*y / EI, (y, 0, yRange))
Eq3 = integrate(diff(W3, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W3 * F*y / EI, (y, 0, yRange))

Eq4 = integrate(diff(W4, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W4 * F*y / EI, (y, 0, yRange))
Eq5 = integrate(diff(W5, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W5 * F*y / EI, (y, 0, yRange))
Eq6 = integrate(diff(W6, y) * diff(xApprox, y), (y, 0, yRange)) - integrate(W6 * F*y / EI, (y, 0, yRange))


 
print(Eq1)
print(Eq2)
print(Eq3)
print(Eq4)
print(Eq5)
print(Eq6)
solution = solve([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6], [a1, a2, a3, a4, a5, a6])

solution['a0'] = Eq1
func = xApprox.subs(solution) 
print(func)
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

ax.plot(stressPoints[:,0], stressPoints[:,1])

# p1.show()
plt.show()
