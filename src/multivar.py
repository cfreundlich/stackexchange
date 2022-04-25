import matplotlib.pyplot as plt
import math
import itertools
import numpy as np
import tabulate
from scipy.optimize import fsolve



# solving https://math.stackexchange.com/q/4430779/143884


s = math.sin
c = math.cos
def c2(g):
    return math.cos(g) ** 2
def s2(g):
    return math.sin(g) ** 2


def f(x, y, z):
    A = s(z)
    B = c(y) - s(x)
    C = c(y) * s(z) + c(x)*s(y)*c(z) - s(x)*c(y)
    return A * B * C


def dfx(x, y, z):
    a1 = s(z)*  c(z) 
    a2 = s(x)-c(y)
    a3 = s(x) * s(y)+c(x)
    A = a1 * a2 * a3
    b1 =  -c(x)
    b2 = c(x) * s(y) * c(z) - s(x) * c(z) + c(y) * s(z)
    B = b1 * b2 
    return A + B


def dfy(x, y, z):
    A = s(z)
    
    a21 = c(y)-s(x)
    a22 = c(x) * c(y) * c(z) - s(y) * s(z)
    B = a21 * a22

    b121 = c(z)  
    b222 = c(x) * s(y) - s(x)
    c2 = b121 * b222

    c3 = c(y) * s(z)

    C = s(y) * (c2 + c3)
    
    return A * (B - C)

def dfz(x, y, z):
    
    A = c(y) - s(x)

    b1 = c2(z) * (c(x)*s(y) - s(x))
    b2 = s2(z) * (s(x) - c(x)*s(y))
    b3 =2 * c(y)*c(z)*s(z)

    B = b1 + b2 + b3

    return A * B


def gradient(vars):
    return [dfx(*vars), dfy(*vars), dfz(*vars)]


eps = 1e-6
decimals = 3
table = [['x', 'y', 'z']]

p = math.pi
interesting_points = {'0': 0, 'pi/2': p/2, 'pi': p, '3pi/2': 3*p/2}
permutations = itertools.permutations(interesting_points.items(), 3)
for (kx, x), (ky, y), (kz, z) in permutations:
    if np.linalg.norm(gradient([x, y ,z])) < eps:
        table.append([kx, ky, kz])
        continue
    
solutions = set(
    tuple(round(interesting_points[k], decimals) for k in row) 
    for row in table[1:])

more_interesting_points = np.arange(0, 2 * math.pi, 0.25)
more_permutations = itertools.permutations(more_interesting_points, 3)
for x, y, z in more_permutations:

    solution = fsolve(gradient, (x, y, z))
    if not all(
        0 < solution[i] < 2 * math.pi
        for i in [0, 1, 2]
    ):
        continue
    solution = tuple(round(v, decimals) for v in solution)
    if solution in solutions:
        continue
    solutions.add(solution)
    # table.append(solution)

# print(tabulate.tabulate(table[1:], headers=table[0]))


class Points:
    def __init__(self, solutions_) -> None:
        self._solutions_ = solutions_
    def x(self):
        return [r[0] for r in self._solutions_]
    def y(self):
        return [r[1] for r in self._solutions_]
    def z(self):
        return [r[2] for r in self._solutions_]



points = Points(solutions)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(points.x(), points.y(), points.z())

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim([0, 2*math.pi])
ax.set_ylim([0, 2*math.pi])
ax.set_zlim([0, 2*math.pi])

plt.show()