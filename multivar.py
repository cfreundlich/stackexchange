# import matplotlib.pyplot as plt
import math
import itertools
import numpy as np
import tabulate


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



p = math.pi
interesting_points = {'0': 0, 'pi/2': p/2, 'pi': p, '3pi/2': 3*p/2}
permutations = itertools.permutations(interesting_points.items(), 3)
table = [['x', 'y', 'z']]
for (kx, x), (ky, y), (kz, z) in permutations:
    # print(f'{x}, {y}, {z}')
    g = np.array([dfx(x, y ,z), dfy(x, y ,z), dfz(x, y ,z)])
    if np.linalg.norm(g) < 1e-6:
        print(f'{kx}, {ky}, {kz}')
        table.append([kx, ky, kz])

print(tabulate.tabulate(table[1:], headers=table[0]))