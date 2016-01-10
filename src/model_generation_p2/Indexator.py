import numpy


def indexator_new((a, b, c, d)):
    u = abs(a - d) * 5.28
    v = abs(a - c) * 5.28
    w = abs(a - b) * 5.28
    U = abs(b - c) * 5.28
    V = abs(b - d) * 5.28
    W = abs(c - d) * 5.28
    M= numpy.matrix([[0, u, v, w, 1], [u, 0, W, V, 1], [v, W, 0, U, 1], [w, V, U, 0, 1], [1, 1, 1, 1, 0]])
    volume = (numpy.linalg.det(M)/288)**.5
    M=[]
    volume = round(volume, 7)
    return volume


def indexator((a, b, c, d)):
    u = abs(a - d)
    v = abs(a - c)
    w = abs(a - b)
    U = abs(b - c)
    V = abs(b - d)
    W = abs(c - d)
    M= numpy.matrix([[0, u, v, w, 1], [u, 0, W, V, 1], [v, W, 0, U, 1], [w, V, U, 0, 1], [1, 1, 1, 1, 0]])
    volume = (numpy.linalg.det(M)/288)**.5
    M=[]
    volume = round(volume, 7)
    return volume

