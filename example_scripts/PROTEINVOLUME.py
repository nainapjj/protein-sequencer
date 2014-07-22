def side(((xa, ya, za), (xb, yb, zb), (xc, yc, zc), (xd, yd, zd))):
    vectorAa= (xa-xb)
    vectorAb= (ya-yb)
    vectorAc= (za-zb)
    vectorBa= (xa-xc)
    vectorBb= (ya-yc)
    vectorBc= (za-zc)
    vectorCa= (xb-xc)
    vectorCb= (yb-yc)
    vectorCc= (zb-zc)
    vectorDa= (xa-xd)
    vectorDb= (ya-yd)
    vectorDc= (za-zd)
    vectorEa= (xb-xd)
    vectorEb= (yb-yd)
    vectorEc= (zb-zd)
    vectorFa= (xc-xd)
    vectorFb= (yc-yd)
    vectorFc= (zc-zd)
    baseA= float(np.power(vectorAa, 2) + np.power(vectorAb, 2))
    U= float(baseA + np.power(vectorAc, 2))
    baseB= float(np.power(vectorBa, 2) + np.power(vectorBb, 2))
    v= float(baseB + np.power(vectorBc, 2))
    baseC= float(np.power(vectorCa, 2) + np.power(vectorCb, 2))
    w= float(baseC + np.power(vectorCc, 2))
    baseD= float(np.power(vectorDa, 2) + np.power(vectorDb, 2))
    W= float(baseD + np.power(vectorDc, 2))
    baseE= float(np.power(vectorEa, 2) + np.power(vectorEb, 2))
    V= float(baseE + np.power(vectorEc, 2))
    baseF= float(np.power(vectorFa, 2) + np.power(vectorFb, 2))
    u= float(baseF + np.power(vectorFc, 2))
    M= numpy.matrix([[0, u, v, w, 1], [u, 0, W, V, 1], [v, W, 0, U, 1], [w, V, U, 0, 1], [1, 1, 1, 1, 0]])
    Z= M
    volume = (det(Z)/288)**.5
    return volume