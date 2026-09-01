import numpy as np


## se usa la matriz de rotacion para rotar
def rota(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])

## la matriz escala la componenete i de un vector Rn en un factor s[i] ==> son las componentes de la diagonal, las reemplazo en un range de ceros
def escala(s):
    dim = np.shape(s)[0]
    m = np.zeros((dim,dim))
    for i in np.arange(dim):
        m[i,i] = s[i]
    return m 

'''
def rota_y_escala1(theta, s):
    dim = np.shape(s)[0]
    m = np.zeros((dim,dim))
    for i in np.arange(dim):
        m = np.array([[s[0]*np.cos(theta), -s[0]*np.sin(theta)],s[1]*[np.sin(theta),  s[1]*np.cos(theta)]])
        m[i,i] = s[i]
    return m
'''


## transformcion de matrices es composicion ==> 1°roto 2° escalo <=> E(R(v))
def rota_y_escala2(theta,s):
    m = np.array([[s[0]*np.cos(theta), -s[0]*np.sin(theta)],
                  [s[1]*np.sin(theta),  s[1]*np.cos(theta)]])
    return m

print(rota_y_escala2(0, [2,3]))


def rota_y_escala(theta, s):
    r = rota(theta)
    e = escala(s) 
    m = e@r
    return m




## 1°rota 2°escalo 3° mueve es otra composicion de matrices --> M(E(R(v))) , uso la matriz de rot escalada y le agrego la columna 3 para moverla 
def afin(theta, s, b):
    m = np.array([[s[0]*np.cos(theta), -s[0]*np.sin(theta), b[0]],
                 [s[1]*np.sin(theta),  s[1]*np.cos(theta), b[1]], 
                 [0                  , 0                    , 1]])
    return m

## devuelvo el vector w que sale de resolver m.v = w, luego w0 = m0.v0 y w1= m1.v1                  
def trans_afin(v, theta, s, b):
    w = np.array([
        s[0]*(v[0]*np.cos(theta) - v[1]*np.sin(theta)) + b[0],
        s[1]*(v[0]*np.sin(theta) + v[1]*np.cos(theta)) + b[1]
    ])
    return w






