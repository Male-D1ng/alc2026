import numpy as np

def norma(x,p):
    suma= 0
    if p != 'inf':
        for i in x:
            suma = suma+(absoluto(i))**p 
        return suma**(1/p)   
    else:
        maximo = 0
        for i in x:
            if absoluto(i)>= maximo:
                maximo = absoluto(i)
        return maximo



def absoluto(x):
    if x <0:
        return -x
    else:
        return x




def normaliza(x,p):
    lista = []
    for vector in x:
        res = norma(vector,p)
        lista.append(res)
    return lista

print(norma([1,1],2))
print(normaliza([[1,1],[1,1],[1,1]],2))
