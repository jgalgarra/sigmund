'''
Created on 31/05/2014

@author: algarra

    creainput lee una matriz de interacciones y crea las matrices de coeficientes
'''
import numpy as np
import os

def completamatrix(mat_nopesada,ncols,n0base,bijbase,alfabase,cibase,rdbase,nchunks=0):
    newfila=np.empty(ncols-1)
    mat_completa = mat_nopesada * bijbase
    # Multiplica los bij por un factor inverso a la distancia al core generalista
    if (nchunks>0):
        rangecols = np.array(range(0,ncols-1))
        chunkcols = np.array_split(rangecols,nchunks)
        for j in range(1,nchunks):
            for k in chunkcols[j]:
                mat_completa[:,k] = mat_completa[:,k]*((j+np.log(j))+1)
    # newfila=np.empty(ncols-1); newfila.fill(n0base)
    # Poblaciones iniciales
    mat_completa = np.vstack([mat_completa,n0base])
    # alfas
    newfila.fill(alfabase)
    mat_completa = np.vstack([mat_completa,newfila])
    # ci
    newfila.fill(cibase)
    mat_completa = np.vstack([mat_completa,newfila])
    # rb
    newfila.fill(0)
    mat_completa = np.vstack([mat_completa,newfila])
    # rd
    newfila.fill(rdbase)
    mat_completa = np.vstack([mat_completa,newfila])
    return(mat_completa)

def construyematrices(nomfich,bijbase = 0.000008,cibase = 0.000015,avgbase=0.05,rdbase = 0.04,randomizenest=0,randomizeall=0,nchunks=0):
    my_data = np.genfromtxt(nomfich, delimiter=',')
    nfilas, ncols = my_data.shape
    # Creamos un array vacio con una linea menos y una columna menos para eliminar
    # los nombres de las especies
    mat_inter = np.zeros((nfilas-1,ncols-1))
    for i in range(1,nfilas):
        mat_inter[i-1,] = my_data[i,1:ncols]
    if (randomizenest):
        for j in range(0,ncols-1):
            np.random.shuffle(mat_inter[:,j])
        if (randomizeall):
            for k in range (0,nfilas-1):
                np.random.shuffle(mat_inter[k])
    #print(mat_inter.shape)
    alfabase = cibase*10
    sustrato = avgbase/cibase
    n0base = np.round(0.2*sustrato+np.random.normal(sustrato,0.3*sustrato,ncols-1))
    mat_A = completamatrix(mat_inter,ncols,n0base,bijbase,alfabase,cibase,rdbase,nchunks)
    n0base = np.round(0.2*sustrato+np.random.normal(sustrato,0.3*sustrato,nfilas-1)) 
    mat_B = completamatrix(mat_inter.transpose(),nfilas,n0base,bijbase,alfabase,cibase,rdbase,nchunks)
    return(mat_A,mat_B)

nfich='M_SD_004.csv'
mat_A, mat_B = construyematrices(nfich,bijbase = 0.0000011,avgbase=0.025,rdbase = 0.045,cibase = 0.000028,randomizenest=0,randomizeall=0,nchunks=80)
prefix = nfich.split('.')[0]
#print(prefix)
np.savetxt(prefix+'_a.txt',mat_A,fmt='%0.08f',delimiter='\t',newline=os.linesep)
np.savetxt(prefix+'_b.txt',mat_B,fmt='%0.08f',delimiter='\t',newline=os.linesep)