
import numpy as np
import matplotlib.pyplot as plt

nH = 2.32 # ZnS
nL = 1.38 # MgF2
ns = 1.5 # Índice de refracción del substrato
na = 1 # Aire, n4
N = 5 # Número de pares DBR
lo = 1500e-9 # Longitud de onda central


i = 0
r= []   # Almacenará los valores de reflectancia
        # para cada longitud de onda
             
l1 = [] # Almacenará los valores de cada longitud 
        # de onda a evaluar

for l in np.arange(250e-9, 2100e-9, 1e-9):
    r.append([])
    l1.append([]) 
    
    kH = (2*np.pi*nH)/l # Número de onda de la capa H
    dH = lo/(4*nH) # Grosor de la capa H
    pH = kH*dH # Propagación en el material H
    
    kL = (2*np.pi*nL)/l # Número de onda de la capa L
    dL = lo/(4*nL) # Grosor de la capa L
    pL = kL*dL # Propagación en el material L
    
    # Coeficientes matrices de transmisión DHL
    A1 = (nL+nH)/(2*nL)
    A2 = (nL-nH)/(2*nL)
    # Coeficientes matrices de transmisión DLH
    A3 = (nH+nL)/(2*nH)
    A4 = (nH-nL)/(2*nH)
    
    # Matrices de transmisión DHL y DLH
    DHL = np.array([[A1, A2], [A2, A1]])
    DLH = np.array([[A3, A4], [A4, A3]])
    
    # Matrices de propagación PH Y PL
    PH = np.array([[np.exp(1j*pH), 0], [ 0, np.exp(-1j*pH)]])
    PL = np.array([[np.exp(1j*pL), 0] , [0, np.exp(-1j*pL)]])

    # Obtención de la ecuación M = PL*DHL*PH*DLH
    dot1 = np.dot(PL,DHL)
    dot2 = np.dot(dot1,PH)
    M = np.dot(dot2,DLH)
    
    # Matriz M^N
    Mn = np.linalg.matrix_power(M,N)
    
    #nL y ns
    kL = (2*np.pi*nL)/l # Número de onda de la capa L
    dL = lo/(4*nL) # Grosor de la capa L
    pL = kL*dL # Propagación en el material L

    # Coeficientes de la matriz de transmisión DLs
    A3 = (ns+nL)/(2*ns)
    A4 = (ns-nL)/(2*ns)
    DLs = np.array([[A3,A4] , [A4, A3]])
    PL = np.array([[np.exp(1j*pL), 0] , [0, np.exp(-1j*pL)]])

    # n_H y na

    # Coeficientes de la matriz de transmisión DaH
    A3 = (nH+na)/(2*nH)
    A4 = (nH-na)/(2*nH)
    DaH = np.array([ [A3,A4],  [A4,A3]])
    
    # Obtención de la matriz de transferencia MT
    dot3 = np.dot(DLs,Mn)
    dot4 = np.dot(dot3,np.linalg.inv(DLH))
    MT = np.dot(dot4,DaH)

    ''' Aquí se almacenan valores de reflectancia r[i] y su 
                 correspondiente longitud de onda '''
  
    r[i].append(abs(MT[1,0]/MT[0,0]))
    l1[i].append(l)
    i = i+1

# Gráfico de reflectancia vs Longitud de onda

plt.plot(l1,r)
plt.title("Filtro infrarrojo")
plt.xlabel("Longitud de onda [m]")
plt.ylabel("Reflectancia")
plt.savefig("Filtro_IR_5_Capas.png")
plt.show()