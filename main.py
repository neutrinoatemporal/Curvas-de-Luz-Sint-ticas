import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.optimize import fsolve
from sympy import *
from sympy.plotting import plot

from Funciones import newton_raphson, areaseg, area1, area2, ra, phiuno, phidos

#Función para graficar la curva de luz sintética // Parametros: (ángulo de inclinación, razon de masas,
# radio relativo de estrella 1, radio relativo de estrella 2 , razon de luminosidades)

def curva( i, q, r1, r2,l):
    l1=1
    l2=l
    R=1
    fase= np.linspace(0, 1, 600) # 600 puntos
    theta = 2 * np.pi * fase
    y_values = []

    for th in theta:
        x = R * np.sin(th)
        y = R * np.cos(i) * np.cos(th)
        z = R * np.sin(i) * np.cos(th)

        x1 = -x / (1 + (1 / q))
        y1 = -y / (1 + (1 / q))
        z1 = -z / (1 + (1 / q))

        x2 = x / (1 + q)
        y2 = y / (1 + q)
        z2 = z / (1 + q)

        ro = ra(x1, x2, y1, y2)
        phi1 = phiuno(r1, r2, ro)
        phi2 = phidos(r1, r2, ro)
        a1 = area1(ro, r1, r2, z1, z2, phi1, phi2)
        a2 = area2(ro, r1, r2, z1, z2, phi1, phi2)

        y_val = (l1 * a1 / r1**2 + l2 * a2 / r2**2)  #Valores de flujo relativo
        y_values.append(y_val)

    y_values = np.array(y_values)

    index_025 = np.argmin(np.abs(fase - 0.25))  # Encontrar el índice más cercano a fase=0.25
    index_075 = np.argmin(np.abs(fase - 0.75))  # Encontrar el índice más cercano a fase=0.75
    Lmax = (y_values[index_025] + y_values[index_075]) / 2  # Promedio
    # Normalizar flujo dividiendo todos los valores por Lmax (flujo máximo)
    y_values = y_values / Lmax

    plt.plot(fase, y_values, 'o', markersize=1.5, color='black') #Graficar
    plt.xlabel(r'Fase')
    plt.ylabel(r'Flujo Relativo')
    plt.xlim(0, 1)
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.show()