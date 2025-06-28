
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 17:47:45 2025

@author: juanc
"""
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
# Defining C structs in Python
class Nodo(ctypes.Structure):
    _fields_ = [
        ('vecinos', ctypes.POINTER(ctypes.c_int)),
        ('Q', ctypes.c_double),
        ('NVecinos', ctypes.c_int),
        ('opinion', ctypes.c_int),
        ('prob', ctypes.c_double),
    ]

class Lattice(ctypes.Structure):
    _fields_ = [
        ('NNodos', ctypes.c_int),
        ('NNodos0', ctypes.c_int),
        ('nodos', ctypes.POINTER(Nodo)),
        ('beta', ctypes.c_double),
    ]

# Load C lib.
lib = ctypes.CDLL('./Lib_red.dll') 
lib.ini_ran.argtypes = [ctypes.c_long]
lib.Simulacion.argtypes = [ctypes.c_int, ctypes.POINTER(Lattice), ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double) ]
lib.Simulacion.restype = None

lib.TriadicClosure.argtypes = [ ctypes.c_int,
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double
]
lib.TriadicClosure.restype = ctypes.POINTER(Lattice)

lib.FreeLatice.argtypes = [ctypes.POINTER(Lattice)]
lib.FreeLatice.restype = None

lib.ExtraerComponenteMayor.argtypes = [ctypes.POINTER(Lattice)]
lib.ExtraerComponenteMayor.restype = ctypes.POINTER(Lattice)

#%%
#Simulación de Q-learning



#resultados_por_lambda = {}



def correr_simulacion(args):
    #Parametros del modelo de expresion

    coste = 0.2

    #Parametros del modelo de formacion de red
    sa = 0.6525 #Choice homophily de la minoría
    sb = 0.65 #Choice homophily de la mayoría
    c= 0.785

    TTotal = 5000
    nnodos = 10000
    n0 = 0.5
    k = 20
    
    beta, sim_id = args
    lib.ini_ran(123456789+sim_id)
    # Crear red
    lattice_ptr = lib.TriadicClosure(3000, nnodos, n0, k, c, sa, sb, 1-sa, 1-sb)
    lattice_ptr = lib.ExtraerComponenteMayor(lattice_ptr)
    if not lattice_ptr:
        return None

    # Crear buffer Q local
    Q_local = np.random.rand(2).astype(np.float64)
    Q_ptr = Q_local.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    lib.Simulacion(TTotal, lattice_ptr, coste, beta, Q_ptr)

    lib.FreeLatice(lattice_ptr)
    return (beta, Q_local[0], Q_local[1])

if __name__ == "__main__":
    resultados_por_lambda = {}
    lista_de_beta = np.linspace(0, 30, 31)
    # Armar lista de tareas
    tareas = []
    num_simulaciones = 100
    for beta in lista_de_beta:
        for sim_id in range(num_simulaciones):
            tareas.append((beta, sim_id))

    # Ejecutar con barra de progreso
    with mp.Pool(processes=7) as pool:
        resultados_brutos = list(tqdm(pool.imap(correr_simulacion, tareas), total=len(tareas)))


    # Agrupar por c
    for beta in lista_de_beta:
        filtrados = [ (r[1], r[2]) for r in resultados_brutos if r and r[0] == beta ]
        resultados_por_lambda[beta] = np.array(filtrados)


    
    lambda_vals_sorted = sorted(resultados_por_lambda.keys())

    print("Simulacion terminada.")
    # Crear lista de filas
    rows = []

    for lambda_val in lambda_vals_sorted:
        Q_vals = resultados_por_lambda[lambda_val]
        for q1, q2 in Q_vals:
            rows.append({'beta': lambda_val, 'Q1': q1, 'Q2': q2})

# Crear DataFrame
    df_Q = pd.DataFrame(rows)
    df_Q.to_csv('datos_Q (beta, sa =0.6525, sb=0.65, c=0.785).csv', index=False)

    plt.figure(figsize=(10, 6))
 # # ---------- Dispersión de Q2 ----------
    for lambda_val in lambda_vals_sorted:
          Q_vals = resultados_por_lambda[lambda_val]
          c_vals = [lambda_val] * len(Q_vals)  # Eje x: mismo valor para todos los Q2
          plt.scatter(c_vals, Q_vals[:, 1], alpha=0.3, color='red', label='_nolegend_')  # '_nolegend_' evita repetición en leyenda
   

# # ---------- Etiquetas y leyenda ----------
    plt.xlabel('T$^{-1}$', fontweight='bold', fontsize=20)
    plt.ylabel('Q$_2$', fontweight='bold', fontsize=20)
    plt.xlim(min(lista_de_beta) - 1, max(lista_de_beta) + 1)
    plt.ylim(-1.01, 1.01)

    plt.grid(False)
    # plt.legend(
    #     frameon=True,
    #     facecolor='white',
    #     edgecolor='black',
    #     fontsize=12,
    #     loc='best',
    #     )

    plt.tight_layout()
    plt.savefig('Campo_Medio_Q2_sobre_curvas.pdf', format='pdf', bbox_inches='tight')
    plt.show()

    plt.figure(figsize=(10, 6))

#     # # ---------- Dispersión de Q1 ----------
    for lambda_val in lambda_vals_sorted:
        Q_vals = resultados_por_lambda[lambda_val]
        c_vals = [lambda_val] * len(Q_vals)  # Eje x: lambda usado como c
        plt.scatter(c_vals, Q_vals[:, 0], alpha=0.3, color='blue', label='_nolegend_')  # '_nolegend_' evita duplicar en leyenda
        
    #plt.scatter(df['c'], df['Q1'], color = 'red', label='_nolegend_')

 # ---------- Etiquetas y leyenda ----------
    plt.xlabel('T$^{-1}$', fontweight='bold', fontsize=20)
    plt.ylabel('Q$_1$', fontweight='bold', fontsize=20)
    
    plt.xlim(min(lista_de_beta) - 1, max(lista_de_beta) + 1)
    plt.ylim(-1.01, 1.01)
    
    plt.grid(False)
    # plt.legend(
    #     frameon=True,
    #     facecolor='white',
    #     edgecolor='black',
    #     fontsize=12,
    #     loc='best',
    #     )

    plt.tight_layout()
    plt.savefig('Campo_Medio_Q1_sobre_curvas.pdf', format='pdf', bbox_inches='tight')
    plt.show()
 

