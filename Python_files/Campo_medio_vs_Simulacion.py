# -*- coding: utf-8 -*-
"""
Created on Sat Jun 14 17:44:10 2025

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
lib = ctypes.CDLL('./Lib_red.dll.dll')  # Reemplaza con el nombre de tu biblioteca compilada
lib2 = ctypes.CDLL('./bias_solver.dll') 
lib.ini_ran.argtypes = [ctypes.c_long]
lib.Simulacion.argtypes = [ctypes.c_int, ctypes.POINTER(Lattice), ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double) ]
lib.Simulacion.restype = None
# Definir firmas de funciones
lib2.IntegracionEcBias.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
]

lib2.IntegracionEcBias.restype = None
lib.TriadicClosure.argtypes = [ ctypes.c_int,
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double
]
lib.TriadicClosure.restype = ctypes.POINTER(Lattice)

lib.FreeLatice.argtypes = [ctypes.POINTER(Lattice)]
lib.FreeLatice.restype = None

lib.ExtraerComponenteMayor.argtypes = [ctypes.POINTER(Lattice)]
lib.ExtraerComponenteMayor.restype = ctypes.POINTER(Lattice)

# Función logística (sigmoide)
def sigmoid(x, beta):
    return 1 / (1 + np.exp(-beta * x))

def Plot(T1, T2, beta, xi):
    df = pd.read_csv('datos_Q.csv') 
    T11 = T1.value
    T22 = T2.value
    q1 = np.linspace(-1, 1, 20)
    q2 = np.linspace(-1, 1, 20)
    Q1, Q2 = np.meshgrid(q1, q2)

    U = T11 * sigmoid(Q1, beta) - (1.0 - T11) * sigmoid(Q2, beta) - Q1 - xi
    V = T22 * sigmoid(Q2, beta) - (1.0 - T22) * sigmoid(Q1, beta) - Q2 - xi

    plt.figure(figsize=(5, 5))
    plt.grid(False)
    plt.yticks(np.arange(-1, 1.01, 0.5))
    plt.xticks(np.arange(-1, 1.01, 0.5))
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)

    # Colores para cada tamaño
    colores = {10**3: 'red', 10**4: 'green', 10**5: 'orange'}
    etiquetas = {10**3: 'N = $10^3$', 10**4: 'N = $10^4$', 10**5: 'N = $10^5$'}

    for N, grupo in df.groupby('N'):
        plt.scatter(grupo['Q1'], grupo['Q2'], color=colores[N], label=etiquetas[N], alpha=0.7)

    plt.quiver(Q1, Q2, U, V, color='blue', alpha=0.5)
    plt.xlabel('Q$_{1}$', fontweight='bold', fontsize=20)
    plt.ylabel('Q$_{2}$', fontweight='bold', fontsize=20)
    plt.legend(
         frameon=True,
         facecolor='white',
         edgecolor='black',
         fontsize=12,
         loc='best',
         )
    plt.savefig('Campo_medioVSsimulacion.pdf', format='pdf', bbox_inches='tight')
    plt.show()


def correr_simulacion(args):
    
    c, k, sa, sb, nnodos, n0, TTotal, coste, beta, sim_id = args
    lib.ini_ran(123456789 + nnodos*sim_id)  # para que cada sim tenga una semilla distinta

    # Crear red
    lattice_ptr = lib.TriadicClosure(5000, nnodos, n0, k, c, sa, sb, 1-sa, 1-sb)
    lattice_ptr = lib.ExtraerComponenteMayor(lattice_ptr)
    if not lattice_ptr:
        return None

    # Crear buffer Q local
    Q_local = np.random.rand(2).astype(np.float64)
    Q_ptr = Q_local.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    lib.Simulacion(TTotal, lattice_ptr, coste, beta, Q_ptr)

    lib.FreeLatice(lattice_ptr)
    return (c, Q_local[0], Q_local[1])

def sim_wrapper(args):
    res = correr_simulacion(args)
    if res is not None:
        return [res[0], res[1], res[2], args[4]]  # args[4] es nnodos (N)


if __name__ == "__main__":
    mp.freeze_support()
    
    # Parámetros fijos
    c = 0.785
    sa = 0.6525
    sb =0.65
    k = 25
    n0 = 0.5
    coste = 0.2
    beta = 25
    TTotal = 3000

    tamaños = [10**3, 10**4, 10**5]
    resultados_totales = []

    for N in tamaños:
        print(f"Simulando para N = {N}")

        if N == 10**5:
            n_sim=5
            for red_id in range(10):  # 10 redes distintas
                print(f"Red {red_id + 1}/10 para N = {N}")
        
                # Semilla distinta para cada red
                lib.ini_ran(987654321 + red_id)
                lattice_ptr = lib.TriadicClosure(5000, N, n0, k, c, sa, sb, 1-sa, 1-sb)
                lattice_ptr = lib.ExtraerComponenteMayor(lattice_ptr)
                if ctypes.cast(lattice_ptr, ctypes.c_void_p).value is None:
                    print(f"No se pudo construir la red #{red_id} para N = {N}")
                    continue

                for i in tqdm(range(n_sim), desc=f"Simulaciones en red {red_id + 1}"):
                    Q_local = np.random.rand(2).astype(np.float64)
                    Q_ptr = Q_local.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    lib.ini_ran(123456789 + red_id * 100 + i)  # Semilla distinta por sim
                    lib.Simulacion(TTotal, lattice_ptr, coste, beta, Q_ptr)
                    resultados_totales.append([c, Q_local[0], Q_local[1], N])

                lib.FreeLatice(lattice_ptr)


        else:
            n_sim = 100
            # Red nueva por cada simulación
            args_list = [(c, k, sa, sb, N, n0, TTotal, coste, beta, i) for i in range(n_sim)]

          
            with mp.Pool(processes=6) as pool:
                for res in tqdm(pool.imap(sim_wrapper, args_list), total=n_sim, desc=f"N = {N}"):
                    if res:
                        resultados_totales.append(res)

    # Guardar resultados
    df_resultados = pd.DataFrame(resultados_totales, columns=["c", "Q1", "Q2", "N"])
    df_resultados.to_csv('datos_Q.csv', index=False)
    Ta = ctypes.c_double(n0)
    Tb = ctypes.c_double(1 - n0)
    lib2.IntegracionEcBias(c, sa, sb, n0, ctypes.byref(Ta), ctypes.byref(Tb))
    # Plot
    Plot(Ta, Tb, beta, coste)
    print(f"Simulaciones exitosas: {len(df_resultados)}")



