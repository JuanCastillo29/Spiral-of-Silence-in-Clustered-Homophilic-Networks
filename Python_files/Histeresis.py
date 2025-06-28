# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 13:59:28 2025

@author: juanc
"""
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import os


#%% ========= CARGA DE DLLs ==========

# Campo medio
lib = ctypes.CDLL(os.path.abspath("bias_solver.dll"))
lib.IntegracionEcBias.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
]
lib.QLearningdynamics.argtypes = [
    ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_double, ctypes.c_double
]
lib.IntegracionEcBias.restype = None
lib.QLearningdynamics.restype = None



# Red compleja
class Nodo(ctypes.Structure):
    _fields_ = [('vecinos', ctypes.POINTER(ctypes.c_int)),
                ('Q', ctypes.c_double),
                ('NVecinos', ctypes.c_int),
                ('opinion', ctypes.c_int),
                ('prob', ctypes.c_double)]

class Lattice(ctypes.Structure):
    _fields_ = [('NNodos', ctypes.c_int),
                ('NNodos0', ctypes.c_int),
                ('nodos', ctypes.POINTER(Nodo)),
                ('beta', ctypes.c_double)]

lib_net = ctypes.CDLL('./Lib_red.dll.dll')
lib_net.ini_ran.argtypes = [ctypes.c_long]
lib_net.Simulacion.argtypes = [ctypes.c_int, ctypes.POINTER(Lattice), ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double)]
lib_net.Simulacion.restype = None
lib_net.TriadicClosure.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                   ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                   ctypes.c_double, ctypes.c_double]
lib_net.TriadicClosure.restype = ctypes.POINTER(Lattice)
lib_net.ExtraerComponenteMayor.argtypes = [ctypes.POINTER(Lattice)]
lib_net.ExtraerComponenteMayor.restype = ctypes.POINTER(Lattice)
lib_net.FreeLatice.argtypes = [ctypes.POINTER(Lattice)]
lib_net.FreeLatice.restype = None
lib_net.Upgrading.argtypes = [ctypes.POINTER(Lattice), ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
lib_net.Upgrading.restype = ctypes.POINTER(Lattice)
#%% ==== Funciones generales ====
def simular_histeresis(beta, c_directo, c_inverso):
    Qa = ctypes.c_double(0.8)
    Qb = ctypes.c_double(0.8)
    Ta = ctypes.c_double(n)
    Tb = ctypes.c_double(1 - n)

    QA_d = np.zeros_like(c_directo)
    QB_d = np.zeros_like(c_directo)
    QA_i = np.zeros_like(c_inverso)
    QB_i = np.zeros_like(c_inverso)

    # Recorrido directo
    for i, c in enumerate(c_directo):
        lib.IntegracionEcBias(c, sa, sb, n, ctypes.byref(Ta), ctypes.byref(Tb))
        lib.QLearningdynamics(Ta, Tb, ctypes.byref(Qa), ctypes.byref(Qb), beta, coste)
        QA_d[i] = Qa.value
        QB_d[i] = Qb.value

    # Recorrido inverso
    for i, c in enumerate(c_inverso):
        lib.IntegracionEcBias(c, sa, sb, n, ctypes.byref(Ta), ctypes.byref(Tb))
        lib.QLearningdynamics(Ta, Tb, ctypes.byref(Qa), ctypes.byref(Qb), beta, coste)
        QA_i[i] = Qa.value
        QB_i[i] = Qb.value

    return QA_d, QB_d, QA_i, QB_i

def graficar(c_directo, c_inverso, resultados, componente, nombre_fichero):
    colores = {10.0: 'blue', 7.5: 'red', 5.0: 'green'}
    
    plt.figure(figsize=(8, 6))
    ax = plt.axes()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlim(0.3, 0.9)
    plt.ylim(-2.0, 1.0)

    for beta in betas:
        plt.plot(c_directo, resultados[beta][f"{componente}_directo"],
                 label=fr'$ \beta = {beta} $', color=colores[beta])
        plt.plot(c_inverso, resultados[beta][f"{componente}_inverso"],
                 color=colores[beta])

    plt.xlabel("c", fontweight='bold', fontsize=20)
    plt.ylabel(f"Q$_{1 if componente == 'QA' else '2'}$", fontweight='bold', fontsize=20)
    plt.legend(frameon=True, facecolor='white', edgecolor='black', fontsize=12, loc='best')
    plt.grid(False)
    plt.savefig(nombre_fichero, format='pdf', bbox_inches='tight')
    plt.show()


def ciclo_histeresis_red(beta, c_vals, sim_id, coste=0.2, sa=0.6525, sb=0.65, TTotal=5000,
                         nnodos=100000, n0=0.5, k=20 ):
  
    lib_net.ini_ran(987654321 + sim_id)

    c_hist = list(c_vals) + list(c_vals[::-1])  # ida + vuelta
    lattice_ptr = lib_net.TriadicClosure(5001, nnodos, n0, k, c_hist[0], sa, sb, 1 - sa, 1 - sb)
    Q1_vals = []
    Q2_vals = []
    Q_local = np.random.rand(2).astype(np.float64)
    for c in c_hist:
        for i in range(3001):
            lattice_ptr = lib_net.Upgrading(lattice_ptr, sa, sb, 1-sa, 1-sb, c)
        lattice_ptr2 = lib_net.ExtraerComponenteMayor(lattice_ptr)
        if not lattice_ptr2:
            Q1_vals.append(np.nan)
            Q2_vals.append(np.nan)
            continue

        
        lib_net.Simulacion(TTotal, lattice_ptr2, coste, beta, Q_local.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        lib_net.FreeLatice(lattice_ptr2)
        
        Q1_vals.append(Q_local[0])
        Q2_vals.append(Q_local[1])
        
    lib_net.FreeLatice(lattice_ptr)
    return np.array(c_hist), np.array(Q1_vals), np.array(Q2_vals)

# Tu función individual de simulación
def simulacion_histeresis(args):
    beta, c_vals, coste, sa, sb, TTotal, n0, sim_id = args
    c_hist, Q1, Q2 = ciclo_histeresis_red(beta, c_vals, sim_id, coste=coste, sa=sa, sb=sb, TTotal=TTotal, n0=n0)
    return Q1, Q2


#%% ==== Parámetros generales ==== 
n = 0.2
sa = 0.755
sb = 0.65
coste = 0.2
betas = [10.0, 7.5, 5.0]
c_directo = np.linspace(0.9, 0.7, 71)
c_inverso = np.linspace(0.7, 0.9, 71)


if __name__ == "__main__":
    modo = 'histeresis_red'  # Nuevo modo para hacer ciclo en la red

    if modo == 'histeresis_red':

        beta = 10
        resultados_dict = {}  # Diccionario para guardar resultados por beta
        c_vals = np.linspace(0.9, 0.7, 21)
        N_runs = 5

        # Empaquetar argumentos para cada simulación
        args_list = [(beta, c_vals, coste, sa, sb, 5000, n, sim_id) for sim_id in range(N_runs)]
        c_hist = np.concatenate([c_vals, c_vals[::-1]])

        # Multiprocessing con barra de progreso
        with Pool(processes=cpu_count() - 1) as pool:
            resultados_lista = list(tqdm(pool.imap(simulacion_histeresis, args_list), total=N_runs, desc='Simulando'))

        # Acumular resultados
        Q1_total = np.array([res[0] for res in resultados_lista])
        Q2_total = np.array([res[1] for res in resultados_lista])
        Q1_avg = np.mean(Q1_total, axis=0)
        Q2_avg = np.mean(Q2_total, axis=0)
        Q1_std = np.std(Q1_total, axis=0)
        Q2_std = np.std(Q2_total, axis=0)

        # Guardar resultados en CSV
        df = pd.DataFrame({
            'c': c_hist,
            'Q1_avg': Q1_avg,
            'Q1_std': Q1_std,
            'Q2_avg': Q2_avg,
            'Q2_std': Q2_std
        })
        df.to_csv(f'Ciclo_Histeresis_beta_{beta}_avg.csv', index=False)

        # Ejecutar simulación directa e inversa
        QA_d, QB_d, QA_i, QB_i = simular_histeresis(beta, c_directo, c_inverso)
        resultados_dict[beta] = {
            "QA_directo": QA_d,
            "QB_directo": QB_d,
            "QA_inverso": QA_i,
            "QB_inverso": QB_i
            }

        # Gráfico
        plt.figure(figsize=(8, 6))
        ax = plt.axes()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.plot(c_directo, resultados_dict[beta]["QA_directo"], color='blue',  alpha=0.4)
        plt.plot(c_inverso, resultados_dict[beta]["QA_inverso"], color='blue',  alpha=0.4)
        plt.plot(c_directo, resultados_dict[beta]["QB_directo"], color='red', alpha=0.4)
        plt.plot(c_inverso, resultados_dict[beta]["QB_inverso"], color='red', alpha=0.4)

        plt.errorbar(c_hist, Q1_avg, yerr=Q1_std, fmt='o', markersize=4, capsize=3,
             label='Q$_{1}$', color='blue', ecolor='blue', alpha=0.8)

        plt.errorbar(c_hist, Q2_avg, yerr=Q2_std, fmt='s', markersize=4, capsize=3,
             label='Q$_{2}$', color='red', ecolor='red', alpha=0.8)

        plt.xlabel('Triadic closure (c)', fontweight='bold', fontsize=20)
        plt.ylabel('Q$_{i}$', fontweight='bold', fontsize=20)
        plt.legend(frameon=True, facecolor='white', edgecolor='black', fontsize=12, loc='best')
        plt.tight_layout()
        plt.savefig(f'Ciclo_Histeresis_Q1_beta_{beta}_avg.pdf')
        plt.show()


        
    if modo == 'histeresis_campo_medio':
        resultados = {}
        for beta in betas:
            QA_d, QB_d, QA_i, QB_i = simular_histeresis(beta, c_directo, c_inverso)
            resultados[beta] = {
                "QA_directo": QA_d,
                "QB_directo": QB_d,
                "QA_inverso": QA_i,
                "QB_inverso": QB_i
            }


        # ==== Graficar resultados ====
        graficar(c_directo, c_inverso, resultados, componente='QB', nombre_fichero='Histeresis_Q2.pdf')
        graficar(c_directo, c_inverso, resultados, componente='QA', nombre_fichero='Histeresis_Q1.pdf')
