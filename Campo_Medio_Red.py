# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 17:04:47 2025

@author: juanc
"""
import os
import numpy as np
import pandas as pd
from itertools import product
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap
import ctypes

# Cargar la DLL
lib = ctypes.CDLL(os.path.abspath("bias_solver.dll"))

# Configurar argumentos y tipo de retorno
lib.IntegracionEcBias.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
]
lib.IntegracionEcBias.restype = None

def integracion_ec_bias(c, sa, sb, n):
    Ta = ctypes.c_double(n)
    Tb = ctypes.c_double(1-n)
    lib.IntegracionEcBias(c, sa, sb, n, ctypes.byref(Ta), ctypes.byref(Tb))
    return Ta.value, Tb.value


# Función que ejecuta una simulación para un conjunto de parámetros
def worker(params):
    c, sa, sb, n = params
    Ta, Tb = integracion_ec_bias(c, sa, sb, n)
    return {'c': c, 'sa': sa, 'sb': sb, 'n': n, 'Ta': Ta, 'Tb': Tb}

# Redondea una tupla de floats para comparaciones robustas
def rounded_tuple(p, decimals=8):
    return tuple(round(x, decimals) for x in p)

# Carga los resultados existentes desde archivo, si existe
def load_existing_results(filename):
    if os.path.exists(filename):
        df = pd.read_csv(filename)
        existing_keys = set(rounded_tuple(p) for p in zip(df['c'], df['sa'], df['sb'], df['n']))
        return df, existing_keys
    else:
        return pd.DataFrame(), set()

def PhaseSpace(data,  c_values, s_values, save_path="phase_space_plot.pdf"):
  C, S = np.meshgrid(c_values, s_values)

  # Initialize color matrix with NaNs
  color_matrix = np.full(C.shape, None, dtype=object)  # Use object dtype to store strings

  # Assign string colors from DataFrame
  for index, row in data.iterrows():
    c_idx = np.argmin(np.abs(c_values - row['c']))
    s_idx = np.argmin(np.abs(s_values - row['sa']))
    color_matrix[s_idx, c_idx] = row['color']  # Store the color name

  # Convert unique colors into a colormap
  unique_colors = data['color'].unique()  # Extract unique color names
  cmap = ListedColormap(unique_colors)  # Create colormap

  # Convert the color matrix to indices in the colormap
  color_indices = np.vectorize(lambda x: np.where(unique_colors == x)[0][0] if x in unique_colors else -1)(color_matrix)

  # Plot using the colormap
  plt.figure(figsize=(6, 5))
  plt.pcolormesh(C, S, color_indices, cmap=cmap, shading='auto')
  plt.xlabel("Triadic closure (c)", fontsize=12, fontweight='bold')
  plt.ylabel("Choice homophily G$_{1}$ (s$_{1}$)", fontsize=12, fontweight='bold')

  plt.savefig(save_path, bbox_inches='tight', dpi=300)

  plt.show()
  
def get_color(T_m, T_M, c):
    if T_m > (1+c)/2 and T_M > (1+c)/2:
        return 'green'
    if T_m < c and T_M < c:
        return 'purple'
    if T_m < c:
        return 'red'
    if T_M < c:
        return 'blue'
    return 'white'  # Default color for other cases

if __name__ == '__main__':
    output_file = "n=0.2, Bias.csv"

    num_s = 100
    num_c = 99
    coste = 0.2
    
    s_edges = np.linspace(0.5, 1.0, num_s + 1)
    c_edges = np.linspace(0, 0.99, num_c + 1)

    s_values = (s_edges[:-1] + s_edges[1:]) / 2
    c_values = (c_edges[:-1] + c_edges[1:]) / 2
    n = 0.2
    sM = 0.65
    
    # Solo combinaciones donde sa == sb
    param_grid = [(c, s, sM, n) for c in c_values for s in s_values]

    # Cargar resultados previos
    existing_df, existing_keys = load_existing_results(output_file)

    # Filtrar combinaciones nuevas
    new_grid = [p for p in param_grid] #if rounded_tuple(p) not in existing_keys]

    if new_grid:
        with Pool(cpu_count()-1) as pool:
            new_results = list(tqdm(pool.imap_unordered(worker, new_grid), total=len(new_grid)))

        new_df = pd.DataFrame(new_results)
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        #combined_df.to_csv(output_file, index=False)
        print(f"✅ Resultados actualizados en: {output_file}")
    else:
        print("✅ No hay nuevas combinaciones para calcular.")
    
    new_df['color'] = new_df.apply(lambda row: get_color(row['Ta'], row['Tb'], coste), axis=1)
    
    PhaseSpace(new_df, c_values, s_values)