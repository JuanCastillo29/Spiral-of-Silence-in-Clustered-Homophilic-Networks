# Spiral-of-Silence-in-Clustered-Homophilic-Network

## Abstract

Public discourse emerges from the interplay between individuals willingness to voice their opinions and the structural features of the social networks in which they are embedded. In this work, we investiate how choice homphily and triadic closure shape the emergence of the spiral of silence, the phenomenon whereby minority views are progressively silenced due to fear isolation. We advanced the state of the art in three ways. First, we integrate a realisitc network formation model, where homophily and triadic closure co-evolve, with a mean-field model of opinion expression. Second, we perform a bifurcation analysis of the associated Q-learning dynamics, revealing conditions for hysteresis cycles and path dependence in collective expression. Third, we validate our theoretical prediction through Monte Carlo simulatios, which highlight the role of fine-size effects and structural noise. Our results show that moderate triadic closure can foster minority expression by reinforcing local cohesion, whereas excesive closure amplifies assymetries and entreches majority dominance. These findings provide new insights into how algorithmic reinforcement of clustering in online platforms can either sustain diversity of opinion or accelerate its suppression.

---

## Code repository structure

The repository is divided into two different directories.

- 'C_codes'. Contains the computational expensive operations such as Monte Carlo simulations and integration algorithms. Two .c files are included, both ready to compile into different .dll libraries. The directory contains also a python script to compile the libraries (in Windows).
- 'Python_files'. This directory contains four python scripts used to simulate and integrate the system evolution to obtain the paper figures. Minor changes might be necessary to correctly reproduce the figures.


## Get started

This project includes C extensions that need to be compiled to generate dynamic libraries (.dll files). Follow the steps below to build and install the required extensions using Python.

### Requirements

- Python 3.x
- A compatible C compiler (e.g., Microsoft Visual C++ Build Tools)
- Windows (the build script uses the /O2 flag, which is specific to MSVC)

### Compiling the DLLs

1. Open a terminal in the root directory of the project.
2. Run the following command to build the extensions and generate optimized .dll files:
   ```bash python build_dlls.py ```

**Notes**: 
- The extensions are compiled with the /O2 flag to optimize for speed.
- If you're on Linux or macOS, you'll need to adapt the script (e.g., change compile flags and output file extensions like .so instead of .dll).

---

*Project by Juan Castillo y Emanuele Cozzo*

References:
This study is primarly based on:
-  A. Asikainen, G. Iñiguez, J. Ureña-Carrion, K. Kaski, and M. Kivelä, "Cumulative effects of triadic closure and homophily in social networks", Science Advances, vol 6, no 19, eaax7310, 2020.
-  F. Gaisbauer, E. Olbrich, and S. Banisch, "Dynamics of opinion expression", Physical Review E, vol 102, no 4, p. 042303, 2020
