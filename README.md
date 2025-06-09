# Microcomb Simulation – LLE Solver

**Caitlin E. Murray, Chawaphon Prayoonyong, and Bill Corcoran, 2025**  
*Written in MATLAB version 2022b*

A split-step symmetrical Fourier solver for the Lugiato–Lefever Equation (LLE), designed to simulate microcomb formation in Kerr microresonators. This simple, well-commented version is intended as a starting point for students or researches unfamiliar with the LLE.

The current setup uses a constant-power continuous-wave (CW) input and generates solitons by varying the frequency detuning between the pump laser and the microring resonator.

---

## Files

Currently, this simulation contains four files:

- **`LLE_simple_commented.m`** – The main simulation file where data is generated and saved.  
- **`parsave.m`** – A saving function that can be used within a `parfor` loop, called from the main simulation.  
- **`SaveClass.m`** – A simple class definition used for organizing the saved data.

Both `parsave.m` and `SaveClass.m` are currently required to run `LLE_simple_commented.m`, though the code can be easily modified to remove this dependency if desired.

A simple script for plotting the output from a single run, **`Read_LLE_data.m`**, has also been included.



## Conventions

Variable names follow the conventions used in:  

**Pasquazi et al.,** *"Micro-combs: A novel generation of optical sources"* (2018)

Which is a fantastic resource.


## Split-Step Fourier Solver

The symmetrical split-step method was chosen for its simplicity and robustness. It is fairly robust to under and overflow, particulary when dealing with unusual initial conditions or feedback loops. This method is straightforward to troubleshoot and relatively fast to run, making it a practical choice, of course more accurate solving methods do exist, such as the Runge Kutta methods.



## MATLAB Implementation

This solver was written in MATLAB R2022b. MATLAB was chosen over Python as it currently computes FFTs faster, and with the heavy reliance on FFTs, a single run is faster in MATLAB.
The code is compatible with the Parallel Computing Toolbox for parallel execution using parfor. Python remains advantageous for large batch processing and more scalable parallelism.


## Dispersion
Includes third-order dispersion (D3) and optional avoided mode crossings (AMX). Can be easily extended to include higher-order terms or loaded from measured dispersion data. Dispersion is formated in a frequency shift per mode basis, which can be converted from or into the material dispersion. See **Kovach et al.,** *"Emerging material systems for integrated optical Kerr frequency combs"* (2020), for conversion equations. 

## Not Included

This version does **not** include:

- Raman scattering
- Shot noise
- Thermal drift
---

## Notes
- AI was used to spell-check comments, and minor formatting across the code and repository.


## Other useful References

- **Kovach et al.,** *"Emerging material systems for integrated optical Kerr frequency combs"* (2020)  
- **Kippenberg et al.,** *"Dissipative Kerr solitons in optical microresonators"* (2018)
