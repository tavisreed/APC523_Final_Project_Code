import Helper_Functions as hf
import numpy as np

# Note, These (and the Run_X_40_V2 files) are loaded differnetly because there were run using an older method of the save function which did not
# keep the tau list

# Load in the the concentration matrixes
run_0_40 =hf.load_object('./Run_0_40_V2')
run_24_40_inv =hf.load_object('./Run_24_40_V2_inv')
run_48_40_inv =hf.load_object('./Run_48_40_V2_inv')
run_72_40_inv=hf.load_object('./Run_72_40_V2_inv')
run_96_40_inv =hf.load_object('./Run_96_40_V2_inv')


# Get the steady state concentrations
Pyr = [run_0_40[0,-1], run_24_40_inv [0,-1], run_48_40_inv [0,-1], run_72_40_inv[0,-1], run_96_40_inv [0,-1]]
NAD = [run_0_40[1,-1], run_24_40_inv [1,-1], run_48_40_inv [1,-1], run_72_40_inv[1,-1], run_96_40_inv [1,-1]]
OAA = [run_0_40[2,-1], run_24_40_inv [2,-1], run_48_40_inv [2,-1], run_72_40_inv[2,-1], run_96_40_inv [2,-1]]
AcCoA = [run_0_40[3,-1], run_24_40_inv [3,-1], run_48_40_inv [3,-1], run_72_40_inv[3,-1], run_96_40_inv [3,-1]]
Cit = [run_0_40[4,-1], run_24_40_inv [4,-1], run_48_40_inv [4,-1], run_72_40_inv[4,-1], run_96_40_inv [4,-1]]
KG = [run_0_40[5,-1], run_24_40_inv [5,-1], run_48_40_inv [5,-1], run_72_40_inv [5,-1], run_96_40_inv [5,-1]]
ATP = [run_0_40[6,-1], run_24_40_inv [6,-1], run_48_40_inv [6,-1], run_72_40_inv [6,-1], run_96_40_inv [6,-1]]
psi = [run_0_40[7,-1], run_24_40_inv [7,-1], run_48_40_inv [7,-1], run_72_40_inv [7,-1], run_96_40_inv [7,-1]]

# Plot one of the concentrations curves for different HPIs
hf.plot_hpi_concentration(Cit, title='Citrate')

# Calculate tau
tau = [0]
for i in range(np.shape(run_24_40_inv)[1]-1):
    tau.append(tau[-1]+38/0.140*0.0000001)

# Plot the Concentration matrix (Note: It make take some time to plot, as there are a lot of points)
hf.plot_concentration_matrix(run_24_40_inv, tau, title='Concentration Curves, 24 HPI')


# Note: This is how to do the same things as above, but for runs saved using the new save function
run_0_40 = hf.load_object('./Run_0_40_V2')
run_24_40_inv_steady = hf.load_object('./Run_24_40_V2_inv_steady')['concentration_matrix']
run_24_40_inv_steady_tau = hf.load_object('./Run_24_40_V2_inv_steady')['tau']
run_48_40_inv_steady = hf.load_object('./Run_48_40_V2_inv_steady')['concentration_matrix']
run_48_40_inv_steady_tau = hf.load_object('./Run_48_40_V2_inv_steady')['tau']
run_72_40_inv_steady = hf.load_object('./Run_72_40_V2_inv_steady')['concentration_matrix']
run_72_40_inv_steady_tau = hf.load_object('./Run_72_40_V2_inv_steady')['tau']
run_96_40_inv_steady = hf.load_object('./Run_96_40_V2_inv_steady')['concentration_matrix']
run_96_40_inv_steady_tau = hf.load_object('./Run_96_40_V2_inv_steady')['tau']

# Get the steady state concentrations
Pyr = [run_0_40[0,-1], run_24_40_inv_steady [0,-1], run_48_40_inv_steady [0,-1], run_72_40_inv_steady[0,-1], run_96_40_inv_steady [0,-1]]
NAD = [run_0_40[1,-1], run_24_40_inv_steady [1,-1], run_48_40_inv_steady [1,-1], run_72_40_inv_steady[1,-1], run_96_40_inv_steady [1,-1]]
OAA = [run_0_40[2,-1], run_24_40_inv_steady [2,-1], run_48_40_inv_steady [2,-1], run_72_40_inv_steady[2,-1], run_96_40_inv_steady [2,-1]]
AcCoA = [run_0_40[3,-1], run_24_40_inv_steady [3,-1], run_48_40_inv_steady [3,-1], run_72_40_inv_steady[3,-1], run_96_40_inv_steady [3,-1]]
Cit = [run_0_40[4,-1], run_24_40_inv_steady [4,-1], run_48_40_inv_steady [4,-1], run_72_40_inv_steady[4,-1], run_96_40_inv_steady [4,-1]]
KG = [run_0_40[5,-1], run_24_40_inv_steady [5,-1], run_48_40_inv_steady [5,-1], run_72_40_inv_steady [5,-1], run_96_40_inv_steady [5,-1]]
ATP = [run_0_40[6,-1], run_24_40_inv_steady [6,-1], run_48_40_inv_steady [6,-1], run_72_40_inv_steady [6,-1], run_96_40_inv_steady [6,-1]]
psi = [run_0_40[7,-1], run_24_40_inv_steady [7,-1], run_48_40_inv_steady [7,-1], run_72_40_inv_steady [7,-1], run_96_40_inv_steady [7,-1]]

hf.plot_hpi_concentration(Cit, title='Citrate')

# Plot the Concentration matrix (Note: It make take some time to plot, as there are a lot of points)
hf.plot_concentration_matrix(run_24_40_inv_steady, run_24_40_inv_steady_tau, title='Concentration Curves, 24 HPI')
