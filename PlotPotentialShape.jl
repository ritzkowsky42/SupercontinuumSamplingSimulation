using PyPlot
PyPlot.matplotlib[:use]("TkAgg")
const cm = 1/2.54


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 8

using Colors
using ProgressBars
using BenchmarkTools

# dir = pwd()
# cd("/Users/felix/Coding/TDSEjulia/")
# include("../TDSEjulia/TDSEfcn.jl")
# cd(dir)
using Keldysh

using CSV
using DataFrames
using FFTW
using Interpolations
using DSP
using JLD2

# Simulation Parameters

N = 500 # Number of points in spatial domain
num_steps = 50000 # Number of time steps

k_B = 8.617e-5 # Boltzmann constant in eV/K

# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-150/t0,150/t0,num_steps)) # Time domain (not used

Δt = t[2]-t[1] # Time step size



# Define Step Potential 


E_fermi = 5.5 #0.05 # incident energy Fermi level
a = 20/x0 # barrier width in x0
cB = 1 # half width of barrier bottom smoothed Region
Wf = 5.10 # work function in eV
Bias_field = 0 # bias field in V/nm

stepHeight = E_fermi + Wf

# Define the width and strength of the absorbing layer
d_abl = 0.5/x0 # or some other appropriate width
V_max = -100 # adjust this value as needed

# Put into a struct

pot1 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max,Type="step", ROC = 15/x0, FE = 15)
pot2 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max,Type="step", ROC = 0, FE = 15)



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 80 / t0 # fwhm of pulse in fs
yc = 1690 / x0 # central wavelength in nm
phi_ce = 0#π # carrier envelope phase in radians
F = 1 # Field strength in V/nm

pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = F)

field1 = pulsefromStruct(t, pulse1)

#V = zeros(Complex{Float64},num_steps, N)

Vosc = genPotential(pot1,field1,x, num_steps, N)
Vosc2 = genPotential(pot2,field1,x, num_steps, N)

# Plotting the Potential

fig, ax = subplots(1,1, figsize = (8.5cm, 6cm))
ax.plot(x*x0, real(Vosc[25000,:]), label = "Real")
ax.plot(x*x0,real(Vosc2[25000,:]), label = "Real2")
ax.set_xlim([-2,5])
ax.set_ylim([-20,10])
ax.set_xlabel("x (nm)")
ax.set_ylabel("Potential Energy (eV)")