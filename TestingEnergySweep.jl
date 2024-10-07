using PyPlot
PyPlot.matplotlib[:use]("TkAgg")
const cm = 1/2.54


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 12

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
num_steps = 5000 # Number of time steps

k_B = 8.617e-5 # Boltzmann constant in eV/K

# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-30/t0,30/t0,num_steps)) # Time domain (not used

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



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 10 / t0 # fwhm of pulse in fs
yc = 2700 / x0 # central wavelength in nm
phi_ce = 0#π # carrier envelope phase in radians
F = 1 # Field strength in V/nm

pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = F)

field1 = pulsefromStruct(t, pulse1)

#V = zeros(Complex{Float64},num_steps, N)

Vosc = genPotential(pot1,field1,x, num_steps, N)

u0 = statSolve(x, Vosc[1,:], E_fermi)



# pulsePump = pulseParams(fwhm = 15/t0, t0 =0/t0, yc = yc, phi_ce = 0, F = 10)

# signalDuration = 8/t0

# pumpPulseField = pulsefromStruct(t,pulsePump)

# pulseList = []

# ysignal = 1690/x0

# for val in delayVals
#     #push!(pulseList,tukeyPulse(t,1,val,0.5))
#     push!(pulseList,pulsefromStruct(t,pulseParams(fwhm = signalDuration, t0=val, yc = ysignal, phi_ce = 0, F = 0.1)))
# end

cepVals = collect(LinRange(0,4*π,20))

fieldList = []


for (i,val) in enumerate(cepVals)
    push!(fieldList,opticalField(pulsefromStruct(t,pulseParams(fwhm = fwhm, yc = yc, phi_ce = val, F = F))))
end

# print("Single Energy Calculation")

# @time jCEP = pulseSweep(u0,pot1,fieldList,E_fermi,5/x0,x,Δx,Δt,N,num_steps)


# @time jSweep = pulseSweepSweep(pot1,fieldList,E,T,xr,Δx,Δt,N,num_steps)

# Define the ground state energies to be sampled

nEnergy= 20

E = collect(LinRange(0,6,nEnergy)) # Energy in eV
T = 293.15 #293.15 # Room temperature in K



n = length(fieldList)
q = zeros(Float64,n,nEnergy+1)
stepHeight = Wf + E_fermi

a =[]

for i in ProgressBar(1:n)
    VoscTemp = genPotential(pot1,fieldList[i] ,x, num_steps, N)
    jTemp = groundStateSweep(VoscTemp,E,xr,x,Δx,Δt,N,num_steps)
    qOut = sum(jTemp.*Δt,dims=1)'
    q[i,2:end] = qOut
    q[i,1],Fout,DistOut = gs_integrated_charge(qOut,stepHeight.-E , E_fermi ,stepHeight,T)
    VoscTemp = []
    jTemp = []
    push!(a,i*2)
end


plot(cepVals,q[:,1])
show()

q[:,1]