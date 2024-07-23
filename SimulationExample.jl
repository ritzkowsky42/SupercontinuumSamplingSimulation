#using Plots
using PyPlot
PyPlot.matplotlib[:use]("TkAgg")
const cm = 1/2.54


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 12

using Colors

using BenchmarkTools
include("TDSEfcn.jl")


# Simulation Parameters

N = 2000 # Number of points in spatial domain
num_steps = 10000 # Number of time steps


# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-50/t0,50/t0,num_steps)) # Time domain (not used

Δt = t[2]-t[1] # Time step size



# Define Step Potential 


E_fermi = 5 #0.05 # incident energy Fermi level
a = 20/x0 # barrier width in x0
cB = 1 # half width of barrier bottom smoothed Region
Wf = 5.10 # work function in eV
Bias_field = 0 # bias field in V/nm

# Define the width and strength of the absorbing layer
d_abl = 0.5/x0 # or some other appropriate width
V_max = -100 # adjust this value as needed

# Put into a struct

pot1 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max, Type="step", ROC = 15/x0, FE = 15)



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 15 / t0 # fwhm of pulse in fs
yc = 2060 / x0 # central wavelength in nm
phi_ce = 0#π # carrier envelope phase in radians
F = 15 # Field strength in V/nm

pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = 1)
pulse2 = pulseParams(fwhm = fwhm, yc = yc/2, phi_ce = phi_ce, F = 0)

field1 = pulsefromStruct(t, pulse1)
field2 = pulsefromStruct(t, pulse2)

w2w = addPulse(field1, field2)

plot(t*t0,w2w.E)
show()
#V = zeros(Complex{Float64},num_steps, N)

Vosc = genCoulombPotentialfromStruct(pot1,field1, x, num_steps, N)





u0 = statSolve(x, Vosc[1,:], E_fermi)



fig5,ax5 = subplots(1,1,figsize=(16*cm,16*cm))

ax5.plot(x*x0,Vosc[1,:],color="black")
ax5.set_xlabel("Position in (nm)")
ax5.set_ylabel("Potential in (eV)")
ax5.set_title("Potential Barrier")
ax5.set_xlim([-2,5])

ax6 = ax5.twinx()
ax6.semilogy(x*x0,abs2.(u0),color="red")
ax6.set_ylabel("Probability Density in (arb.u.)")
ax6.set_xlim([-2,5])
ax6.set_ylim([1e-40,5])
show()

# Solve the TDSE
ψ_stored = zeros(Complex{Float64},num_steps, N)

@time runSimulation2!(u0,ψ_stored,E_fermi, Vosc, Δx, Δt, N, num_steps)

# Calculate the current
j = current(x,ψ_stored)
jxr = currentSingle(x,ψ_stored,xr)
# Plot the results

# Here i add a two row Makie plot showing the probability density and the current j in a heatmap.





skipSteps = 10

fig1,(ax1,ax2) = subplots(1,2,figsize=(16*cm,8*cm),sharex=true)

pcm = ax1.pcolormesh(t[1:skipSteps:end]*t0,x*x0,abs2.(ψ_stored[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-1e-4, vmax=1e-4), shading="auto")

fig1.colorbar(pcm, ax=ax1, extend="max")

ax1.set_ylabel("Position in (nm)")
ax1.set_xlabel("Time in (fs)")
ax1.set_title("Probability Density")
ax1.set_ylim([0,5])
ax1.set_xlim([-20,20])
ax1.invert_xaxis()



pcm2 = ax2.pcolormesh(t[1:skipSteps:end]*t0,x*x0,real(j[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-1e-4, vmax=1e-4))
ax2.set_ylabel("Position in (nm)")
ax2.set_xlabel("Time in (fs)")
ax2.set_title("Probability Current")
ax2.set_ylim([0,5])
ax2.set_xlim([-20,20])
ax2.invert_xaxis()

fig1.colorbar(pcm2, ax=ax2, extend="max")


fig1.tight_layout(pad=0.2)
show()



fig2,ax3 = subplots(1,1,figsize=(16*cm,12*cm),sharex=true)

ax3.plot(t*t0,jxr )
ax3.set_xlabel("Time in (fs)")
ax3.set_ylabel("Current in (arb.u.)")
ax3.set_ylim([-1e-3,1e-3])
#ax2.set_xlim([4,15])}
ax3.set_title("Instantaneous Current at x = $(xr*x0) nm")

ax4 = ax3.twinx()

ax4.plot(t*t0,w2w.E,color="red")
ax4.set_ylim([-1,1])

fig2.savefig("InstantaneousCurrent.png",dpi=600)
show()