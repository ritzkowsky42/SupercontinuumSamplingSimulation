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

N = 500 # Number of points in spatial domain
num_steps = 5000 # Number of time steps

k_B = 8.617e-5 # Boltzmann constant in eV/K

# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-50/t0,50/t0,num_steps)) # Time domain (not used

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

pot1 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max,Type="step", ROC = 10/x0, FE = 15)



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 15 / t0 # fwhm of pulse in fs
yc = 1000 / x0 # central wavelength in nm
phi_ce = 0#π # carrier envelope phase in radians
F = 1 # Field strength in V/nm

pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = F)

field1 = pulsefromStruct(t, pulse1)

#V = zeros(Complex{Float64},num_steps, N)

Vosc = genPotential(pot1,field1, x, num_steps, N)


# Define the ground state energies to be sampled
nEnergy= 500

E = collect(LinRange(0,8,nEnergy)) # Energy in eV
T = 293.15 #293.15 # Room temperature in K

function fermiDistribution(E::Float64, E_fermi::Float64, T::Float64)
    return 1/(1+exp((E-E_fermi)/(k_B*T)))
end

u0 = statSolve(x, Vosc[1,:], E_fermi)

fermi = fermiDistribution.(E,E_fermi,T)


fermi2 = (1/(2*π).^2) * log.(1 .+exp.((E_fermi.-E)./(k_B*T)))

plot(E,fermi2)
show()



# groundStates = zeros(Complex{Float64},nEnergy, N)

# for i in 1:nEnergy
#     groundStates[i,:] = statSolve(x, Vosc[1,:], E[i])
# end



# fig5,ax5 = subplots(1,1,figsize=(16*cm,16*cm))

# ax5.plot(x*x0,Vosc[1,:],color="black")
# ax5.plot(x*x0,Vosc[Int(num_steps/2),:],color="black",linestyle="--")

# ax5.set_xlabel("Position in (nm)")
# ax5.set_ylabel("Potential in (eV)")
# ax5.set_title("Potential Barrier")
# ax5.set_xlim([-2,5])

# ax6 = ax5.twinx()
# ax6.semilogy(x*x0,abs2.(u0),color="red")
# ax6.set_ylabel("Probability Density in (arb.u.)")
# ax6.set_xlim([-2,5])
# ax6.set_ylim([1e-40,5])
# show()

# Solve the TDSE
ψ_stored = zeros(Complex{Float64},num_steps, N)

# running the ground state Sweep

@time jSweep = groundStateSweep(Vosc,E,xr,Δx,Δt,N,num_steps)

totalCharge = sum(jSweep.*Δt,dims=1)'
plot(E,totalCharge)


densityofStates = (2*1.1*511e3/c^2).^(3/2) * sqrt.(E) ./ (2*π^2 *ħ^3)
densityofStates = densityofStates./maximum(densityofStates)

probabilityDensity = densityofStates .* fermi


emissionDistribution = probabilityDensity .* totalCharge./sum(probabilityDensity .* totalCharge)



#plot(E,totalCharge)
# plot(E,probabilityDensity .* totalCharge)


# plot(E,emissionDistribution)
plot(E,totalCharge.*fermi2)
#plot(E,totalCharge)
show()

# rescaled Instantaneous current

jTemp = jSweep .* (probabilityDensity)'


jrenorm = sum(jTemp.*(E[2]-E[1]),dims=2)

plot(t*t0,jrenorm./maximum(jrenorm))
plot(t*t0,jSweep[:,Int(nEnergy/2)]./maximum(jSweep[:,Int(nEnergy/2)]))
show()

jscaled = jSweep * probabilityDensity


plot(t*t0,jscaled./maximum(jscaled))
plot(t*t0,jSweep[:,Int(nEnergy/2)]./maximum(jSweep[:,Int(nEnergy/2)]))
show()




@time runSimulation2!(u0,ψ_stored,E_fermi, Vosc, Δx, Δt, N, num_steps)

# Calculate the current
j = current(x,ψ_stored)
jxr = currentSingle(x,ψ_stored,xr)
# Plot the results

# Here i add a two row Makie plot showing the probability density and the current j in a heatmap.





skipSteps = 1

fig1,(ax1,ax2) = subplots(1,2,figsize=(16*cm,8*cm))

pcm = ax1.pcolormesh(t[1:skipSteps:end]*t0,x*x0,abs2.(ψ_stored[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-1e-15, vmax=1e-15), shading="auto")

fig1.colorbar(pcm, ax=ax1, extend="min")

ax1.set_ylabel("Position in (nm)")
ax1.set_xlabel("Time in (fs)")
ax1.set_title("Probability Density")
# ax1.set_ylim([0,10])
# ax1.set_xlim([-20,20])
ax1.invert_xaxis()

pcm2 = ax2.pcolormesh(t[1:skipSteps:end]*t0,x*x0,real(j[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-3e-15, vmax=3e-15))
ax2.set_ylabel("Position in (nm)")
ax2.set_xlabel("Time in (fs)")
ax2.set_title("Probability Current")
# ax2.set_ylim([0,10])
# ax2.set_xlim([-20,20])
ax2.invert_xaxis()

fig1.colorbar(pcm2, ax=ax2, extend="max")


fig1.tight_layout(pad=0.2)



fig2,ax3 = subplots(1,1,figsize=(16*cm,12*cm))

ax3.plot(t*t0,jxr )
ax3.set_xlabel("Time in (fs)")
ax3.set_ylabel("Current in (arb.u.)")
#ax2.set_xlim([4,15])}
ax3.set_title("Instantaneous Current at x = $(xr*x0) nm")
fig2.savefig("InstantaneousCurrent.png",dpi=600)
show()
