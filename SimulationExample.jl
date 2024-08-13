#using Plots
using PyPlot
PyPlot.matplotlib[:use]("TkAgg")
const cm = 1/2.54


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 8
rcParams["pdf.fonttype"] = 42


using Colors

using BenchmarkTools
using Keldysh
using CSV
using DataFrames
using FFTW
using DSP
using Interpolations
# Simulation Parameters

N = 2000 # Number of points in spatial domain
num_steps = 200000 # Number of time steps


# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-300/t0,300/t0,num_steps)) # Time domain (not used

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

pot1 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max , Type = "step", ROC = 15/x0, FE = 15)



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 15 / t0 # fwhm of pulse in fs
yc = 2060 / x0 # central wavelength in nm
phi_ce = 0#π # carrier envelope phase in radians
F = 15 # Field strength in V/nm

# pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = 1)
# pulse2 = pulseParams(fwhm = fwhm, yc = yc/2, phi_ce = phi_ce, F = 0)

# #field1 = pulsefromStruct(t, pulse1)
# field2 = pulsefromStruct(t, pulse2)

# w2w = addPulse(field1, field2)

# plot(t*t0,w2w.E)
# show()
#V = zeros(Complex{Float64},num_steps, N)


pump = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

interpPump = extrapolate(interpolate((real.(pump[:,1]),),real.(pump[:,2]), Gridded(Linear())),0 )

signal = CSV.read("FilteredWaveformSCG0403.csv", DataFrame)

interpSignal = extrapolate(interpolate((reverse(real.(signal[:,1])),),reverse(real.(signal[:,2])), Gridded(Linear())),0 ) 





tukeyWindow = DSP.Windows.tukey(length(t), 0.3)

tempField = tukeyWindow.*(interpPump(t*t0))
tempField = DSP.Util.hilbert(tempField)
tempField .*= exp(-1im*0*pi)

tempField2 = tukeyWindow.*(interpSignal(t*t0))
tempField2 = DSP.Util.hilbert(tempField2)
tempField2 .*= 1 ./maximum(abs.(tempField2))
tempField2 .*= exp(-1im*0.5*pi)


plot(t*t0,tempField)
plot(t*t0,-tempField2.+2)
plot(t*t0,tukeyWindow)
show()

field1 = opticalField(-real.(tempField),abs.(tempField),0,1)

Vosc = genPotential(pot1,field1, x, num_steps, N)



field2 = opticalField(-real.(tempField2),abs.(tempField2),0,1)

Vosc2 = genPotential(pot1,field2, x, num_steps, N)





u0 = statSolve(x, Vosc[1,:], E_fermi)
u0 = statSolve(x, Vosc2[1,:], E_fermi)



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
ψ_stored2 = zeros(Complex{Float64},num_steps, N)

@time runSimulation2!(u0,ψ_stored,E_fermi, Vosc, Δx, Δt, N, num_steps)
@time runSimulation2!(u0,ψ_stored2,E_fermi, Vosc2, Δx, Δt, N, num_steps)

# Calculate the current
j = current(x,ψ_stored)
jxr = currentSingle(x,ψ_stored,2/x0)

j2 = current(x,ψ_stored2)
jxr2 = currentSingle(x,ψ_stored2,2/x0)




# Plot the results

skipSteps = 50

fig1,(ax1,ax3,ax2,ax4) = subplots(2,2,figsize=(16*cm,8*cm),sharex=true)

pcm = ax1.pcolormesh(t[1:skipSteps:end]*t0,x*x0,abs2.(ψ_stored[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-1e-2, vmax=1e-2), shading="auto",rasterized=true)

fig1.colorbar(pcm, ax=ax1, extend="max", label= "\$|\\psi(x,t)|^2\$")

ax1.set_ylabel("Position in (nm)")
# ax1.set_xlabel("Time in (fs)")
ax1.set_title("Probability Density")
ax1.set_ylim([0,5])
ax1.set_xlim([-50,50])
ax1.invert_xaxis()



pcm2 = ax2.pcolormesh(t[1:skipSteps:end]*t0,x*x0,real(j[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-5e-2, vmax=5e-2),rasterized=true)
# ax2.set_ylabel("Position in (nm)")
# ax2.set_xlabel("Time in (fs)")
ax2.set_title("Probability Current")
ax2.set_ylim([0,5])
ax2.set_xlim([-50,50])
ax2.invert_xaxis()
fig1.colorbar(pcm2, ax=ax2, extend="max",label="Current density")


pcm3 = ax3.pcolormesh(t[1:skipSteps:end]*t0,x*x0,abs2.(ψ_stored2[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-1e-2, vmax=1e-2), shading="auto",rasterized=true)

ax3.set_ylabel("Position in (nm)")
ax3.set_xlabel("Time in (fs)")
# ax3.set_title("Probability Density")
ax3.set_ylim([0,5])
ax3.set_xlim([-50,50])
ax3.invert_xaxis()

fig1.colorbar(pcm3, ax=ax3, extend="max", label= "\$|\\psi(x,t)|^2\$")



pcm4 = ax4.pcolormesh(t[1:skipSteps:end]*t0,x*x0,real(j2[end:-skipSteps:1,:])',cmap="RdBu",norm=matplotlib[:colors][:Normalize](vmin=-5e-2, vmax=5e-2),rasterized=true)
# ax4.set_ylabel("Position in (nm)")
ax4.set_xlabel("Time in (fs)")
# ax4.set_title("Probability Current")
ax4.set_ylim([0,5])
ax4.set_xlim([-50,50])
ax4.invert_xaxis()

fig1.colorbar(pcm4, ax=ax4, extend="max",label="Current density")


fig1.tight_layout(pad=0.2)
fig1.savefig("ProbabilityDensityCurrent.png",dpi=600)
fig1.savefig("ProbabilityDensityCurrent.pdf",dpi=600)
show()



fig2,ax3 = subplots(1,1,figsize=(16*cm,6*cm),sharex=true)

ax3.plot(t*t0,jxr./maximum(abs.(jxr)) ,label = "Long Pulse")
ax3.plot(t*t0,jxr2 ./maximum(abs.(jxr2)),label = "Short Pulse")
ax3.set_xlabel("Time in (fs)")
ax3.set_ylabel("Current in (arb.u.)")
ax3.set_ylim([-1.1,1.1])
ax3.set_xlim([-70,70])
ax3.set_title("Instantaneous Current at x = $(2) nm")
ax4 = ax3.twinx()

ax4.plot(t*t0,field1.E,color="k",alpha=0.5,label="Long Pulse")
ax4.plot(t*t0,field2.E,color="k",alpha=0.5,linestyle="--",label="Short Pulse")
ax4.set_ylim([-1,1])
ax4.set_ylabel("Field in (arb.u.)")
ax3.legend()
ax4.legend(loc="lower right")

fig2.tight_layout(pad=0.2)
fig2.savefig("InstantaneousCurrent.png",dpi=600)
fig2.savefig("InstantaneousCurrent.pdf",dpi=600)
show()