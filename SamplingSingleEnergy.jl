#using Plots
using PyPlot
PyPlot.matplotlib[:use]("TkAgg")
const cm = 1/2.54


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 12

using BenchmarkTools
using DSP
# dir = pwd()
# # cd("/Users/felix/Coding/TDSEjulia/")
# # cd("C:\\Users\\ritzk\\Coding\\TDSEjulia.jl")
# # pwd()
# include("C:\\Users\\ritzk\\Coding\\TDSEjulia.jl\\TDSEfcn.jl")
# # cd(dir)
using Keldysh

using CSV
using DataFrames
using FFTW
using DSP
using Interpolations

# loading pump pulse


# Simulation Parameters

N = 500 # Number of points in spatial domain
num_steps = 50000 # Number of time steps


# Define spatial and time domains

x = collect(LinRange(-2/x0, 10/x0, N)) # Spatial domain
xStart = x[1]
xEnd = x[end]

Δx = x[2] - x[1] # Spatial step size

t= collect(LinRange(-150/t0,150/t0,num_steps)) # Time domain (not used

Δt = t[2]-t[1] # Time step size



# Define Step Potential 


E_fermi = 0.05 # incident energy Fermi level
a = 15/x0 # barrier width in x0
cB = 1 # half width of barrier bottom smoothed Region
Wf = 13 # work function in eV
Bias_field = 0 # bias field in V/nm

# Define the width and strength of the absorbing layer
d_abl = 0.5/x0 # or some other appropriate width
V_max = -100 # adjust this value as needed

# Put into a struct

pot1 = PotentialParams(E_fermi = E_fermi, a = a, cB = cB, Wf = Wf, Bias_field = Bias_field, d_abl = d_abl, V_max = V_max , Type = "step", ROC = 15/x0, FE = 15)



xr = 5/x0 #  measurement reference plane in x0


# Define optical pulse

fwhm = 80 / t0 # fwhm of pulse in fs
yc = 1690 / x0 # central wavelength in nm
phi_ce = π # carrier envelope phase in radians
F = 15 # Field strength in V/nm

pulse1 = pulseParams(fwhm = fwhm, yc = yc, phi_ce = phi_ce, F = F)

field1 = pulsefromStruct(t, pulse1)

#V = zeros(Complex{Float64},num_steps, N)
 
Vosc = genPotential(pot1,field1,x, num_steps, N)


u0 = statSolve(x, Vosc[1,:], E_fermi)



# logspaced values of the field strength

delaySteps = 1600


delayVals = collect(LinRange(-400/t0 , 400/t0, delaySteps))


pump = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

interpPump = extrapolate(interpolate((real.(pump[:,1]),),real.(pump[:,2]), Gridded(Linear())),0 ) 



signal = CSV.read("FilteredWaveformSCG.csv", DataFrame,types=Complex{Float64})

interpSignal = extrapolate(interpolate((real.(signal[:,1]),),real.(signal[:,2]), Gridded(Linear())),0 ) 


# pulsePump = pulseParams(fwhm = 15/t0, t0 =0/t0, yc = yc, phi_ce = 0, F = 10)

# signalDuration = 8/t0

# pumpPulseField = pulsefromStruct(t,pulsePump)

# pulseList = []

# ysignal = 1690/x0

# for val in delayVals
#     #push!(pulseList,tukeyPulse(t,1,val,0.5))
#     push!(pulseList,pulsefromStruct(t,pulseParams(fwhm = signalDuration, t0=val, yc = ysignal, phi_ce = 0, F = 0.1)))
# end

fieldList = []

tukeyWindow = DSP.Windows.tukey(length(t), 0.2)

plot(t*t0,tukeyWindow.*interpPump(t*t0))
show()


for (i,val) in enumerate(delayVals)
    tempField = tukeyWindow.*(interpPump(t*t0) - 0.01.*interpSignal(t*t0 .-val*t0))
    push!(fieldList,opticalField(tempField,abs.(tempField),0,1))
end




@time jCEP = pulseSweep(u0,pot1,fieldList,E_fermi,5/x0,x,Δx,Δt,N,num_steps)


plot(t*t0,(jCEP[:,250]))
show()



# function FNemission(inField, criticalField, fieldEnhancement)
#     heaviside(x) = ifelse(real(x) >= 0, 1, 0)
#     out = heaviside.(real.(inField)) .* (real.(inField) .* fieldEnhancement).^2 .* exp.(-criticalField ./ abs.(fieldEnhancement .* real.(inField)))
#     # out = heaviside.(real.(inField)) .* (real.(inField) .* fieldEnhancement).^6
#     return out
#   end


# function FNemissionDerivative(inField, criticalField, fieldEnhancement)
#     derCurrent = -ifelse(real(inField) >= 0, 1, 0) * 2 * ((real(inField) * fieldEnhancement^2) * exp(-criticalField / abs(fieldEnhancement * real(inField))) - fieldEnhancement^2 * criticalField * exp(-criticalField / abs(fieldEnhancement * real(inField))))
#     return derCurrent
# end

# analyticalCurrent = zeros(length(delayVals))

# for (i,val) in enumerate(delayVals)
#     analyticalCurrent[i] = sum(FNemission.(fieldList[i].E.*15, 76, 1))
# end

# analyticalCurrent = analyticalCurrent .- mean(analyticalCurrent)

# analyticalCurrent = analyticalCurrent./maximum(analyticalCurrent)

samplingCurrent = zeros(length(delayVals))
samplingCurrent = (sum(jCEP,dims=1)' .- mean(sum(jCEP,dims=1)[1,1:100]'))./maximum((sum(jCEP,dims=1)' .- mean(sum(jCEP[:,:],dims=1)[1,1:100]')))

# #refpulse = tukeyPulse(delayVals,1,0,0.7) #
# refpulse= pulsefromStruct(delayVals,pulseParams(fwhm = signalDuration, t0=0, yc = ysignal, phi_ce = 0, F = 1))

fig2,ax2 = subplots(1,1,figsize=(16*cm,16*cm))

ax2.plot(delayVals*t0,reverse(samplingCurrent))
ax2.plot(delayVals*t0,interpSignal(delayVals*t0))

# ax2.loglog(fieldVals,0.7e-9 .* fieldVals .^ 10 )
# ax2.set_xlabel("Field strength in (V/nm)")
ax2.set_ylabel("Current in (arb.u.)")
#ax2.set_ylim([1e-10,1e4])
#ax2.set_xlim([-20,20])
ax2.legend(["sampled","Reference Pulse","Analytical"])
fig2.savefig("SamplingVsRef.png",dpi=600)
show()


# specSample = fft(fftshift(samplingCurrent))
# specRef= fft(fftshift((refpulse.E)))

# specAnalytical = fft(fftshift(analyticalCurrent))


# freqs = fftfreq(length(delayVals),1/((delayVals[2]-delayVals[1])*t0))
# freqs2= fftfreq(length(t),1/((t[2]-t[1])*t0))

# fig3,(ax3,ax4) = subplots(2,1,figsize = (16*cm,16*cm),sharex=true)

# ax3.semilogy(freqs[begin:Int(delaySteps/2)],abs.(specSample[begin:Int(delaySteps/2)])./maximum(abs.(specSample[begin:Int(delaySteps/2)])),label="Sampled Pulse")
# ax3.semilogy(freqs[begin:Int(delaySteps/2)],abs.(specRef[begin:Int(delaySteps/2)])./maximum(abs.(specRef[begin:Int(delaySteps/2)])),label="Reference Pulse")
# ax3.semilogy(freqs[begin:Int(delaySteps/2)],abs.(specAnalytical[begin:Int(delaySteps/2)])./maximum(abs.(specAnalytical[begin:Int(delaySteps/2)])),label="Analytical Pulse")
# ax3.set_xlim([0,2])
# ax3.set_ylim([1e-4,1])
# ax3.legend()

# #ax4 = ax3.twinx()
# ax4.plot(freqs[begin:Int(delaySteps/2)],unwrap(angle.(specSample[begin:Int(delaySteps/2)]))/π,label="Sampled Pulse")
# ax4.plot(freqs[begin:Int(delaySteps/2)],unwrap(angle.(specRef[begin:Int(delaySteps/2)]))/π,label="Reference Pulse")
# ax4.plot(freqs[begin:Int(delaySteps/2)],unwrap(angle.(specAnalytical[begin:Int(delaySteps/2)]))/π,label="Analytical Pulse")
# ax4.set_xlabel("Frequency in (PHz)")
# ax4.set_ylabel("Phase in π")
# ax4.set_ylim([-5,5])
# ax4.legend()
# fig3.savefig("SamplingVsRefSpec.png",dpi=600)
# show()



# TransferFCN = specSample./specRef
# TransferFCNanalytical = specAnalytical./specRef

# TransferFCNDerivative = FNemissionDerivative.(pumpPulseField.E, 76, 1)

# semilogy(freqs[begin:Int(delaySteps/2)],abs.(TransferFCN[begin:Int(delaySteps/2)]))
# semilogy(freqs[begin:Int(delaySteps/2)],abs.(TransferFCNanalytical[begin:Int(delaySteps/2)]))
# #semilogy(freqs2[begin:Int(round(num_steps/2))],abs.(TransferFCNDerivative[begin:Int(round(num_steps/2))])./maximum(abs.(TransferFCNDerivative[begin:Int(round(num_steps/2))])))
# # xlim([0,2])
# # ylim([1e-4,1e5])
# show()