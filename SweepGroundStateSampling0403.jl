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

u0 = statSolve(x, Vosc[1,:], E_fermi)


delaySteps = 1400


delayVals = collect(LinRange(-450/t0 , -300/t0, delaySteps))


pump = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

interpPump = extrapolate(interpolate((real.(pump[:,1]),),real.(pump[:,2]), Gridded(Linear())),0 ) 



signal = CSV.read("FilteredWaveformSCG0403.csv", DataFrame)

plot(real.(signal[:,1]),signal[:,2])
plot(reverse(real.(signal[:,1])),reverse(signal[:,2]).+2)
show()


interpSignal = extrapolate(interpolate((reverse(real.(signal[:,1])),),reverse(real.(signal[:,2])), Gridded(Linear())),0 ) 


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

tukeyWindow = DSP.Windows.tukey(length(t), 0.3)

plot(t*t0,tukeyWindow.*interpPump(t*t0))
show()


for (i,val) in enumerate(delayVals)
    tempField = tukeyWindow.*(interpPump(t*t0) + (1 ./sqrt(3300)).*interpSignal(t*t0 .-val*t0))
    push!(fieldList,opticalField(tempField,abs.(tempField),0,1))
end

# print("Single Energy Calculation")

# @time jCEP = pulseSweep(u0,pot1,fieldList,E_fermi,5/x0,x,Δx,Δt,N,num_steps)


# @time jSweep = pulseSweepSweep(pot1,fieldList,E,T,xr,Δx,Δt,N,num_steps)

# Define the ground state energies to be sampled

nEnergy= 20

E = collect(LinRange(2,6,nEnergy)) # Energy in eV
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


# Plotting the results



# samplingCurrent = zeros(length(delayVals))
# samplingCurrent = (sum(jCEP,dims=1)' .- mean(sum(jCEP,dims=1)[1,1:100]'))./maximum((sum(jCEP,dims=1)' .- mean(sum(jCEP[:,:],dims=1)[1,1:100]')))

# # convert array to vector
# samplingCurrent = vec(samplingCurrent)

# samplingCurrent = samplingCurrent .- mean(samplingCurrent[1:20])



# refpulse= pulsefromStruct(delayVals,pulseParams(fwhm = signalDuration, t0=0, yc = ysignal, phi_ce = 0, F = 1))


qAvg = -(q[:,1] .- mean(q[1:20,1]))

qAll = -(q[:,2:end] .- mean(q[1:20,2:end],dims=1))

CSV.write("SampledWaveform15Vnm0403Extension.csv", DataFrame(time = delayVals*t0, field = qAvg./maximum(qAvg)), writeheader = true)




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
#     analyticalCurrent[i] = sum(FNemission.(fieldList[i].E, 76, 15))
# end

# analyticalCurrent = analyticalCurrent .- mean(analyticalCurrent[1:20])

# analyticalCurrent = analyticalCurrent./maximum(analyticalCurrent)










# fig2,ax2 = subplots(1,1,figsize=(16*cm,16*cm))
# ax2.plot(delayVals*t0,qAvg./maximum(qAvg),label="Continuum State")
# ax2.plot(delayVals*t0,samplingCurrent,label="Single Energy E_F")
# ax2.plot(delayVals*t0,refpulse.E,label="Reference Pulse")
# ax2.plot(delayVals*t0,analyticalCurrent,label="Analytical FN")

# # ax2.loglog(fieldVals,0.7e-9 .* fieldVals .^ 10 )
# ax2.set_xlabel("Time (fs)")
# ax2.set_ylabel("Current in (arb.u.)")
# #ax2.set_ylim([1e-10,1e4])
# #ax2.set_xlim([-20,20])
# ax2.legend()
# fig2.savefig("SamplingVsRef.png",dpi=600)
# fig2.savefig("SamplingVsRef.pdf",dpi=600)
# show()



# pad = 1000
# zerosPad = zeros(pad)

# tukeyWin = DSP.Windows.tukey(delaySteps,0.002)

# plot(tukeyWin)
# show()
# # using hamming window to taper the edges of the pulse


# prepend!(qAvg,zerosPad)
# append!(qAvg,zerosPad)



# samplingCurrent = tukeyWin .* samplingCurrent

# prepend!(samplingCurrent,zerosPad)
# append!(samplingCurrent,zerosPad)

# refField = refpulse.E

# prepend!(refField,zerosPad)
# append!(refField,zerosPad)

# prepend!(analyticalCurrent,zerosPad)
# append!(analyticalCurrent,zerosPad)

# specSample = fft(fftshift(qAvg))
# specSample = specSample ./maximum(abs.(specSample))

# specSampleSingle = fft(fftshift(samplingCurrent))
# specSampleSingle = specSampleSingle ./maximum(abs.(specSampleSingle))

# specRef= fft(fftshift((refpulse.E)))
# specRef = specRef./maximum(abs.(specRef))

# specAnalytical = fft(fftshift(analyticalCurrent))
# specAnalytical = specAnalytical./maximum(abs.(specAnalytical))



# freqs = fftfreq(length(delayVals)+2*pad,1/((delayVals[2]-delayVals[1])*t0))
# freqs2= fftfreq(length(t)+2*pad,1/((t[2]-t[1])*t0))




# fig3,(ax3,ax4) = subplots(2,1,figsize = (20*cm,16*cm),sharex=true)

# halfSize = Int((delaySteps+2*pad)/2)

# ax3.semilogy(freqs[begin:halfSize],abs.(specSample[begin:halfSize])./maximum(abs.(specSample[begin:halfSize])),label="Continuum State")
# ax3.semilogy(freqs[begin:halfSize],abs.(specSampleSingle[begin:halfSize])./maximum(abs.(specSampleSingle[begin:halfSize])),label="Single Energy")
# ax3.semilogy(freqs[begin:halfSize],abs.(specRef[begin:halfSize])./maximum(abs.(specRef[begin:halfSize])),label="Reference Pulse")
# ax3.semilogy(freqs[begin:halfSize],abs.(specAnalytical[begin:halfSize])./maximum(abs.(specAnalytical[begin:halfSize])),label="Analytical FN")

# # ax3.semilogy(freqs[begin:Int(delaySteps/2)],abs.(specAnalytical[begin:Int(delaySteps/2)])./maximum(abs.(specAnalytical[begin:Int(delaySteps/2)])),label="Analytical Pulse")
# ax3.set_xlim([0,0.5])
# ax3.set_ylim([1e-8,1])
# ax3.set_ylabel("Amplitude in (arb.u.)")
# ax3.legend()
# ax3.legend(loc="lower right")

# #ax4 = ax3.twinx()
# ax4.plot(freqs[begin:halfSize],unwrap(angle.(specSample[begin:halfSize]))/π,label="Continuum State")
# ax4.plot(freqs[begin:halfSize],unwrap(angle.(specSampleSingle[begin:halfSize]))/π,label="Single Energy")
# ax4.plot(freqs[begin:halfSize],unwrap(angle.(specRef[begin:halfSize]))/π,label="Reference Pulse")
# ax4.plot(freqs[begin:halfSize],unwrap(angle.(specAnalytical[begin:halfSize]))/π,label="Analytical Pulse")
# # ax4.plot(freqs[begin:Int(delaySteps/2)],unwrap(angle.(specAnalytical[begin:Int(delaySteps/2)]))/π,label="Analytical Pulse")
# ax4.set_xlabel("Frequency in (PHz)")
# ax4.set_ylabel("Phase in (π)")
# ax4.set_ylim([-5,10])
# ax4.legend(loc="upper right")
# fig3.savefig("SamplingVsRefSpec.png",dpi=600)
# fig3.savefig("SamplingVsRefSpec.pdf",dpi=600)
# show()


# CSV.write("SampledWaveform.csv", DataFrame(time = delayVals*t0, field = qAvg./maximum(qAvg),ref=refpulse.E, writeheader = false))

# plot(delayVals*t0,-qAll./maximum(-qAll,dims=1))
# show()

# # transferFCN = specSample./specRef
# # trasnferFNCSingle = specSampleSingle./specRef
# # TransferFCNAnalytical = specAnalytical./specRef

# # fig4,ax4 = subplots(1,1,figsize=(16*cm,16*cm))
# # ax4.semilogy(freqs,abs.(transferFCN),label="Continuum State")
# # ax4.semilogy(freqs,abs.(trasnferFNCSingle),label="Single Energy")
# # ax4.semilogy(freqs,abs.(TransferFCNAnalytical),label="Analytical FN")

 
# # show()