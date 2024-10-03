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
using MAT

#create time vector

t = collect(LinRange(-150,150,50000))

vars = matread("Waveforms/nondegenerate.mat")
vars


meas1Time = reverse(vec(vars["nondegen_20230415_10ms_longgateshortsig_time"]))
meas1Field = reverse(vec( vars["nondegen_20230415_10ms_longgateshortsig_total"]))

meas1Cmplx = DSP.Util.hilbert(meas1Field)
meas1Cmplx = meas1Cmplx./maximum(abs.(meas1Cmplx))


# plot(meas1Time,real.(meas1Cmplx))
# show()

interpMeas = extrapolate(interpolate((meas1Time,),real.(meas1Cmplx), Gridded(Linear())),0 )


pump = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

interpPump = extrapolate(interpolate((real.(pump[:,1]),),real.(pump[:,2]), Gridded(Linear())),0 ) 



signal = CSV.read("FilteredWaveformSCG0415.csv", DataFrame)

# plot(real.(signal[:,1]),signal[:,2])
# plot(reverse(real.(signal[:,1])),reverse(signal[:,2]).+2)
# show()


interpSignal = extrapolate(interpolate((reverse(real.(signal[:,1])),),reverse(real.(signal[:,2])), Gridded(Linear())),0 ) 



figure()
plot(t,real.(interpPump.(t)))
plot(t,real.(interpSignal.(t)).+2)
plot(t,real.(interpMeas.(t)).+4)
show()

# Convert all signals to complex signals


measC = DSP.Util.hilbert(real.(interpMeas.(t)))
measC = measC./maximum(abs.(measC))
measC .-= mean(measC)

pumpC = DSP.Util.hilbert(real.(interpPump.(t)))
pumpC = pumpC./maximum(abs.(pumpC))

signalC = DSP.Util.hilbert(real.(interpSignal.(t)))
signalC = signalC./maximum(abs.(signalC))



# create spectra 

pumpSpectrum = fftshift(fft(real.(pumpC)))
pumpSpectrum = pumpSpectrum./maximum(abs.(pumpSpectrum))
signalSpectrum = fftshift(fft(real.(signalC)))
signalSpectrum = signalSpectrum./maximum(abs.(signalSpectrum))
measSpectrum = fftshift(fft(real.(measC)))
measSpectrum = measSpectrum./maximum(abs.(measSpectrum))

# create frequency vector

dt = abs(t[2]-t[1])
f = fftshift(fftfreq(length(t),1/dt))


# Attempt deconvolution

conv = signalSpectrum


figure()
# plot(f,abs.(signalSpectrum),label="Signal")
# plot(f,abs.(measSpectrum), label="Measurement")
semilogy(f,abs.(deconv), label="Decon")
xlim(0,0.3)
show()

# # Laden der CSV-Datei
# MeasurementResult = CSV.read("Waveforms/LongGateSCG.csv", DataFrame)
# ReferenceResult = CSV.read("FilteredWaveformSCG.csv", DataFrame,types=Complex{Float64})
# PumpPulse = CSV.read("Waveforms/Degenerate1690.csv", DataFrame,types=Complex{Float64})
# PumpFiltered = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

# TDSE0329 = CSV.read("SampledWaveform15Vnm0329.csv", DataFrame)
# TDSE0403 = CSV.read("SampledWaveform15Vnm0403Extension.csv", DataFrame)
# TDSE0415 = CSV.read("SampledWaveform15Vnm0415Extension.csv", DataFrame)


# # loading reference Waveforms
# SCG0329 = CSV.read("FilteredWaveformSCG0329.csv",DataFrame)
# SCG0403 = CSV.read("FilteredWaveformSCG0403.csv",DataFrame)
# SCG0415 = CSV.read("FilteredWaveformSCG0415.csv",DataFrame)


