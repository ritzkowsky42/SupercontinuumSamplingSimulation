using FFTW

using PyPlot
using CSV
using DataFrames
using DSP
using Tables
cm = 1/2.54

# Laden der CSV-Datei
MeasurementResult = CSV.read("Waveforms/LongGateSCG.csv", DataFrame)
TDSE3 = CSV.read("SampledWaveform10Vnm.csv", DataFrame)
TDSE = CSV.read("SampledWaveform15Vnm.csv", DataFrame)
TDSE2 = CSV.read("SampledWaveform30Vnm.csv", DataFrame)


tdseTime = TDSE[:,1]
tdseField = (TDSE[:,2])./maximum(abs.(TDSE[:,2]))

tdseField2 = (TDSE2[:,2])./maximum(abs.(TDSE2[:,2]))
tdseField3 = (TDSE3[:,2])./maximum(abs.(TDSE3[:,2]))

fig, ax = subplots(1, 1, figsize=(12, 8))  # Größe in Zoll

ax.plot(MeasurementResult[:, 1], MeasurementResult[:, 2].+2,label="Measurement")
ax.plot(tdseTime, (tdseField),label="15 V/nm")
ax.plot(tdseTime, (tdseField2),label="30 V/nm")
ax.plot(tdseTime, (tdseField3),label="10 V/nm")
ax.set_xlabel("Time (fs)")
ax.set_ylabel("Field (arb.u.)")
ax.legend()
fig.savefig("TDSEvsMeasurement.png",dpi=600)
show()

# Converting time domain to Frequency domain
spectrumMeasurement = fft(MeasurementResult[:, 2]./maximum(abs.(MeasurementResult[:, 2])))
spectrumTDSE = fft(tdseField)
spectrumTDSE2 = fft(tdseField2)
spectrumTDSE3 = fft(tdseField3)

# Creating the frequency axis
freqsMeasurement = fftfreq(length(MeasurementResult[:, 1]), 1 / (MeasurementResult[2, 1] - MeasurementResult[1, 1]))
freqsTDSE = fftfreq(length(tdseTime), 1 / (tdseTime[2] - tdseTime[1]))

fig,ax = subplots(1,1,figsize = (20*cm,16*cm))
ax.semilogy(freqsMeasurement, abs.(spectrumMeasurement).^2,label="Measurement")
ax.semilogy(freqsTDSE, abs.(spectrumTDSE).^2,label="15 V/nm")
ax.semilogy(freqsTDSE, abs.(spectrumTDSE2).^2,label="30 V/nm")
ax.semilogy(freqsTDSE, abs.(spectrumTDSE3).^2,label="10 V/nm")
ax.set_xlim([0.15,0.3])
ax.set_xlabel("Frequency in (PHz)")
ax.set_ylabel("Spectral Intensity in (arb.u.)")
ax.legend()
fig.savefig("TDSEvsMeasurementSpec.png",dpi=600)
show()