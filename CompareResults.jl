using FFTW

using PyPlot
using CSV
using DataFrames
using DSP
using Tables
using Interpolations
using MAT
using CurveFit

cm = 1/2.54


vars = matread("Waveforms/nondegenerate.mat")

vars


# Laden der CSV-Datei
MeasurementResult = CSV.read("Waveforms/LongGateSCG.csv", DataFrame)
ReferenceResult = CSV.read("FilteredWaveformSCG.csv", DataFrame,types=Complex{Float64})
PumpPulse = CSV.read("Waveforms/Degenerate1690.csv", DataFrame,types=Complex{Float64})
PumpFiltered = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

TDSE3 = CSV.read("SampledWaveform10Vnm.csv", DataFrame)
TDSE = CSV.read("SampledWaveform15Vnm.csv", DataFrame)
TDSE2 = CSV.read("SampledWaveform30Vnm.csv", DataFrame)
TDSE4 = CSV.read("SampledWaveform40Vnm.csv", DataFrame)



tdseTime = TDSE[:,1]
tdseField = (TDSE[:,2])./maximum(abs.(TDSE[:,2]))

tdseField2 = (TDSE2[:,2])./maximum(abs.(TDSE2[:,2]))
tdseField3 = (TDSE3[:,2])./maximum(abs.(TDSE3[:,2]))
tdseField4 = (TDSE4[:,2])./maximum(abs.(TDSE4[:,2]))

refField = (ReferenceResult[:,2]./maximum(abs.(ReferenceResult[:,2])))

interpSignal = extrapolate(interpolate((real.(ReferenceResult[:,1]),),real.(refField), Gridded(Linear())),0 ) 



pumpField = (PumpFiltered[:,2]./maximum(abs.(PumpFiltered[:,2])))
interpPump = extrapolate(interpolate((real.(PumpFiltered[:,1]),),real.(pumpField), Gridded(Linear())),0 ) 


function FNemission(inField, criticalField, fieldEnhancement)
    heaviside(x) = ifelse(real(x) >= 0, 1, 0)
    out = heaviside.(real.(inField)) .* (real.(inField) .* fieldEnhancement).^2 .* exp.(-criticalField ./ abs.(fieldEnhancement .* real.(inField)))
    # out = heaviside.(real.(inField)) .* (real.(inField) .* fieldEnhancement).^6
    return out
end


function FNemissionDerivative(inField, criticalField, fieldEnhancement)
    derCurrent = -ifelse(real(inField) >= 0, 1, 0) * 2 * ((real(inField) * fieldEnhancement^2) * exp(-criticalField / abs(fieldEnhancement * real(inField))) - fieldEnhancement^2 * criticalField * exp(-criticalField / abs(fieldEnhancement * real(inField))))
    return derCurrent
end

delaySteps = 1600
timeSteps = 100000 
timeVec =collect(LinRange(-2000, 2000, timeSteps))

delayVals = collect(LinRange(-400, 400, delaySteps))

analyticalCurrent = zeros(length(delayVals))

for (i,val) in enumerate(delayVals)
    analyticalCurrent[i] = sum(FNemission.(interpPump(timeVec)+sqrt(3300).*interpSignal(timeVec.-val), 76, 15))
end

analyticalCurrent = analyticalCurrent .- sum(analyticalCurrent[1:20])./length(analyticalCurrent[1:20])

analyticalCurrent = analyticalCurrent./maximum(analyticalCurrent)



meas1Time =vars["nondegen_20230415_10ms_longgateshortsig_time"]
meas1Field = vars["nondegen_20230415_10ms_longgateshortsig_total"]
meas1Cmplx = DSP.Util.hilbert(meas1Field)
meas1Cmplx = meas1Cmplx./maximum(abs.(meas1Cmplx))

meas2Time =vars["nondegen_20230403_10ms_longgateshortsig_time"]
meas2Field = vars["nondegen_20230403_10ms_longgateshortsig_total"]
meas2Cmplx = DSP.Util.hilbert(meas2Field)
meas2Cmplx = meas2Cmplx./maximum(abs.(meas2Cmplx))

meas3Time =vars["nondegen_20230329_10ms_longgateshortsig_time"]
meas3Field = vars["nondegen_20230329_10ms_longgateshortsig_total"]
meas3Cmplx = DSP.Util.hilbert(meas3Field)
meas3Cmplx = meas3Cmplx./maximum(abs.(meas3Cmplx))

fig, ax = subplots(1, 1, figsize=(12, 8))  # Größe in Zoll
ax.plot(meas1Time,real.(meas1Cmplx)./maximum(abs.(meas1Cmplx)),label="20230415")
ax.plot(meas2Time.+304.5,real.(meas2Cmplx)./maximum(abs.(meas2Cmplx)).+2,label="20230403")
ax.plot(meas3Time.+339,real.(meas3Cmplx)./maximum(abs.(meas3Cmplx)).+4,label="20230329")

# ax.plot(vars["nondegen_20230403_10ms_longgateshortsig_time"],vars["nondegen_20230403_10ms_longgateshortsig_total"].+10,label="20230403")
# ax.plot(vars["nondegen_20230329_10ms_longgateshortsig_time"],vars["nondegen_20230329_10ms_longgateshortsig_total"].+20,label="20230329")
ax.legend()
ax.set_xlabel("Time (fs)")
ax.set_ylabel("Field (arb.u.)")
fig.savefig("DifferentMeasurements.png",dpi=600)
show()






fig, ax = subplots(1, 1, figsize=(12, 8))  # Größe in Zoll

ax.plot(MeasurementResult[:, 1], MeasurementResult[:, 2].+2,label="Measurement")
ax.plot(tdseTime, (tdseField),label="TDSE, 15 V/nm")
ax.plot(meas1Time,real.(meas1Cmplx)./maximum(abs.(meas1Cmplx)).+2,label="20230415")
ax.plot(meas2Time.+304.5,real.(meas2Cmplx)./maximum(abs.(meas2Cmplx)).+4,label="20230403")
ax.plot(meas3Time.+339,real.(meas3Cmplx)./maximum(abs.(meas3Cmplx)).+6,label="20230329")
# ax.plot(tdseTime, (tdseField2),label="30 V/nm")
# ax.plot(tdseTime, (tdseField3),label="10 V/nm")
# ax.plot(tdseTime, (tdseField4),label="40 V/nm")
ax.plot(delayVals,analyticalCurrent.+8,label="Analytical")
ax.plot(ReferenceResult[:, 1], refField.+10,label="Reference",color="black")
ax.set_xlabel("Time (fs)")
ax.set_ylabel("Field (arb.u.)")
ax.legend()
fig.savefig("TDSEvsMeasurement.png",dpi=600)
show()

# Converting time domain to Frequency domain
spectrumMeasurement = fft(MeasurementResult[:, 2]./maximum(abs.(MeasurementResult[:, 2])))
spectrumTDSE = fft(fftshift(tdseField))
spectrumTDSE2 = fft(fftshift(tdseField2))
spectrumTDSE3 = fft(fftshift(tdseField3))
spectrumTDSE4 = fft(fftshift(tdseField4))

spectrumMeas1 = vec(fft(fftshift(real.(meas1Cmplx))))
spectrumMeas2 = vec(fft(fftshift(real.(meas2Cmplx))))
spectrumMeas3 = vec(fft(fftshift(real.(meas3Cmplx))))


spectrumTDSE 
spectrumMeas1

freqsMeas1 = fftfreq(length(meas1Time), 1 / (meas1Time[2] - meas1Time[1]))
freqsMeas2 = fftfreq(length(meas2Time), 1 / (meas2Time[2] - meas2Time[1]))
freqsMeas3 = fftfreq(length(meas3Time), 1 / (meas3Time[2] - meas3Time[1]))


spectrumAnalytical = fft(fftshift(analyticalCurrent))
spectrumReference = fft(fftshift(refField))





# Creating the frequency axis
freqsMeasurement = fftfreq(length(MeasurementResult[:, 1]), 1 / (MeasurementResult[2, 1] - MeasurementResult[1, 1]))
freqsTDSE = fftfreq(length(tdseTime), 1 / (tdseTime[2] - tdseTime[1]))
freqRef = fftfreq(length(ReferenceResult[:, 1]), 1 / (ReferenceResult[2, 1] - ReferenceResult[1, 1]))
freqsAnalytical = fftfreq(length(delayVals), 1 / (delayVals[2] - delayVals[1]))


f0 = 1e-15*2.98e8/1690e-9
fStart = 0.225
fStop = 0.15

# Finding the index of reference wavelength
refIndexMeas1 = argmin(abs.(freqsMeas1 .- f0))
refIndexStart = argmin(abs.(freqsMeas1 .- fStart))
refIndexStop = argmin(abs.(freqsMeas1 .- fStop))

refIndexMeas2 = argmin(abs.(freqsMeas2 .- f0))
refIndexStart2 = argmin(abs.(freqsMeas2 .- fStart))
refIndexStop2 = argmin(abs.(freqsMeas2 .- fStop))

refIndexMeas3 = argmin(abs.(freqsMeas3 .- f0))
refIndexStart3 = argmin(abs.(freqsMeas3 .- fStart))
refIndexStop3 = argmin(abs.(freqsMeas3 .- fStop))



refIndexTDSE = argmin(abs.(freqsTDSE .- f0))
refIndexStopTDSE = argmin(abs.(freqsTDSE .- fStart))
refIndexStartTDSE = argmin(abs.(freqsTDSE .- fStop))


refIndexAnalytical = argmin(abs.(freqsAnalytical .- f0))
refIndexStopAnalytical = argmin(abs.(freqsAnalytical .- fStart))
refIndexStartAnalytical = argmin(abs.(freqsAnalytical .- fStop))


refIndexRef = argmin(abs.(freqRef .- f0))
refIndexStopRef = argmin(abs.(freqRef .- fStart))
refIndexStartRef = argmin(abs.(freqRef .- fStop))


# Calculating spectral phases

phaseTDSE1 = unwrap(angle.(spectrumTDSE))/π
phaseTDSE1 = phaseTDSE1 .- phaseTDSE1[refIndexTDSE]
linFitTDSE1 = curve_fit(LinearFit, freqsTDSE[refIndexStartTDSE:refIndexStopTDSE], phaseTDSE1[refIndexStartTDSE:refIndexStopTDSE])
phase1TDSE = phaseTDSE1 .- linFitTDSE1.(freqsTDSE)




phaseTDSE2 = unwrap(angle.(spectrumTDSE2))/π
phaseTDSE2 = phaseTDSE2 .- phaseTDSE2[refIndexTDSE]
linFitTDSE2 = curve_fit(LinearFit, freqsTDSE[refIndexStartTDSE:refIndexStopTDSE], phaseTDSE2[refIndexStartTDSE:refIndexStopTDSE])
phase2TDSE = phaseTDSE2 .- linFitTDSE2.(freqsTDSE)



phaseTDSE3 = unwrap(angle.(spectrumTDSE3))/π
phaseTDSE3 = phaseTDSE3 .- phaseTDSE3[refIndexTDSE]
linFitTDSE3 = curve_fit(LinearFit, freqsTDSE[refIndexStartTDSE:refIndexStopTDSE], phaseTDSE3[refIndexStartTDSE:refIndexStopTDSE])
phase3TDSE = phaseTDSE3 .- linFitTDSE3.(freqsTDSE)




phaseTDSE4 = unwrap(angle.(spectrumTDSE4))/π
phaseTDSE4 = phaseTDSE4 .- phaseTDSE4[refIndexTDSE]
linFit4 = curve_fit(LinearFit, freqsTDSE[refIndexStartTDSE:refIndexStopTDSE], phaseTDSE4[refIndexStartTDSE:refIndexStopTDSE])
phase4TDSE = phaseTDSE4 .- linFit4.(freqsTDSE)




phaseMeas1 = unwrap(angle.(spectrumMeas1))/π
phaseMeas1 = phaseMeas1 .- phaseMeas1[refIndexMeas1]
linFitMeas1 = curve_fit(LinearFit, freqsMeas1[refIndexStart:refIndexStop], phaseMeas1[refIndexStart:refIndexStop])
phaseMeas1 = phaseMeas1 .- linFitMeas1.(freqsMeas1)


phaseMeas2 = unwrap(angle.(spectrumMeas2))/π
phaseMeas2 = phaseMeas2 .- phaseMeas2[refIndexMeas2]
linFitMeas2 = curve_fit(LinearFit, freqsMeas2[refIndexStart2:refIndexStop2], phaseMeas2[refIndexStart2:refIndexStop2])
phaseMeas2 = phaseMeas2 .- linFitMeas2.(freqsMeas2)


phaseMeas3 = unwrap(angle.(spectrumMeas3))/π
phaseMeas3 = phaseMeas3 .- phaseMeas3[refIndexMeas3]
linFitMeas3 = curve_fit(LinearFit, freqsMeas3[refIndexStart3:refIndexStop3], phaseMeas3[refIndexStart3:refIndexStop3])
phaseMeas3 = phaseMeas3 .- linFitMeas3.(freqsMeas3)


phaseAnalytical = unwrap(angle.(spectrumAnalytical))/π
phaseAnalytical = phaseAnalytical .- phaseAnalytical[refIndexAnalytical]
linFitAnalytical = curve_fit(LinearFit, freqsAnalytical[refIndexStartAnalytical:refIndexStopAnalytical], phaseAnalytical[refIndexStartAnalytical:refIndexStopAnalytical])
phaseAnalytical = phaseAnalytical .- linFitAnalytical.(freqsAnalytical)



phaseReference = unwrap(angle.(spectrumReference))/π
phaseReference = phaseReference .- phaseReference[refIndexRef]
linFitRef = curve_fit(LinearFit, real.(freqRef[refIndexStartRef:refIndexStopRef]), phaseReference[refIndexStartRef:refIndexStopRef])
phaseReference = phaseReference .- linFitRef.(freqRef)



fig,(ax,ax2) = subplots(2,1,figsize = (20*cm,16*cm),sharex=true)
# ax.semilogy(freqsMeasurement, abs.(spectrumMeasurement).^2,label="Measurement")
ax.semilogy(freqsTDSE, 1e2.*abs.(spectrumTDSE).^2,linestyle="--",label="TDSE, 15 V/nm")
# ax.semilogy(freqsTDSE, 1e2.*abs.(spectrumTDSE2).^2,label="30 V/nm")
# ax.semilogy(freqsTDSE, 1e2.*abs.(spectrumTDSE3).^2,label="10 V/nm")
# ax.semilogy(freqsTDSE, 1e2.*abs.(spectrumTDSE4).^2,label="40 V/nm")
ax.semilogy(freqsMeas1, abs.(spectrumMeas1).^2,label="20230415")
ax.semilogy(freqsMeas2, abs.(spectrumMeas2).^2,label="20230403")
ax.semilogy(freqsMeas3, abs.(spectrumMeas3).^2,label="20230329")

ax.semilogy(freqsAnalytical, 1e2.*abs.(spectrumAnalytical).^2,linestyle="--",label="Analytical")
ax.semilogy(freqRef, abs.(spectrumReference).^2,label="Reference",color="black")
ax.set_xlim([0.13,0.275])
ax.set_ylim([1e1,1e8])
ax.set_xlabel("Frequency in (PHz)")
ax.set_ylabel("Spectral Intensity in (arb.u.)")
ax.legend()

ax2.plot(freqsTDSE, phase1TDSE,linestyle="--",label="TDSE, 15 V/nm")
ax2.plot(freqsMeas1, phaseMeas1,label="20230415")
ax2.plot(freqsMeas2, phaseMeas2,label="20230403")
ax2.plot(freqsMeas3, phaseMeas3,label="20230329")
ax2.plot(freqsAnalytical, phaseAnalytical,linestyle="--",label="Analytical")
ax2.plot(freqRef,phaseReference,label="Reference",color="black")


# ax2.plot(freqsTDSE[refIndexStartTDSE:refIndexStopTDSE], phase1TDSE[refIndexStartTDSE:refIndexStopTDSE],linestyle="--",label="TDSE, 15 V/nm")
# ax2.plot(freqsMeas1[refIndexStart:refIndexStop], phaseMeas1[refIndexStart:refIndexStop],label="20230415")
# ax2.plot(freqsMeas2[refIndexStart:refIndexStop], phaseMeas2[refIndexStart:refIndexStop],label="20230403")
# ax2.plot(freqsMeas3[refIndexStart:refIndexStop], phaseMeas3[refIndexStart:refIndexStop],label="20230329")
# ax2.plot(freqsAnalytical[refIndexStartAnalytical:refIndexStopAnalytical], phaseAnalytical[refIndexStartAnalytical:refIndexStopAnalytical],linestyle="--",label="Analytical")
# ax2.plot(freqRef[refIndexStartRef:refIndexStopRef],phaseReference[refIndexStartRef:refIndexStopRef],label="Reference",color="black")

ax2.set_xlabel("Frequency in (PHz)")
ax2.set_ylabel("Phase in π")
ax2.set_ylim([-5,5])

fig.savefig("TDSEvsMeasurementSpec.png",dpi=600)
show()
