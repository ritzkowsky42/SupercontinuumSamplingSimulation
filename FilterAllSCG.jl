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

# Laden der CSV-Datei

vars = matread("Waveforms/scgmeasurement.mat")

vars



scg20230403 =vars["degenerateSCG_20230403_10ms_total"][:]
scg20230403Cmplx = DSP.Util.hilbert(scg20230403)
scg20230403Cmplx = scg20230403Cmplx./maximum(abs.(scg20230403Cmplx))

scg20230415 = vars["degenerateSCG_20230415_10ms_total"][:]
scg20230415Cmplx = DSP.Util.hilbert(scg20230415)
scg20230415Cmplx = scg20230415Cmplx./maximum(abs.(scg20230415Cmplx))

scg20230329 = vars["degenerateSCG_20230329_30ms_total"][:]
scg20230329Cmplx = DSP.Util.hilbert(scg20230329)
scg20230329Cmplx = scg20230329Cmplx./maximum(abs.(scg20230329Cmplx))


scg20230403Time = vars["degenerateSCG_20230403_10ms_time"][:]
scg20230415Time = vars["degenerateSCG_20230415_10ms_time"][:]
scg20230329Time = vars["degenerateSCG_20230329_30ms_time"][:]



# Creating windows for the time domain using Tukey Window

tukeyWindow0403 = DSP.Windows.tukey(length(scg20230403), 0.3)
tukeyWindow0415 = DSP.Windows.tukey(length(scg20230415), 0.3)
tukeyWindow0329 = DSP.Windows.tukey(length(scg20230329), 0.3)


# Erstellen der Figur und der Achsen
fig, ax = subplots(figsize=(16, 12))  # Größe in Zoll
# Hinzufügen der Datenlinien zum Plot
ax.plot(scg20230403Time, real.(scg20230403Cmplx), label="20230403")
ax.plot(scg20230403Time, real.(scg20230403Cmplx).*tukeyWindow0403, label="20230403 - Windowed")


ax.plot(scg20230415Time, real.(scg20230415Cmplx).+2, label="20230415")
ax.plot(scg20230415Time, real.(scg20230415Cmplx).*tukeyWindow0415.+2, label="20230415 - Windowed")

ax.plot(scg20230329Time, real.(scg20230329Cmplx).+4, label="20230329")
ax.plot(scg20230329Time, real.(scg20230329Cmplx).*tukeyWindow0329.+4, label="20230329 - Windowed")


# Beschriftungen hinzufügen
ax.set_xlabel("Time (fs)")
ax.set_ylabel("Field (arb.u.)")
ax.set_title("Waveform Plot")
ax.legend()
show()


# Anzeigen des Plots


# Converting time domain to Frequency domain

spectrumSCG0403 = fft(scg20230403Cmplx.*tukeyWindow0403)
spectrumSCG0415 = fft(scg20230415Cmplx.*tukeyWindow0415)
spectrumSCG0329 = fft(scg20230329Cmplx.*tukeyWindow0329)

# Creating the frequency axis
freqsSCG0403 = fftfreq(length(spectrumSCG0403), 1 / (scg20230403Time[2] - scg20230403Time[1]))
freqsSCG0415 = fftfreq(length(spectrumSCG0415), 1 / (scg20230415Time[2] - scg20230415Time[1]))
freqsSCG0329 = fftfreq(length(spectrumSCG0329), 1 / (scg20230329Time[2] - scg20230329Time[1]))


# Define Window Function with respect to specific frequencies
freqStart = -0.12
freqStop = -0.27

# Find closest index of the frequency
indexStop0403 = argmin(abs.(freqsSCG0403 .- freqStart))
indexStart0403 = argmin(abs.(freqsSCG0403 .- freqStop))

indexStop0415 = argmin(abs.(freqsSCG0415 .- freqStart))
indexStart0415 = argmin(abs.(freqsSCG0415 .- freqStop))

indexStop0329 = argmin(abs.(freqsSCG0329 .- freqStart))
indexStart0329 = argmin(abs.(freqsSCG0329 .- freqStop))


# create tukey window
tukeyWindowFreq0403 = DSP.Windows.tukey(length(freqsSCG0403[indexStart0403:indexStop0403]), 0.4)
tukeyWindowFreq0415 = DSP.Windows.tukey(length(freqsSCG0415[indexStart0415:indexStop0415]), 0.4)
tukeyWindowFreq0329 = DSP.Windows.tukey(length(freqsSCG0329[indexStart0329:indexStop0329]), 0.4)

window0403 = zeros(length(freqsSCG0403))
window0403[indexStart0403:indexStop0403] = tukeyWindowFreq0403

window0415 = zeros(length(freqsSCG0415))
window0415[indexStart0415:indexStop0415] = tukeyWindowFreq0415

window0329 = zeros(length(freqsSCG0329))
window0329[indexStart0329:indexStop0329] = tukeyWindowFreq0329




# Creating the figure and axes
fig2, ax2 = subplots(1, 1, figsize=(12, 8))

# ax2.plot(freqsSCG[1:Int(length(freqsSCG) / 2)], abs.(spectrumSCG[1:Int(length(freqsSCG) / 2)]))
# ax2.plot(freqDeg[1:Int(length(freqDeg) / 2)], abs.(spectrumDegenerate[1:Int(length(freqDeg) / 2)]))

ax2.semilogy(freqsSCG0403, abs.(spectrumSCG0403).^2,label="Measurement")
ax2.semilogy(freqsSCG0403, window0403.*abs.(spectrumSCG0403).^2,color="red",label="Filtered" )

#ax2.semilogy(freqDeg, abs.(spectrumDegenerate).^2)
 ax2.set_xlim([-0.1, -0.3])
ax2.set_xlabel("Frequency (PHz)")
ax2.set_ylabel("Spectr. Intensity (arb.u.)")
ax2.set_title("Spectrum Plot")
ax2.legend()

fig2.savefig("WaveformSpectrumSCG.png", dpi=800)
show()



# filter the spectrum and backconvert to time domain

filteredSpectrumSCG = spectrumSCG .* window.*2
filterSCG = ifft(filteredSpectrumSCG)

scg[:,1] .-= scg[findmax(abs.(filterSCG))[2],1]


# Creating the figure and axes

fig3, ax3 = subplots(1, 1, figsize=(12, 8))
ax3.plot(scg[:, 1], scg[:, 2])
ax3.plot(scg[:, 1], filterSCG)
ax3.plot(scg[:, 1], filterSCG.*DSP.Windows.tukey(length(filterSCG), 0.4))
ax3.set_xlabel("Time (fs)")
ax3.set_ylabel("Field (arb.u.)")
ax3.set_title("Filtered Waveform Plot")


fig3.savefig("FilteredWaveformSCG.png", dpi=800)
show()

CSV.write("FilteredWaveformSCG.csv", DataFrame(time = reverse(scg[:, 1]), field = reverse((filterSCG./maximum(abs.(filterSCG))).*DSP.Windows.tukey(length(filterSCG), 0.4))),writeheader=true)