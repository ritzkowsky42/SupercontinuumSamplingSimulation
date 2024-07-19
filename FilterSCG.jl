using FFTW

using PyPlot
using CSV
using DataFrames
using DSP
using Tables

# Laden der CSV-Datei
scg = CSV.read("Waveforms/ShortGateSCG.csv", DataFrame)
degen = CSV.read("Waveforms/Degenerate1690.csv", DataFrame)

# # Erstellen der Figur und der Achsen
# fig, ax = subplots(figsize=(12, 8))  # Größe in Zoll

# # Hinzufügen der Datenlinien zum Plot
# ax.plot(scg[:, 1], scg[:, 2])
# ax.plot(degen[:, 1], degen[:, 2])

# # Beschriftungen hinzufügen
# ax.set_xlabel("Time (fs)")
# ax.set_ylabel("Field (arb.u.)")
# ax.set_title("Waveform Plot")

# Anzeigen des Plots


# Converting time domain to Frequency domain

spectrumSCG = fft(scg[:, 2])
spectrumDegenerate = fft(degen[:, 2])

# Creating the frequency axis
freqsSCG = fftfreq(length(scg[:, 1]), 1 / (scg[2, 1] - scg[1, 1]))
freqDeg = fftfreq(length(degen[:, 1]), 1 / (degen[2, 1] - degen[1, 1]))


# Define Window Function with respect to specific frequencies
freqStart = 0.12
freqStop = 0.27

# Find closest index of the frequency
indexStop = argmin(abs.(freqsSCG .- freqStart))
indexStart = argmin(abs.(freqsSCG .- freqStop))




# create tukey window
tukeyWindow = DSP.Windows.tukey(length(freqsSCG[indexStart:indexStop]), 0.4)

window = zeros(length(freqsSCG))
window[indexStart:indexStop] = tukeyWindow


indexStart
indexStop

# Creating the figure and axes
fig2, ax2 = subplots(1, 1, figsize=(12, 8))

# ax2.plot(freqsSCG[1:Int(length(freqsSCG) / 2)], abs.(spectrumSCG[1:Int(length(freqsSCG) / 2)]))
# ax2.plot(freqDeg[1:Int(length(freqDeg) / 2)], abs.(spectrumDegenerate[1:Int(length(freqDeg) / 2)]))

ax2.semilogy(freqsSCG, abs.(spectrumSCG).^2)
ax2.semilogy(freqsSCG, window.*abs.(spectrumSCG).^2,color="red" )

#ax2.semilogy(freqDeg, abs.(spectrumDegenerate).^2)
 ax2.set_xlim([0, 8])
ax2.set_xlabel("Frequency (PHz)")
ax2.set_ylabel("Spectr. Intensity (arb.u.)")
ax2.set_title("Spectrum Plot")

fig2.savefig("WaveformSpectrum.png", dpi=800)
show()

# filter the spectrum and backconvert to time domain

filteredSpectrumSCG = spectrumSCG .* window.*2
filterSCG = ifft(filteredSpectrumSCG)


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

CSV.write("FilteredWaveformSCG.csv", DataFrame(time = scg[:, 1], field = filterSCG.*DSP.Windows.tukey(length(filterSCG), 0.4)), writeheader = false)