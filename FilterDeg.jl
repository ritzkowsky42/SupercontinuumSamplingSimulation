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

spectrumDegenerate = fft(degen[:, 2])

# Creating the frequency axis
freqDeg = fftfreq(length(degen[:, 1]), 1 / (degen[2, 1] - degen[1, 1]))


# Define Window Function with respect to specific frequencies
freqStart = 0.15
freqStop = 0.2

# Find closest index of the frequency
indexStop = argmin(abs.(freqDeg .- freqStart))
indexStart = argmin(abs.(freqDeg .- freqStop))




# create tukey window
tukeyWindow = DSP.Windows.tukey(length(freqDeg[indexStart:indexStop]), 0.4)

window = zeros(length(freqDeg))
window[indexStart:indexStop] = tukeyWindow


indexStart
indexStop

# Creating the figure and axes
fig2, ax2 = subplots(1, 1, figsize=(12, 8))

# ax2.plot(freqsSCG[1:Int(length(freqsSCG) / 2)], abs.(spectrumSCG[1:Int(length(freqsSCG) / 2)]))
# ax2.plot(freqDeg[1:Int(length(freqDeg) / 2)], abs.(spectrumDegenerate[1:Int(length(freqDeg) / 2)]))

ax2.semilogy(freqDeg, abs.(spectrumDegenerate).^2)
ax2.semilogy(freqDeg, window.*abs.(spectrumDegenerate).^2,color="red" )

#ax2.semilogy(freqDeg, abs.(spectrumDegenerate).^2)
ax2.set_xlim([0, 8])
ax2.set_xlabel("Frequency (PHz)")
ax2.set_ylabel("Spectr. Intensity (arb.u.)")
ax2.set_title("Spectrum Plot")

fig2.savefig("WaveformSpectrum.png", dpi=800)
show()

# filter the spectrum and backconvert to time domain

filteredSpectrumDeg = spectrumDegenerate .* window.*2
filterDeg = ifft(filteredSpectrumDeg)

degen[:,1] .-= degen[findmax(abs.(filterDeg))[2],1]


# Creating the figure and axes

fig3, ax3 = subplots(1, 1, figsize=(12, 8))
ax3.plot(degen[:, 1], degen[:, 2])
ax3.plot(degen[:, 1], filterDeg)
ax3.plot(degen[:, 1], filterDeg.*DSP.Windows.tukey(length(filterDeg), 0.4))
ax3.set_xlabel("Time (fs)")
ax3.set_ylabel("Field (arb.u.)")
ax3.set_title("Filtered Waveform Plot")
 

fig3.savefig("FilteredWaveformDeg.png", dpi=800)
show()

CSV.write("FilteredWaveformDeg.csv", DataFrame(time = reverse(degen[:, 1]), field = reverse((filterDeg./maximum(abs.(filterDeg))).*DSP.Windows.tukey(length(filterDeg), 0.4))),writeheader=true)