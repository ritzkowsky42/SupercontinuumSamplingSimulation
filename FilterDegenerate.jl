using FFTW

using PyPlot
using CSV
using DataFrames

# Laden der CSV-Datei
scg = CSV.read("Waveforms/ShortGateSCG.csv", DataFrame)
degen = CSV.read("Waveforms/Degenerate1690.csv", DataFrame)

# Erstellen der Figur und der Achsen
fig, ax = subplots(figsize=(12, 8))  # Größe in Zoll

# Hinzufügen der Datenlinien zum Plot
ax.plot(scg[:, 1], scg[:, 2])
ax.plot(degen[:, 1], degen[:, 2])

# Beschriftungen hinzufügen
ax.set_xlabel("Time (fs)")
ax.set_ylabel("Field (arb.u.)")
ax.set_title("Waveform Plot")

# Anzeigen des Plots


# Converting time domain to Frequency domain

spectrumSCG = fft(scg[:, 2])
spectrumDegenerate = fft(degen[:, 2])

# Creating the frequency axis
freqsSCG = fftfreq(length(scg[:, 1]), 1 / (scg[2, 1] - scg[1, 1]))
freqDeg = fftfreq(length(degen[:, 1]), 1 / (degen[2, 1] - degen[1, 1]))

# Creating the figure and axes
fig2, ax2 = subplots(1, 1, figsize=(12, 8))

# ax2.plot(freqsSCG[1:Int(length(freqsSCG) / 2)], abs.(spectrumSCG[1:Int(length(freqsSCG) / 2)]))
# ax2.plot(freqDeg[1:Int(length(freqDeg) / 2)], abs.(spectrumDegenerate[1:Int(length(freqDeg) / 2)]))

ax2.semilogy(freqsSCG, abs.(spectrumSCG).^2)
ax2.semilogy(freqDeg, abs.(spectrumDegenerate).^2)
ax2.set_xlim([0, 2])
ax2.set_xlabel("Frequency (PHz)")
ax2.set_ylabel("Spectr. Intensity (arb.u.)")
ax2.set_title("Spectrum Plot")

fig2.savefig("WaveformSpectrum.png", dpi=800)
show()