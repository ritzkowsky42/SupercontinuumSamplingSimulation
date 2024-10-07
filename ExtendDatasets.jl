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

# TDSE0329 = CSV.read("SampledWaveform15Vnm0329.csv", DataFrame)
TDSE0403 = CSV.read("SampledWaveform15Vnm0403.csv", DataFrame)
TDSE0415 = CSV.read("SampledWaveform15Vnm0415.csv", DataFrame)


TDSE0403Extension = CSV.read("SampledWaveform15Vnm0403Extension.csv", DataFrame)
TDSE0415Extension = CSV.read("SampledWaveform15Vnm0415Extension.csv", DataFrame)

# Skalierungsfaktoren berechnen
scaling0403 = TDSE0403[1,2] / TDSE0403Extension[end,2]
scaling0415 = TDSE0415[1,2] / TDSE0415Extension[end,2]

# Erweiterungsdaten skalieren
scaled_TDSE0403Extension = hcat(TDSE0403Extension[:,1], scaling0403 .* TDSE0403Extension[:,2])
scaled_TDSE0415Extension = hcat(TDSE0415Extension[:,1], scaling0415 .* TDSE0415Extension[:,2])


# Plotten der Daten

tdse0403Time = vcat(TDSE0403Extension[:,1], TDSE0403[:,1])
tdse0403Field = vcat(scaling0403.*TDSE0403Extension[:,2], TDSE0403[:,2])

tdse0415Time = vcat(TDSE0415Extension[:,1], TDSE0415[:,1])
tdse0415Field = vcat(scaling0415.*TDSE0415Extension[:,2], TDSE0415[:,2])


plot(tdse0403Time, tdse0403Field .+2, label="TDSE0403")
plot(tdse0415Time, tdse0415Field, label="TDSE0415") 
show()

CSV.write("SampledWaveform15Vnm0403Extension.csv", DataFrame(time = tdse0403Time, field = tdse0403Field, writeheader = true))
CSV.write("SampledWaveform15Vnm0415Extension.csv", DataFrame(time = tdse0415Time, field = tdse0415Field, writeheader = true))