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
using Interpolations

SCG = CSV.read("FilteredWaveformSCG.csv", DataFrame,types=Complex{Float64})
deg = CSV.read("FilteredWaveformDeg.csv", DataFrame,types=Complex{Float64})

# The files are sorted in the wrong order so we need to sort them

# Center both pulses at zeros
# SCG[:,1] .-= SCG[findmax(abs.(SCG[:,2]))[2],1]
# deg[:,1] .-= deg[findmax(abs.(deg[:,2]))[2],1]

plot(SCG[:,1], real.(SCG[:,2])./maximum(abs.(SCG[:,2])), label="SCG")
plot(deg[:,1], real.(deg[:,2])./maximum(abs.(deg[:,2])), label="Deg")
show()

scgField = (SCG[:,2])
degField = (deg[:,2])

interpSCG= extrapolate(interpolate((real.(SCG[:,1]),),real.(SCG[:,2]), Gridded(Linear())),0 ) #extrapolate((interpolate(real.(scgField), BSpline(Linear()))), 0.0)

interpDeg = extrapolate(interpolate((real.(deg[:,1]),),real.(deg[:,2]), Gridded(Linear())),0 ) #extrapolate((interpolate(real.(degField), BSpline(Linear()))), 0.0)

# dtSCG = abs(SCG[2,1] - SCG[1,1])
# dtDEG = abs(deg[2,1] - deg[1,1])

# timearray = collect(range(0,length(scgField)-1,step=1))

delaySCG = interpSCG(real.(SCG[:,1]).+ 100)

plot(SCG[:,1], real.(scgField), label="SCG")
plot(SCG[:,1], delaySCG, label="SCG Delayed")
show()



t= collect(LinRange(-300/t0,300/t0,50000)) # Time domain (not used

DegRes = interpDeg(t*t0)
plot(t*t0, DegRes, label="Degenerate")
show()