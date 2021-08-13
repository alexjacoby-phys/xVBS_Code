using DelimitedFiles
using LinearAlgebra
using LaTeXStrings
using Plots

path = "../"
include(join([path,"Phase_Reconstruction/BWALGS.jl"]))
using .BWALGS

K = 15
θ = 0.59
dat = open("Phase_Reconstruction/XVBSESPEC/DMRG_EE_ThetaIt=9_Sites=15.txt") do file
    return readdlm(file, skipstart = 3)
end
dm = normalize(Vector{Float64}(dat[8,1:K]),1)


notdm = BW(7,θ,K,1500)
somewhatdm = makeBWRHO(notdm,1.44*2. *π)
TRdist(dm,somewhatdm,0.001)


ansatz = []
distances = []
theta = []
for i in 0:8
    thet = 0.02*i-0.08 + θ
    notdm = BW(7,thet,K,1500)
    somewhatdm = makeBWRHO(notdm,1.44*2. *π)
    append!(distances,TRdist(dm,somewhatdm,0.001))
    append!(ansatz,notdm)
    append!(theta,thet)
end

distances

savearray = [Array{Float64}(theta), Array{Float64}(distances)]
filename = "BW-theta=0.58-beta=1.44.txt"
touch(filename)
open(filename) do io
    writedlm(filename,savearray)
end

plot(theta,distances, xlabel = L"\theta",ylabel = "Relative Entropy" ,color = :blue,linewidth = 8,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Computer Modern", label = L"\theta = 0.73")
savefig("DMRGvsEDvsAnalytic_fivesite.pdf")
#savefig("NumericalvsExact.png")
