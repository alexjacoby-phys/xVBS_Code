using LinearAlgebra
using DelimitedFiles
using Plots

phase = []
theta = []
for i in 1:19
    dat =  open(string("/Users/alexjacoby/Documents/xVBS_Code/Phase_Reconstruction/11-Site-xVBS_Phase/DMRG_EE_Theta=",i,".txt")) do file
         readdlm(file)
    end
    EE = dat[3,:]
    thet = dat[1]
    push!(phase,EE)
    push!(theta,thet)
end

GRAD = []
K=19
for i in 0:K-1
    push!(GRAD,RGB(1. -(i/K),0.,Float64(i/K)))
end

for i in 1:length(phase)
    phase[i] = phase[i]-(sum(phase[i])/length(phase[i]))*ones(Float64,length(phase[i]))
    normalize!(phase[i])
end

phase
agr = []
for i in 1:length(phase)
    push!(agr,(dot(phase[4],phase[i]))^2)
end
agr

theta
plot(theta, agr, seriestype = :scatter)
