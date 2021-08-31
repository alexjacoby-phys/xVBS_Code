using Plots
using DelimitedFiles
using LinearAlgebra

dat1 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat/corrmat-theta=0.51.txt") do file
    readdlm(file)
end
dat2 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat/corrmat-theta=0.59.txt") do file
    readdlm(file)
end
dat3 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat/corrmat-theta=0.63.txt") do file
    readdlm(file)
end
dat4 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat/corrmat-theta=0.73.txt") do file
    readdlm(file)
end

M1 = [dat1[1] dat1[3] ; dat1[3] dat1[2]]
M2 = [dat2[1] dat2[3] ; dat2[3] dat2[2]]
M3 = [dat3[1] dat3[3] ; dat3[3] dat3[2]]
M4 = [dat4[1] dat4[3] ; dat4[3] dat4[2]]

Hvec1 = eigvecs(M1)[1,:]
Hvec2 = eigvecs(M2)[1,:]
Hvec3 = eigvecs(M3)[1,:]
Hvec4 = eigvecs(M4)[1,:]

thetarecon = [atan(4*(Hvec1[2]/Hvec1[1])),atan(4*(Hvec2[2]/Hvec2[1])),atan(4*(Hvec3[2]/Hvec3[1])),atan(4*(Hvec4[2]/Hvec4[1]))]
