using Plots
using DelimitedFiles
using LinearAlgebra

dat1 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.5.txt") do file
    readdlm(file)
end
dat2 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.52.txt") do file
    readdlm(file)
end
dat3 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.54.txt") do file
    readdlm(file)
end
dat4 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.56.txt") do file
    readdlm(file)
end
dat5 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.58.txt") do file
    readdlm(file)
end
dat6 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.6.txt") do file
    readdlm(file)
end
dat7 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.62.txt") do file
    readdlm(file)
end
dat8 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.64.txt") do file
    readdlm(file)
end
dat9 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.66.txt") do file
    readdlm(file)
end
dat10 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.6799999999999999.txt") do file
    readdlm(file)
end
dat11 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/CORMAT/Cormat-dat2/corrmat-theta=0.7.txt") do file
    readdlm(file)
end


M1 = [dat1[1] dat1[3] ; dat1[3] dat1[2]]
M2 = [dat2[1] dat2[3] ; dat2[3] dat2[2]]
M3 = [dat3[1] dat3[3] ; dat3[3] dat3[2]]
M4 = [dat4[1] dat4[3] ; dat4[3] dat4[2]]
M5 = [dat5[1] dat5[3] ; dat5[3] dat5[2]]
M6 = [dat6[1] dat6[3] ; dat6[3] dat6[2]]
M7 = [dat7[1] dat7[3] ; dat7[3] dat7[2]]
M8 = [dat8[1] dat8[3] ; dat8[3] dat8[2]]
M9 = [dat9[1] dat9[3] ; dat9[3] dat9[2]]
M10 = [dat10[1] dat10[3] ; dat10[3] dat10[2]]
M11 = [dat11[1] dat11[3] ; dat11[3] dat11[2]]


Hvec1 = eigvecs(M1)[1,:]
Hvec2 = eigvecs(M2)[1,:]
Hvec3 = eigvecs(M3)[1,:]
Hvec4 = eigvecs(M4)[1,:]
Hvec5 = eigvecs(M5)[1,:]
Hvec6 = eigvecs(M6)[1,:]
Hvec7 = eigvecs(M7)[1,:]
Hvec8 = eigvecs(M8)[1,:]
Hvec9 = eigvecs(M9)[1,:]
Hvec10 = eigvecs(M10)[1,:]
Hvec11 = eigvecs(M11)[1,:]

vecvec = [Hvec1,Hvec2,Hvec3,Hvec4,Hvec5,Hvec6,Hvec7,Hvec8,Hvec9,Hvec10,Hvec11]


thetarecon = []
for i in 1:11
    vec0 = vecvec[i]
    push!(thetarecon,atan(4*(vec0[2]/vec0[1])))
end

theta = [0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7]

plot(theta,theta, color = :red, linewidth = 6,xlabel = L"\theta",ylabel = L"\theta, \ {\rm Reconstructed}",label = "Exact",fontfamily = "Times",legend = :right)

plot!(theta,thetarecon,label = "Reconstructed",marker = :cross,markercolor = :black,markersize =8, seriestype = :scatter)

savefig("Thetarecon.pdf")
