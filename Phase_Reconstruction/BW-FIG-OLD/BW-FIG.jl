using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


A =open("/Users/alexjacoby/Documents/xVBS_Code/Phase_Reconstruction/BW-FIG-OLD/BW-theta=0.51-beta=0.53.txt") do file
     readdlm(file)
end

B =open("/Users/alexjacoby/Documents/xVBS_Code/Phase_Reconstruction/BW-FIG-OLD/BW-theta=0.59-beta=1.44.txt") do file
     readdlm(file)
end

C =open("/Users/alexjacoby/Documents/xVBS_Code/Phase_Reconstruction/BW-FIG-OLD/BW-theta=0.63-beta=0.84.txt") do file
     readdlm(file)
end

D = open("/Users/alexjacoby/Documents/xVBS_Code/Phase_Reconstruction/BW-FIG-OLD/BW-theta=0.73-beta=0.71.txt") do file
     readdlm(file)
end



GRAD = []
K=4
for i in 0:K-1
    push!(GRAD,RGB(1. -(i/K),0.,Float64(i/K)))
end


X = [A[1,:],B[1,:],C[1,:],D[1,:]]
Y = [A[2,:],B[2,:],C[2,:],D[2,:]]
Labels = [L"\theta = 0.51, \;  \beta = 0.53 ",L"\theta = 0.59, \beta = 1.44 \; ",L"\theta  = 0.63, \;\beta = 0.84",L"\theta  = 0.73, \; \beta = 0.71 "]

plot(X[1],Y[1], xlabel = L"\theta",ylabel = L"{\rm Relative \ Entropy,} \ \mathcal{S}_{R}" ,#=seriestype = :scatter,=#color = GRAD[1] ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Times", label = Labels[1], legend = :right, markersize = 6, markeralpha = 0.8,linewidth = 6)


for i in 2:4
    plot!(X[i],Y[i],  linewidth = 6,color = GRAD[i],label = Labels[i]#=,seriestype = :scatter=#,legend = :right, markersize = 6, markeralpha = 0.8)
end
plot!()
savefig("BWANSATZ.pdf")
