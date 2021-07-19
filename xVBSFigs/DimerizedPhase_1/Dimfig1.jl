using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


Theta1 =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.0.txt") do file
     readdlm(file)
end

Theta2 =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.47123889803846897.txt") do file
     readdlm(file)
end

Theta3 =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.5497787143782138.txt") do file
     readdlm(file)
end

Theta4 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.6283185307179586.txt") do file
     readdlm(file)
end



GRAD = []
K=4
for i in 0:K-1
    push!(GRAD,RGB(1. -(i/K),0.,Float64(i/K)))
end


X = Vector{Int64}(Theta1[1,:])
Y = [Theta1[2,:],Theta2[2,:],Theta3[2,:],Theta4[2,:]]
Labels = [L"\theta = 0",L"\theta  \approx 0.47",L"\theta  \approx 0.55",L"\theta  \approx 0.63"]

plot(X,Y[1], xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,#=seriestype = :scatter,=#color = GRAD[1] ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Times", label = Labels[1], legend = :right, markersize = 6, markeralpha = 0.8,linewidth = 6)


for i in 2:4
    plot!(X,Y[i],  linewidth = 6,color = GRAD[i],label = Labels[i]#=,seriestype = :scatter=#,legend = :right, markersize = 6, markeralpha = 0.8)
end
plot!()
savefig("PhaseTransition.pdf")
