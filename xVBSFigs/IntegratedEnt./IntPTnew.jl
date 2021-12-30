using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting

gr()
Theta1 =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.0.txt") do file
     readdlm(file)
end


Theta2 =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.5497787143782138.txt") do file
     readdlm(file)
end

Theta3 = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DimerizedPhase_1/CBREAK1-DATA/DMRG_EE_theta=0.6283185307179586.txt") do file
     readdlm(file)
end

Analytic =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fifteen-Site/Analytic_EE_15Site-FIG.txt") do file
     readdlm(file)
end
DMRG = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fifteen-Site/DMRG_EE_15_SITE-FIG.txt") do file
     readdlm(file)
end


GRAD = []
K=3
for i in 0:K-1
    push!(GRAD,RGB(0.7 -0.6(i/K),0.7 -0.6(i/K),0.7 -0.6(i/K)))
end


X = Vector{Int64}(Theta1[1,:])
Y = [Theta1[2,:],Theta2[2,:],Theta3[2,:]]
Labels = [L"\theta = 0",L"\theta  \approx 0.55",L"\theta  \approx 0.63"]

X1 = Vector{Int64}(Analytic[1,:])
X2 = Vector{Int64}(DMRG[1,:])
Y1 = Analytic[2,:]
Y2 = DMRG[2,:]


plot(X,Y[1], xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,color = GRAD[3] ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Times", label = Labels[1], legend = :bottom, markersize = 5, marker = :rect,linewidth = 2, legendfontsize = 10)



plot!(X,Y[2],  linewidth = 2,color = GRAD[2],label = Labels[2], marker = :hex, markersize = 6)
plot!(X,Y[3],  linewidth = 2,color = GRAD[1],label = Labels[3], marker = :star5, markersize = 8)

plot!()
plot!(X2,Y2,seriestype = :scatter,color = :black , marker = :xcross, markersize = 8,label = L"{\rm XVBS, \ DMRG}")
plot!(X1,Y1,  linewidth = 12,color = :black,label = L"{\rm XVBS, \ Analytic}",seriestype = :scatter, markersize = 6, legend = :bottomright)

savefig("IntPT.pdf")
