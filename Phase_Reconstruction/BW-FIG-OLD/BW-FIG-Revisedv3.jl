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
    push!(GRAD,RGB(0.7 -0.65(i/K),0.7 -0.65(i/K),0.7 -0.65(i/K)))
end


X = [A[1,:],B[1,:],C[1,:],D[1,:]]
Y = [A[2,:],B[2,:],C[2,:],D[2,:]]
Labels = [L"\theta_{0} = 0.51, \;  \beta = 0.53 ",L"\theta_{0} = 0.59, \; \beta = 1.44 \; ",L"\theta_{0}  = 0.63, \;\beta = 0.84",L"\theta_{0}  = 0.73, \; \beta = 0.71 "]

plot(X[1],Y[1], label = false, color = GRAD[1],linewidth = 2)
plot!(X[1],Y[1], xlabel = L"\theta",ylabel = L"{\rm Relative \ Entropy,} \ \mathcal{S}_{\left.\rho\right|\sigma^{{\rm BW}}}" ,seriestype = :scatter,color = GRAD[1] ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 20,yguidefontsize =15, fontfamily = "Times", label = Labels[1], legendfontsize = 15, legend = :right, markersize = 6)

for i in 2:4
    plot!(X[i],Y[i], label = false, color = GRAD[i],linewidth = 2)
    plot!(X[i],Y[i],color = GRAD[i],label = Labels[i],seriestype = :scatter,legend = :topright, markersize = 6)
end
plot!()
savefig("BWANSATZrevisedv3.pdf")
