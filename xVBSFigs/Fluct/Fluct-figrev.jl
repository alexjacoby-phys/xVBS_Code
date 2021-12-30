using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


A =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fluct/Fluctuation-Dat/Fluctuation-theta=0.51.txt") do file
     readdlm(file)
end

B =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fluct/Fluctuation-Dat/Fluctuation-theta=0.59.txt") do file
     readdlm(file)
end

C =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fluct/Fluctuation-Dat/Fluctuation-theta=0.63.txt") do file
     readdlm(file)
end

D = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fluct/Fluctuation-Dat/Fluctuation-theta=0.73.txt") do file
     readdlm(file)
end



GRAD = []
K=4
for i in 0:K-1
    push!(GRAD,RGB(0.7 -0.65*(i/K),0.7 -0.65*(i/K),0.7 -0.65*(i/K)))
end


X = [A[1,:],B[1,:],C[1,:],D[1,:]]
Y = [A[2,:],B[2,:],C[2,:],D[2,:]]
Labels = [L"\theta = 0.51",L"\theta = 0.59",L"\theta  = 0.63",L"\theta  = 0.73"]

plot(X[1],Y[1],color = GRAD[1],linewidth = 2,label = false)
plot!(X[1],Y[1], xlabel = L"\theta",ylabel = L"{\rm Fluctuation} \ \left<H^{2}\right>-\left<H\right>^{2}" ,seriestype = :scatter,color = GRAD[1] ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 20,yguidefontsize =15, fontfamily = "Times", label = Labels[1], legendfontsize = 15, legend = :right, markersize = 6)



for i in 2:4
    plot!(X[i],Y[i],  linewidth = 2,color = GRAD[i],label = false,legend = :right, markersize = 6)
    plot!(X[i],Y[i],color = GRAD[i],label = Labels[i],seriestype = :scatter,legend = :topright, markersize = 6)
end
plot!()
savefig("FluctuationFigRevised.pdf")
