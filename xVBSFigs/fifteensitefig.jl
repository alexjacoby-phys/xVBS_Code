using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


Analytic =open("/Users/alexjacoby/Documents/Julia_Code/xVBSFigs/Analytic_EE_15Site-FIG.txt") do file
     readdlm(file)
end
DMRG = open("/Users/alexjacoby/Documents/Julia_Code/xVBSFigs/DMRG_EE_15_SITE-FIG.txt") do file
     readdlm(file)
end
X1 = Vector{Int64}(Analytic[1,:])
X2 = Vector{Int64}(DMRG[1,:])
Y1 = Analytic[2,:]
Y2 = DMRG[2,:]
plot(X2,Y2, xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,color = :blue,linewidth = 8,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Computer Modern", label = "DMRG", legend = :right)
plot!(X1,Y1,  linewidth = 12,color = :red,label = "Analytic",seriestype = :scatter,legend = :right)

savefig("NumericalvsExact.pdf")
#savefig("NumericalvsExact.png")
