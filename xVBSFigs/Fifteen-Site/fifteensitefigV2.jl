using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


Analytic =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fifteen-Site/Analytic_EE_15Site-FIG.txt") do file
     readdlm(file)
end
DMRG = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Fifteen-Site/DMRG_EE_15_SITE-FIG.txt") do file
     readdlm(file)
end
X1 = Vector{Int64}(Analytic[1,:])
X2 = Vector{Int64}(DMRG[1,:])
Y1 = Analytic[2,:]
Y2 = DMRG[2,:]
plot(X2,Y2, xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,seriestype = :scatter,color = :black ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Times", label = "DMRG", legend = :right, marker = :xcross, markersize = 8)
plot!(X1,Y1,  linewidth = 12,color = :black,label = "Analytic",seriestype = :scatter, markeralpha = 0.5,legend = :right, markersize = 6)

savefig("NumericalvsExact.pdf")
#savefig("NumericalvsExact.png")
