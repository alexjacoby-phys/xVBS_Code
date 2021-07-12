using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


ED =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/ED_EE_5Site-FIG.txt") do file
     readdlm(file)
end

DMRG = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DMRG_EE_5_SITE-FIG.txt") do file
     readdlm(file)
end

Analytic = open("/users/alexjacoby/Documents/xVBS_Code/xVBSFigs/DMRG_EE_5_SITE-FIG.txt") do file
    readdlm(file)
end

X = Vector{Int64}(ED[1,:])
Y1 = ED[2,:]
Y2 = DMRG[2,:]
Y3 = Analytic[2,:]
plot(X,Y1, xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,color = :blue,linewidth = 8,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Computer Modern", label = "ED")
plot!(X,Y2,  linewidth = 8,color = :red,label = "DMRG", linestyle = :dot)
plot!(X,Y1,  linewidth = 12,color = :green,label = "Analytic",seriestype = :scatter,legend = :right)
savefig("DMRGvsEDvsAnalytic_fivesite.pdf")
#savefig("NumericalvsExact.png")
