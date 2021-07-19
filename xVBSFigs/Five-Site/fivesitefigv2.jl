using DelimitedFiles
using Plots
using LaTeXStrings
using Formatting


ED =open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Five-Site/ED_EE_5Site-FIG.txt") do file
     readdlm(file)
end

DMRG = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Five-Site/DMRG_EE_5_SITE-FIG.txt") do file
     readdlm(file)
end

Analytic = open("/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/Five-Site/DMRG_EE_5_SITE-FIG.txt") do file
    readdlm(file)
end

X = Vector{Int64}(ED[1,:])
Y1 = ED[2,:]
Y2 = DMRG[2,:]
Y3 = Analytic[2,:]
plot(X,Y1, xlabel = "Site After Bipartition",ylabel = L"\mathbf{S}_{vN}" ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15,marker = :cross,markercolor = :black,markersize =8, seriestype = :scatter, fontfamily = "Times", label = "ED",legend = :right)
plot!(X,Y2,label = "DMRG",marker = :xcross,markercolor = :black,markersize =8, seriestype = :scatter)
plot!(X,Y3,marker = :circle,markercolor = :black,markersize =6, markeralpha = 0.5, markerstrokealpha = 1,label = "Analytic",seriestype = :scatter)
savefig("DMRGvsEDvsAnalytic_fivesite_V2.3.pdf")
#savefig("NumericalvsExact.png")
