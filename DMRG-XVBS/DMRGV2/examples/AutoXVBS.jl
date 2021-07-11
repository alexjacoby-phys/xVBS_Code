#########################################################################
#
#  Density Matrix Renormalization Group (and other methods) in julia (DMRjulia)
#                              v0.1
#
#########################################################################
# Made by Thomas E. Baker (2018)
# See accompanying license with this program
# This code is native to the julia programming language (v1.1.0) or (v1.5)
#

path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia
using Plots


Ns = 5
MXM = 1000 #6^Ns

psi = makePsi0(Ns)
H = XVBSmake(Int(Ns),Float64(1/2),Float64(1/12))
#@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
nD = 10 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-5, storeD=dummy_D,allSvNbond=true)
YAX = params.SvNvec


plot(1:Ns,YAX,title = string("SvN of AKLT Chain of Length ",Ns), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()



#removes zero elements to fix an annoying bug
BiQuadLeftPop = copy(BiQuadLeft)
BiQuadRightPop = copy(BiQuadRight)
for i in 1:length(BiQuadLeftPop)
    if iszero(BiQuadLeftPop[1+length(BiQuadLeftPop)-i]) | iszero(BiQuadRightPop[1+length(BiQuadLeftPop)-i])
        popat!(BiQuadRight,1+length(BiQuadLeftPop)-i)
        popat!(BiQuadLeft,1+length(BiQuadLeftPop)-i)
    end
end
