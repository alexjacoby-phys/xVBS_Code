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
import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")

using LinearAlgebra
using DelimitedFiles

path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia

for i in 1:60
    θ = 0.41+(1/60)*i*1.2
    Ns = 15
    MXM = 40 #6^Ns
    psi = makePsi0(Ns)
    H = XVBSmake(Int(Ns),cos(θ),0.25*sin(θ))
    #@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
    nD = 100 # number of elements of sqrt(rho) to save
    dummy_D = zeros(nD)
    @time params = dmrg(psi,H,maxm=MXM,sweeps=100,cutoff=1E-8, storeD=dummy_D,allSvNbond=true)
    Y1 = params.SvNvec
    Y2 = params.Dvec
    popfirst!(Y1)
    popfirst!(Y2)
    savearray = [θ,Vector(2:Ns),Y1,Y2]
    filename = string("DMRG_EE_ThetaIt=",i,"_Sites=",Ns,".txt")
    touch(filename)
    open(filename,"w") do io
        writedlm(filename, savearray)
    end
end

#plot(2:Ns,YAX,title = string("SvN of xVBS Chain of Length ",Ns), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig(string(Ns,"_site"))
