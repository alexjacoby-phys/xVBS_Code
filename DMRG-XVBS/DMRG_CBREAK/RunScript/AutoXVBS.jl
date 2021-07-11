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

for i in 0:10
    θ = 0.25*π*i*0.1
    Ns = 15
    MXM = 35 #6^Ns
    psi = makePsi0(Ns)
    H = XVBSmake(Int(Ns),cos(θ),0.25*sin(θ))
    #@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
    nD = 25 # number of elements of sqrt(rho) to save
    dummy_D = zeros(nD)
    @time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-3, storeD=dummy_D,allSvNbond=true)
    YAX = params.SvNvec

    popat!(YAX,1)

    savearray = [Vector(2:Ns),YAX]
    filename = string("DMRG_EE_theta=",θ,".txt")
    touch(filename)
    open(filename,"w") do io
        writedlm(filename, savearray)
    end
end

#plot(2:Ns,YAX,title = string("SvN of xVBS Chain of Length ",Ns), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig(string(Ns,"_site"))
