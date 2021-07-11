

import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")

using LinearAlgebra
using DelimitedFiles
path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia


θ = 0
Ns = 3
MXM = 2 #6^Ns
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
