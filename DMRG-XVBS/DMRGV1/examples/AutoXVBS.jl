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





Ns = 4
MXM = 6^Ns

psi = makePsi0(Ns)
H = XVBSmake(Int(Ns),Float64(π))
@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
