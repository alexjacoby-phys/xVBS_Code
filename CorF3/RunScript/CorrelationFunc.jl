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
Pkg.add("Combinatorics")
using Combinatorics
using LinearAlgebra
using DelimitedFiles
path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia
for i in 0:10
    θ = 0.5+0.02*i
    Ns = 10
    MXM = 40 #6^Ns
    psi = makePsi0(Ns)
    H = XVBSmake(Int(Ns),#=Float64(1/2)=#cos(θ),#=Float64(1/12)=#0.25*sin(θ))
    #@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
    nD = 100 # number of elements of sqrt(rho) to save
    dummy_D = zeros(nD)
    @time params = dmrg(psi,H,maxm=MXM,sweeps=100,cutoff=1E-6, storeD=dummy_D,allSvNbond=true)

    bilin = XVBSmake(Int(Ns),1.,0.)
    biquad = XVBSmake(Int(Ns),0.,1.)

    bldis = expect(psi,bilin)
    bqdis = expect(psi,biquad)
    A = expect(psi,bilin,bilin)- bldis^2
    B = expect(psi,biquad,biquad)-bqdis^2
    C = expect(psi,bilin,biquad)-bldis*bqdis
    savearray = [A,B,C]
    filename = string("corrmat-theta=",θ,".txt")
    touch(filename)
    open(filename, "w") do io
        writedlm(filename, savearray)
    end
end















# function linear(m::Int64,n::Int64,Ns::Int64,psi)
#     Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
#     base = Array{Float64,2}[Id for i = 1:Ns]
#     lin = zeros(Float64,16)
#     lintp = zeros(Float64,16)
#     for i in 1:16
#         lin[i]= expect(psi,mpoterm(1.0,[Generators[i],Id],[m,n],base))
#     end
#     for i in 1:16
#         lintp[i] = expect(psi,mpoterm(1.0,[Id,transpose(Generators[i])],[m,n],base))
#     end
#     bilin = zeros(Float64,16)
#     for i in 1:16
#         bilin[i] = expect(psi,mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base))
#     end
#     return sum(bilin)-dot(lin,lintp)
# end
#
#
#
# function bilinear(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
#     Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
#     base = Array{Float64,2}[Id for i = 1:Ns]
#     bilinop1  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[m,n],base)
#     for i in 2:16
#         bilinop1 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base)
#     end
#     bilinop2  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[k,l],base)
#     for i in 2:16
#         bilinop2 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[k,l],base)
#     end
#     return expect(psi,bilinop1,bilinop2)-expect(psi,bilinop1)*expect(psi,bilinop2)
# end
#
# function biquadratic(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
#     Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
#     base = Array{Float64,2}[Id for i = 1:Ns]
#     biquadop1 = mpoterm(1.0,[BiQuadLeft[1],BiQuadRight[1]],[m,n],base)
#     for i in 2:196
#         biquadop1 += mpoterm(1.0,[BiQuadLeft[i],BiQuadRight[i]],[m,n],base)
#     end
#     biquadop2 = mpoterm(1.0,[BiQuadLeft[1],BiQuadRight[1]],[k,l],base)
#     for i in 2:196
#         biquadop2 += mpoterm(1.0,[BiQuadLeft[i],BiQuadRight[i]],[k,l],base)
#     end
#     return expect(psi,biquadop1,biquadop2)-expect(psi,biquadop1)*expect(psi,biquadop2)
# end
#
#
#
# function BQBLIN(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
#     Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
#     base = Array{Float64,2}[Id for i = 1:Ns]
#     bilinop  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[m,n],base)
#     for i in 2:16
#         bilinop += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base)
#     end
#     biquadop = mpoterm(1.0,[BiQuadLeft[1],BiQuadRight[1]],[k,l],base)
#     for i in 2:196
#         biquadop += mpoterm(1.0,[BiQuadLeft[i],BiQuadRight[i]],[k,l],base)
#     end
#     return expect(psi,bilinop,biquadop)-expect(psi,bilinop)*expect(psi,biquadop)
# end
#
#
#
# l=5
# A = biquadratic(l,l+1,l,l+1,Ns,psi)
# B = bilinear(l,l+1,l,l+1,Ns,psi)
# C = BQBLIN(l,l+1,l,l+1,Ns,psi)
# expect(psi,H,H)-expect(psi,H)^2



#plot(2:Ns,YAX,title = string("SvN of xVBS Chain of Length ",Ns), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig(string(Ns,"_site"))
