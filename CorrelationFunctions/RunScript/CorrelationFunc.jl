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

using LinearAlgebra
using DelimitedFiles

path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia

#for i in 1:60
θ = 0.58
Ns = 10
MXM = 6 #6^Ns
psi = makePsi0(Ns)
H = XVBSmake(Int(Ns),cos(θ),0.25*sin(θ))
#@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
nD = 100 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=100,cutoff=1E-4, storeD=dummy_D,allSvNbond=true)
psi

function linear(m::Int64,n::Int64,Ns::Int64,psi)
    Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
    base = Array{Float64,2}[Id for i = 1:Ns]
    lin = zeros(Float64,16)
    lintp = zeros(Float64,16)
    for i in 1:16
        lin[i]= expect(psi,mpoterm(1.0,[Generators[i],Id],[m,n],base))
    end
    for i in 1:16
        lintp[i] = expect(psi,mpoterm(1.0,[Id,transpose(Generators[i])],[m,n],base))
    end
    bilin = zeros(Float64,16)
    for i in 1:16
        bilin[i] = expect(psi,mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base))
    end
    return sum(bilin)-dot(lin,lintp)
end



function bilinear(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
    Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
    base = Array{Float64,2}[Id for i = 1:Ns]
    bilinop1  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[m,n],base)
    for i in 2:16
        bilinop1 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base)
    end
    bilinop2  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[k,l],base)
    for i in 2:16
        bilinop2 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[k,l],base)
    end
    return expect(psi,bilinop*bilinop2)-expect(psi,bilinop)*expect(psi,bilinop2)
end

function biquadratic(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
    Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
    base = Array{Float64,2}[Id for i = 1:Ns]
    bilinop1  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[m,n],base)
    for i in 2:16
        bilinop1 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base)
    end
    bilinop2  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[k,l],base)
    for i in 2:16
        bilinop2 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[k,l],base)
    end
    biquad1 = bilinop1*bilinop1
    biquad2 = bilinop2*bilinop2
    return expect(psi,biquad1*biquad2)-expect(psi,biquad1)*expect(psi,biquad2)
end

function BQBLIN(m::Int64,n::Int64,k::Int64,l::Int64,Ns::Int64,psi)
    Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
    base = Array{Float64,2}[Id for i = 1:Ns]
    bilinop1  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[m,n],base)
    for i in 2:16
        bilinop1 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[m,n],base)
    end
    bilinop2  = mpoterm(1.0,[Generators[1],transpose(Generators[1])],[k,l],base)
    for i in 2:16
        bilinop2 += mpoterm(1.0,[Generators[i],transpose(Generators[i])],[k,l],base)
    end
    biquad1 = bilinop1*bilinop1
    return expect(psi,biquad1*bilinop1)-expect(psi,biquad1)*expect(psi,bilinop1)
end

biquadratic(5,6,5,6,Ns,psi)
bilinear(5,6,5,6,Ns,psi)
BQBLIN(5,6,5,6,Ns,psi)


M = [79.98 -13.3 ; -13.3 2.22]
eigvecs(M)
eigvals(M)




atan((0.164*4)/(0.986))


#plot(2:Ns,YAX,title = string("SvN of xVBS Chain of Length ",Ns), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig(string(Ns,"_site"))
