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

Ns = 100
spinmag = 0.5

hereQS = convert(Int64,2*spinmag+1)
QS = cld(hereQS,2)

initTensor = [zeros(1,hereQS,1) for i=1:Ns]
for i = 1:Ns
   initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
end

psi = MPS(initTensor,1) #oc on site 1 ok here?

Sx,Sy,Sz,Sp,Sm,O,Id = spinOps(s=spinmag)

base = Array{Float64,2}[Id for i = 1:Ns]


a = 0.5; # prefactor for Sp/Sm term (1/2)
b = 1.0; # prefactor for Sz^2 term (1)

function makeH()
    baseMPO = mpoterm(a,[Sp,Sm],[1,2],base)
    baseMPO += mpoterm(a,[Sm,Sp],[1,2],base)
    baseMPO += mpoterm(b,[Sz,Sz],[1,2],base)
    # loop over rest of sites, adding the three terms of H
    for i = 2:Ns-1
        baseMPO += mpoterm(a,[Sp,Sm],[i,i+1],base)
        baseMPO += mpoterm(a,[Sm,Sp],[i,i+1],base)
        baseMPO += mpoterm(b,[Sz,Sz],[i,i+1],base)
    end
    return baseMPO
end

mpo = makeH()

println("#############")
println("nonQN version")
println("#############")

nD = 10 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,mpo,maxm=45,sweeps=20,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
params.energy
#params.SvNvec
#params.Dvec
