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

Ns = 10
spinmag = 0.5

hereQS = convert(Int64,2*spinmag+1)
QS = cld(hereQS,2)

initTensor = [zeros(1,hereQS,1) for i=1:Ns]
for i = 1:Ns
   initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
end

psi = MPS(initTensor,1) #oc on site 1 ok here?

Sx,Sy,Sz,Sp,Sm,O,Id = spinOps(s=spinmag)
function H(i::Int64)
    return [Id O O O O;
        -Sp O O O O;
        -Sm O O O O;
        -Sz O O O O;
        O Sm/2 Sp/2 Sz Id]
    end

println("Making qMPO")
@time mpo = convert2MPO(H,hereQS,Ns)

println("#############")
println("nonQN version")
println("#############")

@time energy = dmrg(psi,mpo,maxm=45,sweeps=5,cutoff=1E-9)
