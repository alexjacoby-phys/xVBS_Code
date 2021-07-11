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

Ns = 3


hereQS = convert(Int64,6)
QS = cld(hereQS,2)

initTensor = [zeros(1,hereQS,1) for i=1:Ns]
for i = 1:Ns
   initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
end

psi = MPS(initTensor,1) #oc on site 1 ok here?

display(psi[3][2])

display(makePsi0(3)[1])

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
N = zeros(Float64,6,6)
base = Array{Float64,2}[Id for i = 1:Ns]
theta = 0.5*π

psi == makePsi0(3)

a =  sin(theta); #bilinear coefficient
b = 0.25*cos(theta); # biquadratic coefficient
print("initializations successful")
function makeH()
    baseMPO = mpoterm(a,[Id,Id],[1,1],base)
    # loop over rest of sites, adding the three terms of H
    for i = 1:Ns-1
    @time print("site ", i, " constructed")
        for j in 1:length(Generators)
            baseMPO += mpoterm(a,[Generators[j],transpose(Generators[j])],[i,i+1],base)
        end
        for k in 1: length(BiQuadLeft)
            baseMPO += mpoterm(b,[BiQuadLeft[k],BiQuadRight[k]],[i,i+1],base)
        end
    end
    return baseMPO
end
print("Hamiltonian MPOmaker Defined")
#@time mpo = makeH()



println("#############")
println("nonQN version")
println("#############")
y2=zeros(Ns)
dummy_D = zeros(100)

#@time params = dmrg(psi,mpo,maxm=999,sweeps=250,cutoff=1E-9)




#=function correlationfnctn(i::Integer,j::Integer)
    base = Array{Float64,2}[Id for i = 1:Ns]
    correlator = mpoterm(1,[Generators[1],transpose(Generators[1])], [i,j],base)
    for k in 2:16
        correlator += mpoterm(1,[Generators[k],transpose(Generators[k])], [i,j],base)
    end
    return correlator
end

x = 1:Ns

y = zeros(Float64, Ns)

for i in 1:Ns
    y[i] = expect(psi,correlationfnctn(1,i))-0.5
end

plot(x,y)



size(mpo[2])

 Z= []
for i in 1:16
    push!(Z , convertTensor(mpo[3][18,:,:,i+1]) == Generators[i])
end

mpo[3][1,:,:,18]


display(convertTensor(mpo[3][18,:,:,2]))

display(transpose(Generators[2]))
=#
