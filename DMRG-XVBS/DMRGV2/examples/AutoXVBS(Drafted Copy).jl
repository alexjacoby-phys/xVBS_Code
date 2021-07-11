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



println("#############")
println("nonQN version")
println("#############")
theta = π

print("thing")
a = cos(theta) #bilinear
b = 0.25*sin(theta) #biquadratic
FirstSiteMaker = zeros(Float64,1,6,6,2+length(Generators)+length(BiQuadLeft))
LastSiteMaker = zeros(Float64,2+length(Generators)+length(BiQuadLeft),6,6,1)
MiddleSiteMaker = zeros(Float64, 2+length(Generators)+length(BiQuadLeft),6,6,2+length(Generators)+length(BiQuadLeft))


FirstSiteMaker[1,:,:,2+length(Generators)+length(BiQuadLeft)] = Id
LastSiteMaker[1,:,:,1] = Id
MiddleSiteMaker[1,:,:,1] = Id
MiddleSiteMaker[2+length(Generators)+length(BiQuadLeft),:,:,2+length(Generators)+length(BiQuadLeft)] = Id
for i in 1:length(Generators)
  FirstSiteMaker[1,:,:,i+1] = a*Generators[i]
  MiddleSiteMaker[i+1,:,:,1] = transpose(Generators[i])
  MiddleSiteMaker[2+length(Generators)+length(BiQuadLeft),:,:,i+1] = a*Generators[i]
  LastSiteMaker[i+1,:,:,1] = transpose(Generators[i])
end
for i in 1:length(BiQuadLeft)
  FirstSiteMaker[1,:,:,i+1+length(Generators)] = b*BiQuadLeft[i]
  MiddleSiteMaker[1,:,:,i+1+length(Generators)] = BiQuadRight[i]
  MiddleSiteMaker[2+length(Generators)+length(BiQuadLeft),:,:,i+1+length(Generators)] = b*BiQuadLeft[i]
  LastSiteMaker[i+1+length(Generators),:,:,1] = BiQuadRight[i]
end

 Z1 = tens(FirstSiteMaker)
 Z2 = tens(MiddleSiteMaker)
 Z3 = tens(LastSiteMaker)


localopsarray = Vector{tens{Float64}}(undef, Ns)
localopsarray[1] = Z1
localopsarray[Ns] = Z3
for i in 1:Ns-2
  localopsarray[i+1] = Z2
end

Hamiltonian = MPO(localopsarray)

@time params = dmrg(psi,Hamiltonian,maxm=999,sweeps=25,cutoff=1E-9)
