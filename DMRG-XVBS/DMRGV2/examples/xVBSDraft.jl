path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia

Ns = 5
dimrep = 6
theta = 0.5*π #sets the position in the configuration space; refer to Nuclear Phys. B

hereQS = convert(Int64,dimrep)
QS = cld(hereQS,2)

initTensor = [zeros(1,hereQS,1) for i=1:Ns]
for i = 1:Ns
   initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
end

psi = MPS(initTensor,1) #this is the dividing site in the bipartition

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
N = zeros(Float64,6,6)
base = Array{Float64,2}[Id for i = 1:Ns]
InitMPO = mpoterm(1,[Id,Id],[1,1],base)

print(base)
function makeH()
  baseMPO = InitMPO
    for i = 1:Ns-1
      xpos,ypos = i,i+1
      #test

      #bilinear term
      #for j in 1:16
      #  baseMPO += mpoterm(cos(theta), [Generators[j],transpose(Generators[j])], [xpos,ypos], base)
      #end
      #for j in 1:256
      #  baseMPO += mpoterm(0.25*sin(theta),[BiQuadLeft[j]*BiQuadRight[j]],[xpos,ypos],base)
      #end
    end
    return baseMPO
end

println("Making qMPO")
@time mpo = makeH()

println("#############")
println("nonQN version")
println("#############")

@time energy = dmrg(psi,mpo,maxm=45,sweeps=5,cutoff=1E-3)
