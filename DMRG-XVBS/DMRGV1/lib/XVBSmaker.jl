module XVBSmaker
import LinearAlgebra
using ..tensor
using ..MPutil
using ..Opt

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()




  function XVBSmake(Ns::Integer,theta::Float64)
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

    return Hamiltonian
  end

  function makePsi0(Ns::Integer)
    hereQS = convert(Int64,6)
    QS = cld(hereQS,2)
    initTensor = [zeros(1,hereQS,1) for i=1:Ns]
    for i = 1:Ns
       initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
    end
    psi = MPS(initTensor,1)
    return psi
  end



  export  XVBSmake, makePsi0
end
