module XVBSmaker
import LinearAlgebra
using ..tensor
using ..MPutil
using ..Opt

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()
############################################### OPLIST for ends of chains
  function Z1(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    FirstSite = zeros(Float64,1,6,6,1+G+BQ)
    for i in 1:G
      FirstSite[1,:,:,i] = a*Generators[i]
    end
    for i in 1:BQ
      FirstSite[1,:,:,G+i]= b* BiQuadLeft[i]
    end
    FirstSite[1,:,:,1+G+BQ] = Id
    FirstSiteDone = tens(FirstSite)
    return FirstSiteDone
  end

  function Z5()
    G = length(Generators)
    BQ = length(BiQuadLeft)
    LastSite = zeros(Float64,1+G+BQ,6,6,1)
    LastSite[1,:,:,1] = Id
    for i in 1:G
      LastSite[1+i,:,:,1] = transpose(Generators[i])
    end
    for i in 1:BQ
      LastSite[1+G+i,:,:,1]=BiQuadRight[i]
    end
    LastSiteDone = tens(LastSite)
    return LastSiteDone
  end

############################################### OPLIST for next to ends of chains
  function Z2(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    NextToFirstSite = zeros(Float64,1+G+BQ,6,6,2+G+BQ)
    NextToFirstSite[1+G+BQ,:,:,2+G+BQ] = Id
    for i in 1:G
      NextToFirstSite[i,:,:,1]  = transpose(Generators[i])
    end
    for i in 1:BQ
      NextToFirstSite[G+i,:,:,1]  = BiQuadRight[i]
    end
    for i in 1:G
        NextToFirstSite[1+G+BQ,:,:,1+i] = a*Generators[i]
    end
    for i in 1:BQ
      NextToFirstSite[1+G+BQ,:,:,i+G+1] = b*BiQuadLeft[i]
    end
    NextToFirstSiteDone = tens(NextToFirstSite)
    return NextToFirstSiteDone
  end
  function Z4(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    NextToLastSite = zeros(Float64,2+G+BQ,6,6,1+G+BQ)
    NextToLastSite[1,:,:,1] = Id
    for i in 1:G
      NextToLastSite[i+1,:,:,1] = transpose(Generators[i])
    end
    for i in 1:BQ
      NextToLastSite[i+G+1,:,:,1] = BiQuadRight[i]
    end
    for i in 1:G
      NextToLastSite[G+BQ+2,:,:,i+1] = a*Generators[i]
    end
    for i in 1:BQ
      NextToLastSite[G+BQ+2,:,:,i+G+1] = b*BiQuadLeft[i]
    end
    NextToLastSiteDone = tens(NextToLastSite)
    return NextToLastSiteDone
  end
############################################### OPLIST for middle sites
  function Z3(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    MiddleSite = zeros(Float64,G+BQ+2,6,6,G+BQ+2)
    MiddleSite[1,:,:,1] = Id
    MiddleSite[G+BQ+2,:,:,G+BQ+2] = Id
    for i in 1:G
      MiddleSite[i+1,:,:,1] = transpose(Generators[i])
    end
    for i in 1:BQ
      MiddleSite[i+1+G,:,:,1] = BiQuadRight[i]
    end
    for i in 1:G
      MiddleSite[G+BQ+2,:,:,1+i] = a*Generators[i]
    end
    for i in 1:BQ
      MiddleSite[G+BQ+2,:,:,i+G+1] = b*BiQuadLeft[i]
    end
    MiddleSiteDone = tens(MiddleSite)
    return MiddleSiteDone
  end
############################################### OPLIST for 2 site case
  function Z6(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    TwoSiteL = zeros(Float64,1,6,6,BQ+G)
    for i in 1:G
      TwoSiteL[1,:,:,i] = a*Generators[i]
    end
    for i in 1:BQ
      TwoSiteL[1,:,:,i+G] = b*BiQuadLeft[i]
    end
    TwoSiteLDone = tens(TwoSiteL)
    return TwoSiteLDone
  end
  function Z7()
    G = length(Generators)
    BQ = length(BiQuadLeft)
    TwoSiteR = zeros(Float64,G+BQ,6,6,1)
    for i in 1:G
      TwoSiteR[i,:,:,1] = transpose(Generators[i])
    end
    for i in 1:BQ
      TwoSiteR[i+G,:,:,1] = BiQuadRight[i]
    end
    TwoSiteRDone = tens(TwoSiteR)
    return TwoSiteRDone
  end
############################################### OPLIST for 3 site case
  function Z8(a::Float64,b::Float64)
    G = length(Generators)
    BQ = length(BiQuadLeft)
    ThreeSiteMid = zeros(Float64,G+BQ+1,6,6,G+BQ+1)
    for i in 1:G
      ThreeSiteMid[i,:,:,1]  = transpose(Generators[i])
    end
    for i in 1:BQ
      ThreeSiteMid[i+G,:,:,1] = BiQuadRight[i]
    end
    for i in 1:G
      ThreeSiteMid[BQ+G+1,:,:,1+i] = a*Generators[i]
    end
    for i in 1:BQ
      ThreeSiteMid[BQ+G+1,:,:,1+G+i] = b*BiQuadLeft[i]
    end
    ThreeSiteMidDone = tens(ThreeSiteMid)
    return ThreeSiteMidDone
  end
############################################### Hamiltonian Maker-- Stitches all the ops together
  function XVBSmake(Ns::Integer, a::Float64, b::Float64)
     localopsarray = Vector{tens{Float64}}(undef, Ns)
     if Ns == 1
       return print("This is an error message! You need more than one site friend!")
     elseif Ns >= 5
       localopsarray[1] = Z1(a,b)
       localopsarray[Ns] = Z5()
       ###############################################
       localopsarray[2] = Z2(a,b)
       localopsarray[Ns-1] = Z4(a,b)
       ###############################################
       for i in 3:Ns-2
         localopsarray[i] = Z3(a,b)
       end
     elseif Ns == 4
       localopsarray[1] = Z1(a,b)
       localopsarray[2] = Z2(a,b)
       localopsarray[3] = Z4(a,b)
       localopsarray[4] = Z5()
     elseif Ns == 2
       localopsarray[1] = Z6(a,b)
       localopsarray[2] = Z7()
     elseif Ns == 3
       localopsarray[1] = Z1(a,b)
       localopsarray[2] = Z8(a,b)
       localopsarray[3] = Z5()
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
