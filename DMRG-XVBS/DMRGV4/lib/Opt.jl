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

"""
    Module: Opt

Operator definitions for spin systems, fermions, and t-J model
"""
module Opt
import LinearAlgebra

  """
      spinOps()

  Make operators for a spin model (magnitude of size s)
    + equivalent: spinOps(s=0.5), spinOps(0.5)
  """
  function spinOps(;s=0.5)
    QS = convert(Int64,2*s+1) #number of quantum states

    O = zeros(Float64,QS,QS) #zero matrix
    Id = copy(O) + LinearAlgebra.I #identity matrix
    oz = copy(O) # z operator
    op = copy(O) # raising operator
    for (q,m) in enumerate(s:-1:-s) #counts from m to -m (all states)
      oz[q,q] = m
      if m+1 <= s
        op[q-1,q] = sqrt(s*(s+1)-m*(m+1))
      end
    end
    om = Array(op') # lowering operator
    ox = (op+om)/2 #x matrix
    oy = (op-om)/(2*im) #y matrix

    return ox,oy,oz,op,om,O,Id
  end
  function spinOps(a::Float64)
    return spinOps(s=a)
  end

  """
      fermionOps()
  Make fermion operators
  """
  function fermionOps()
    QS = 4 #fock space size
    O = zeros(Float64,QS,QS) #zero matrix
    Id = copy(O)+LinearAlgebra.I #identity

    Cup = copy(O) #annihilate (up)
    Cup[1,2] = 1.
    Cup[3,4] = 1.

    Cdn = copy(O) #annihilate (down)
    Cdn[1,3] = 1.
    Cdn[2,4] = -1.

    Nup = Cup' * Cup #density (up)
    Ndn = Cdn' * Cdn #density (down)
    Ndens = Nup + Ndn #density (up + down)

    F = copy(Id) #Jordan-Wigner string operator
    F[2,2] *= -1.
    F[3,3] *= -1.

    return Cup,Cdn,Nup,Ndn,Ndens,F,O,Id
  end

  """
      tJOps()
  Operators for a t-J model
  """
  function tJOps()
    #many of the Hubbard operators can be truncated
    Cup,Cdn,Nup,Ndn,Ndens,F,O,Id = fermionOps()
    QS = 3 #fock space size
    Cup = Cup[1:QS,1:QS]
    Cdn = Cdn[1:QS,1:QS]
    Nup = Nup[1:QS,1:QS]
    Ndn = Ndn[1:QS,1:QS]
    Ndens = Ndens[1:QS,1:QS]
    F = F[1:QS,1:QS]
    O = O[1:QS,1:QS]
    Id = Id[1:QS,1:QS]

    Sz = copy(O) #z-spin operator
    Sz[2,2] = 0.5
    Sz[3,3] = -0.5

    Sp = copy(O) #spin raising operator
    Sp[3,2] = 1.
    Sm = Sp' #spin lowering operator

   return Cup,Cdn,Nup,Ndn,Ndens,F,Sz,Sp,Sm,O,Id
  end


  function xVBSOps()
    #SU(4) generators + Relevant Operators for VBS
    QS = 6
    N = zeros(Float64, QS, QS)
    Gen = zeros(Float64,4,4,QS,QS)
    #Cartan Subalgebra (Rnk 3+1-- not linearly independent, but doesn't matter)
    Gen[1,1,1,1] = Gen[1,1,2,2] = Gen[1,1,3,3] = 0.5
    Gen[1,1,4,4] = Gen[1,1,5,5] = Gen[1,1,6,6] = -0.5
    #####
    Gen[2,2,1,1] = Gen[2,2,4,4] = Gen[2,2,5,5] = 0.5
    Gen[2,2,2,2] = Gen[2,2,3,3] = Gen[2,2,6,6] = -0.5
    #####
    Gen[3,3,1,1] = Gen[3,3,3,3] = Gen[3,3,5,5] = 0.5
    Gen[3,3,2,2] = Gen[3,3,4,4] = Gen[3,3,6,6] = -0.5
    #####
    Gen[4,4,3,3] = Gen[4,4,5,5] = Gen[4,4,6,6] = 0.5
    Gen[4,4,1,1] = Gen[4,4,2,2] = Gen[4,4,4,4] = -0.5
    #####
    #The rest of the generators
    ##### first block in c code
    Gen[1,2,4,2] = Gen[2,1,2,4] = Gen[1,2,5,3] = Gen[2,1,3,5] = 1.0
    ##### third block in c code
    Gen[1,4,5,1] = Gen[4,1,1,5] = Gen[1,4,6,2] = Gen[4,1,2,6] = -1.0
    ##### fourth block in c code
    Gen[2,3,2,1] = Gen[3,2,1,2] = Gen[2,3,6,5] = Gen[3,2,5,6] = 1.0
    ##### sixth block in c code
    Gen[3,4,3,2] = Gen[4,3,2,3] = Gen[3,4,5,4] = Gen[4,3,4,5] = 1.0
    ##### second block in c code
    Gen[1,3,4,1] = Gen[3,1,1,4] = 1.0
    Gen[1,3,6,3] = Gen[3,1,3,6] = -1.0
    ##### fifth block in c code
    Gen[2,4,3,1] = Gen[4,2,1,3] = 1.0
    Gen[2,4,6,4] = Gen[4,2,4,6] = -1.0
   ##### Manual reshape of generators-- reshape function didn't work
    indx = [[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4]]
    Generators = fill(N,16)
    for i = 1:16
      Generators[i] = Gen[indx[i][1],indx[i][2],:,:]
    end
    ###### Bilinear is already covered
    ###### Biquadratic
    BiQuadLeft = []
    BiQuadRight = []
    for i in 1:16, j in 1:16
      if Generators[i]*Generators[j] !== zeros(Float64,6,6)
        push!(BiQuadLeft,Generators[i]*Generators[j])
      end
      if transpose(Generators[i])*transpose(Generators[j]) !== zeros(Float64,6,6)
        push!(BiQuadRight,transpose(Generators[i])*transpose(Generators[j]))
      end
    end
    Id = copy(N)+LinearAlgebra.I
    #removes zero elements to fix an annoying bug
    BiQuadLeftPop = copy(BiQuadLeft)
    BiQuadRightPop = copy(BiQuadRight)
    for i in 1:length(BiQuadLeftPop)
        if iszero(BiQuadLeftPop[1+length(BiQuadLeftPop)-i]) | iszero(BiQuadRightPop[1+length(BiQuadLeftPop)-i])
            popat!(BiQuadRight,1+length(BiQuadLeftPop)-i)
            popat!(BiQuadLeft,1+length(BiQuadLeftPop)-i)
        end
    end
  return Gen, Id, Generators, BiQuadLeft, BiQuadRight
end

Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()

  export  spinOps,fermionOps,tJOps,xVBSOps
end
