using LinearAlgebra
using Combinatorics
using SparseArrays
using Arpack
#using Plots
using DelimitedFiles

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
#creates list of operators to pass to makeop#
function MakeOpList(A::Array{Complex{Float64},2},B::Array{Complex{Float64},2},x::Int64, y::Int64,L::Int64)
    Id = [1. 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1]
    OpList = Vector{Array{Complex{Float64},2}}(undef,L)
    for i in 1:L
        if i == x
            OpList[i]= A
        elseif i == y
            OpList[i] = B
        else
            OpList[i] = Id
        end
    end
    return OpList
end
#does tensor product to make operator act on the whole space#
function MakeOp(C::Vector{Array{Complex{Float64},2}})
    HTerm = C[1]
    for i in 2:length(C)
        HTerm = kron(HTerm,C[i])
    end
    return HTerm
end
function BoundaryOp(N::Int64,L::Int64,P::Int64)
    Basis = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    function delta(a::Int64,b::Int64)
        if a == b
            return 1
        else
            return 0
        end
    end
    Id = zeros(Float64,6,6)+I
    BCOps = [Id,Id,Id,Id,Id,Id]
    asm = [zeros(Float64,6,6),zeros(Float64,6,6),zeros(Float64,6,6),zeros(Float64,6,6),zeros(Float64,6,6),zeros(Float64,6,6)]

    for i in 1:6
        for j in 1:6, k in 1:6
            μ = Basis[j][1]
            ν = Basis[j][2]
            γ = Basis[k][1]
            η = Basis[k][2]
            asm[i][j,k] = delta(i,i)*delta(μ,γ)*delta(ν,η)+delta(ν,i)*delta(i,γ)*delta(μ,η)+delta(μ,i)*delta(ν,γ)*delta(i,η)-delta(i,i)*delta(ν,γ)*delta(μ,η)-delta(ν,i)*delta(μ,γ)*delta(i,η)-delta(μ,i)*delta(i,γ)*delta(ν,η)
        end
    end
    Ids = sparse(Id)
    for i in 1:6
        BCOps[i] = BCOps[i]-asm[i]
    end
    kronarray = []
    for i in 1:L
        if i == N
            push!(kronarray,sparse(BCOps[P]))
        else
            push!(kronarray,Ids)
        end
    end
    BCterm = kron(kronarray...)
    return BCterm
end



L = 7
θ = 0.58
CM = zeros(Complex{Float64},36,36)
@time for j in 1:length(Generators)
    CM = CM + MakeOp(MakeOpList((cos(θ)+0*im)*Generators[j],(1. +0*im)*transpose(Generators[j]),1,2,2))
end
@time for j in 1:length(BiQuadLeft)
    CM = CM + MakeOp(MakeOpList((0.25*sin(θ)+0*im)*BiQuadLeft[j],(1+0*im)*BiQuadRight[j],1,2,2))
end
Ids = sparse(Id)
CMs = real(sparse(CM))

Harray = []
@time for i in 1:L-1
    K = []
    for j in 1:L-1
        if i == j
            push!(K,CMs)
        elseif j != i
            push!(K,Ids)
        end
    end
    Hterm = i*kron(K...)
    push!(Harray,Hterm)
end
BC = BoundaryOp(1,L,1)
@time H = sum(Harray, dims = 1)[1]-100*BC
@time states = eigs(H ,which = :SR, nev =10 , maxiter = 1000)
ev = real(states[1])-real(states[1])[1]*ones(10)
diagm(-2π*ev)
notrho = exp.(-ev)/sum(exp.(-ev))




rhoplus = open("Phase_Reconstruction/XVBSESPEC/DMRG_EE_ThetaIt=16_Sites=15.txt") do file
    return readdlm(file, skipstart = 3)
end


rho = Vector{Float64}(rhoplus[8,:])/sum(Vector{Float64}(rhoplus[8,:]))
