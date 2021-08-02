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
function SvN(dmat::Array{Complex{Float64},2},cutoff::Float64)
    schmidtnumbers = filter(!iszero,eigvals(dmat))
    print(schmidtnumbers)
    vne = 0
    for i in 1:length(schmidtnumbers)
        if real(schmidtnumbers[i]) > cutoff
            vne = vne - real(schmidtnumbers[i])*log(real(schmidtnumbers[i]))
        end
    end
    return vne
end
function traceop(k::Int64,n::Int64,L::Int64#=k is which state, n is the site L self explanatory=#)
    Id = sparse([1. 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
    vec = sparse(zeros(Complex{Float64},6))
    vec[k] = 1.
    if n == 1
        trop = vec
    else
        trop = Id
    end
    for i in 2:L
        if i != n
            trop = kron(trop,Id)
        else
            trop = kron(trop,transpose(vec))
        end
    end
    return sparse(trop)
end

function tracefunc(rho#=::Array{Complex{Float64},2}=#,n::Int64,L::Int64)
    tracerho = traceop(1,n,L)*rho*transpose(traceop(1,n,L))+traceop(2,n,L)*rho*transpose(traceop(2,n,L))+traceop(3,n,L)*rho*transpose(traceop(3,n,L))+traceop(4,n,L)*rho*transpose(traceop(4,n,L))+traceop(5,n,L)*rho* transpose(traceop(5,n,L))+traceop(6,n,L)*rho*transpose(traceop(6,n,L))
    return tracerho
end


function partialtrace(rho#=::Array{Complex{Float64},2}=#,n::Int64,L::Int64#= n = how many sites you want to trace out from right to left, L is size as always=#)
    result = rho
    for i in 0:n-1
        result = tracefunc(result,L-i,L-i)
    end
    return result
end

function makerho(statein::Vector{Complex{Float64}})
    rhostate = statein
    rho = Array{Complex{Float64},2}
    rho = kron(sparse(rhostate),sparse(conj(transpose(rhostate))))
    rho = (1/tr(rho))*rho
    return rho
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



L = 5
θ = 0.47
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
    Hterm = kron(K...)
    push!(Harray,Hterm)
end
BC = BoundaryOp(1,L,1)
@time H = sum(Harray, dims = 1)[1]-1000*BC
@time states = eigs(H ,which = :SR, nev =6 #=, maxiter = 1000=#)
state=states[2][:,1]
normalize!(state)
@time dmat = kron(sparse(state),sparse(adjoint(state)))
@time rho = partialtrace(dmat,3,L)

eigen(Array(rho))





YAX = zeros(L-1)

for i in 1:L-1
    YAX[i] = SvN(partialtrace(dmat,i,L),0.01)
end
XAX = 2:L
savearray = [XAX,YAX]
filename  = string("ED_EE_",L,"Site")
touch(filename)
open(filename,"w") do io
    writedlm(filename, savearray)
end


#plot(XAX,YAX,title = string("SvN of xVBS Chain of Length ",L, ", N =",n), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig("5siteED")
