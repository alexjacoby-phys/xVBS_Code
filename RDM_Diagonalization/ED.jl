using LinearAlgebra
using Combinatorics
using Plots
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


#=function makerho(statein::Array{Complex{Float64},2},n::Int64)
    rhostate = statein[n,:]
    rho = Array{Complex{Float64},2}
    rho = kron(rhostate,transpose(rhostate))
    rho = (1/tr(rho))*rho
    return rho
end=#




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

function makerho(statein::Array{Complex{Float64},2},n::Int64)
    rhostate = statein[:,n]
    rho = Array{Complex{Float64},2}
    rho = kron(rhostate,conj(transpose(rhostate)))
    rho = (1/tr(rho))*rho
    return rho
end

function traceop(k::Int64,n::Int64,L::Int64#=k is which state, n is the site L self explanatory=#)
    vec = zeros(Complex{Float64},6)
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
    return trop
end

function tracefunc(rho::Array{Complex{Float64},2},n::Int64,L::Int64)
    tracerho = traceop(1,n,L)*rho*transpose(traceop(1,n,L))+traceop(2,n,L)*rho*transpose(traceop(2,n,L))+traceop(3,n,L)*rho*transpose(traceop(3,n,L))+traceop(4,n,L)*rho*transpose(traceop(4,n,L))+traceop(5,n,L)*rho* transpose(traceop(5,n,L))+traceop(6,n,L)*rho*transpose(traceop(6,n,L))
    return tracerho
end


function partialtrace(rho::Array{Complex{Float64},2},n::Int64,L::Int64#= n = how many sites you want to trace out from right to left, L is size as always=#)
    result = rho
    for i in 0:n-1
        result = tracefunc(result,L-i,L-i)
    end
    return result
end


for i in 2:6
    L=i
    H = zeros(6^L,6^L)
    function bcop(i::Int64)
        returnval = zeros(Complex{Float64},6,6)
        returnval[i,i] = 1.
        return returnval
    end
    projoplist = (1+0*im)*fill(Id,L)
    projoplist[1] = -10*bcop(1)
    projoplist[L]= -10*bcop(2)
    for i in 1:L-1
        for j in 1:length(Generators)
            H = H + MakeOp(MakeOpList((1/2+0*im)*Generators[j],(1. +0*im)*transpose(Generators[j]),i,i+1,L))
        end
        for j in 1:length(BiQuadLeft)
            H = H + MakeOp(MakeOpList((1/12 +0*im)*BiQuadLeft[j],(1+0*im)*BiQuadRight[j],i,i+1,L))
        end
        H = H  +100* MakeOp(MakeOpList(bcop(2),bcop(2),1,1,L))+100* MakeOp(MakeOpList(bcop(1),bcop(1),1,1,L))+100*MakeOp(MakeOpList(bcop(1),bcop(1),L,L,L))+100*MakeOp(MakeOpList(bcop(3),bcop(3),L,L,L))
    end
    states = eigvecs(H)
    #which excited state-- 1-3^L are permitted
    n=1
    dot(states[:,n],H*states[:,n])
    dmat = makerho(states,n)
    rho = partialtrace(dmat,1,L)
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
end


#plot(XAX,YAX,title = string("SvN of xVBS Chain of Length ",L, ", N =",n), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig("5siteED")
