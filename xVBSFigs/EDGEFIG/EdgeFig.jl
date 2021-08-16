using LinearAlgebra
using Combinatorics
using DelimitedFiles
using Plots

##This code is extremely ugly and not very well put together because I wrote two parts of it a few weeks apart. I would not recommend doing anything to change it because it will likely break.




##Number of singlets in chain is Chainlength all distances are measured in the number of singlets


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

function delta(a::Int64,b::Int64)
    if a == b
        return 1
    else
        return 0
    end
end




#A is long B is short



function x_nA(L::Int64)
    if  L!= 0
        coefficients = 3.
        for i in 2:L
            coefficients = 3*(3*coefficients-2)
        end
    elseif L == 0
        return 1
    end

    return coefficients
end


#=term 1 is just δ^{abc]d}_{e[fg]h}=#
function term1A(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return delta(a,c)*delta(β,ν)*delta(b,d)
end


#=term 2 is δ^{ab[cd]}_{e[fg]h}=#
function term2A(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return delta(a,c)*(2*delta(β,d)*delta(ν,b)+delta(β,ν)*delta(b,d))
end
#=term 3 is δ^{[ab]cd}_{e[fg]h}=#
function term3A(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return (2*delta(ν,a)*delta(β,c)+delta(a,c)*delta(β,ν))*delta(d,b)
end
#=term 4 is δ^{[ab]cd}_{ef[gh]}=#
function term4A(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return (delta(a,c)*delta(b,ν)*delta(d,β)+delta(a,ν)*delta(β,c)*delta(b,d)-delta(a,ν)*delta(b,c)*delta(β,d)-delta(β,c)*delta(a,d)*delta(b,ν))
end


function RDMA(a::Int64#=μ_{i}=#,b::Int64#=ν_{f}=#,c::Int64#=α_{i}=#,d::Int64#=β_{f}=#,β::Int64,ν::Int64,N::Int64#=Number of singlets in chain=#,L::Int64#=Number of singlets into chain for bipartition=#)
    c1 = x_nA(L)
    c2 = x_nA(N-L)
    RDMel = 6(c1-1)*(c2-1)*term1A(a,b,c,d,β,ν)+(c1-1)*term2A(a,b,c,d,β,ν)+(c2-1)*term3A(a,b,c,d,β,ν)+term4A(a,b,c,d,β,ν)
    return RDMel
end


function RHOA(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},N::Int64,L::Int64)
    Basis = [1,2,3,4]
    #=structured with the first two indices α and β and the next two μ and ν; final two indices are RDM indices with hilbert space of the boundary site between block and bulk=#
    rho = zeros(Float64,4,4,4,4,4,4)
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4, β in 1:4, ν in 1:4
        rho[a,b,c,d,β,ν] = LeftEdge[a]*RightEdge[b]*LeftEdge[c]*RightEdge[d]*RDMA(a,b,c,d,β,ν,N,L)
    end
    test = dropdims(sum(rho, dims = (1,2,3,4)),dims = (1,2,3,4))
    fnl = test*(1/tr(test))
    return fnl
end

BasisB = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]

function x_nB(L::Int64)
    coefficients = zeros(Float64,L)
    coefficients[1] = 3
    for i in 2:L
        coefficients[i] = 3*(3*coefficients[i-1]-2)
    end
    return coefficients
end


#=term 1 is just δ^{a[bc]d}_{e[fg]h}=#
function term1B(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return 2*delta(a,e)*delta(d,h)*(delta(b,f)*delta(c,g)-delta(b,g)*delta(c,f))
end


#=term 2 is δ^{ab[cd]}_{e[fg]h}=#
function term2B(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return delta(a,e)*(delta(b,f)*delta(c,g)*delta(d,h)+delta(b,g)*delta(d,f)*delta(c,h)-delta(b,g)*delta(c,f)*delta(d,h)-delta(b,f)*delta(d,g)*delta(c,h))
end
#=term 3 is δ^{[ab]cd}_{e[fg]h}=#
function term3B(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return (delta(a,e)*delta(b,f)*delta(c,g)+delta(b,e)*delta(a,g)*delta(c,f)-delta(b,e)*delta(a,f)*delta(c,g)-delta(a,e)*delta(b,g)*delta(c,f))*delta(d,h)
end
#=term 4 is δ^{[ab]cd}_{ef[gh]}=#
function term4B(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return (delta(a,e)*delta(b,f)-delta(b,e)*delta(a,f))*(delta(c,g)*delta(d,h)-delta(c,h)*delta(d,g))
end


function RDMB(a::Int64#=μ_{i}=#,b::Int64#=ν_{f}=#,c::Int64#=α_{i}=#,d::Int64#=β_{f}=#,α::Int64,β::Int64,μ::Int64,ν::Int64,N::Int64#=Number of singlets in chain=#,L::Int64#=Number of singlets into chain for bipartition=#)
    tempcof = x_nB(N)
    if L == 0
        c1 = 1
    elseif L > 0
        c1 = tempcof[L]
    end
    if N-L == 0
        c2  = 1
    elseif N-L > 0
        c2 = tempcof[N-L]
    end
    RDMel = (c1-1)*(c2-1)*term1B(a,α,β,b,c,μ,ν,d)+(c1-1)*term2B(a,b,α,β,c,ν,μ,d)+(c2-1)*term3B(α,β,a,b,c,ν,μ,d)+term4B(α,β,a,b,c,d,μ,ν)
    return RDMel
end


function RHOB(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},N::Int64,L::Int64)
    Basis = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    #=structured with the first two indices α and β and the next two μ and ν; final two indices are RDM indices with hilbert space of the boundary site between block and bulk=#
    rho = zeros(Float64,4,4,4,4,6,6)
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4, e in 1:6, f in 1:6
        μ = Basis[e][1]
        ν = Basis[e][2]
        α = Basis[f][1]
        β = Basis[f][2]
        rho[a,b,c,d,e,f] = LeftEdge[a]*RightEdge[b]*LeftEdge[c]*RightEdge[d]*RDMB(a,b,c,d,α,β,μ,ν,N,L)
    end
    test = dropdims(sum(rho, dims = (1,2,3,4)),dims = (1,2,3,4))
    fnl = test*(1/tr(test))
    return fnl
end


BasisC =([[1,2],[1,2]],[[1,3],[1,2]],[[1,4],[1,2]],[[2,3],[1,2]],[[2,4],[1,2]],[[3,4],[1,2]],[[1,2],[1,3]],[[1,3],[1,3]],[[1,4],[1,3]],[[2,3],[1,3]],[[2,4],[1,3]],[[3,4],[1,3]],[[1,2],[1,4]],[[1,3],[1,4]],[[1,4],[1,4]],[[2,3],[1,4]],[[2,4],[1,4]],[[3,4],[1,4]],[[1,2],[2,3]],[[1,3],[2,3]],[[1,4],[2,3]],[[2,3],[2,3]],[[2,4],[2,3]],[[3,4],[2,3]],[[1,2],[2,4]],[[1,3],[2,4]],[[1,4],[2,4]],[[2,3],[2,4]],[[2,4],[2,4]],[[3,4],[2,4]],[[1,2],[3,4]],[[1,3],[3,4]],[[1,4],[3,4]],[[2,3],[3,4]],[[2,4],[3,4]],[[3,4],[3,4]])

function x_nc(L::Int64)
    if  L!= 0
        coefficients = 3.
        for i in 2:L
            coefficients = 3*(3*coefficients-2)
        end
    elseif L == 0
        return 1
    end

    return coefficients
end


function delta3(μ::Int64,ν::Int64,λ::Int64,α::Int64,β::Int64,ρ::Int64)
    return delta(μ,α)*delta(ν,β)*delta(λ,ρ)+delta(λ,α)*delta(μ,β)*delta(ν,ρ)+delta(ν,α)*delta(λ,β)*delta(μ,ρ)-delta(μ,α)*delta(λ,β)*delta(ν,ρ)-delta(λ,α)*delta(ν,β)*delta(μ,ρ)-delta(ν,α)*delta(μ,β)*delta(λ,ρ)
end

function mainterm(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64,α::Int64,β::Int64,μ::Int64,ν::Int64,N::Int64)
    return 24*(3*x_nc(N-1)-2)*delta(β,ν)*delta(μ,a)*delta(e,α)*delta3(f,g,h,b,c,d)-delta(μ,a)*delta(e,α)*levicivita([f,g,h,ν])*levicivita([b,c,d,β])
end



function RDMC(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64,α::Int64,β::Int64,μ::Int64,ν::Int64,N::Int64)
    return mainterm(a,b,c,d,e,f,g,h,α,β,μ,ν,N)+mainterm(b,a,c,d,f,e,g,h,α,β,μ,ν,N)-mainterm(b,a,c,d,e,f,g,h,α,β,μ,ν,N)-mainterm(a,b,c,d,f,e,g,h,α,β,μ,ν,N)
end

function RHOC(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},L::Int64)
    #=structured with the first two indices α and β and the next two μ and ν; final two indices are RDM indices with hilbert space of the boundary site between block and bulk=#
    rho = zeros(Float64,4,4,4,4,36,36)
    for α in 1:4, β in 1:4, μ in 1:4, ν in 1:4, i in 1:36, j in 1:36
        a = BasisC[i][1][1]
        b = BasisC[i][1][2]
        c = BasisC[i][2][1]
        d = BasisC[i][2][2]
        e = BasisC[j][1][1]
        f = BasisC[j][1][2]
        g = BasisC[j][2][1]
        h = BasisC[j][2][2]
        rho[α,β,μ,ν,i,j] = LeftEdge[μ]*RightEdge[ν]*LeftEdge[α]*RightEdge[β]*RDMC(a,b,c,d,e,f,g,h,α,β,μ,ν,L)
    end
    test = dropdims(sum(rho, dims = (1,2,3,4)),dims = (1,2,3,4))
    fnl = test*(1/tr(test))
    return fnl
end

#This whole function is defined to avoid an error message
function nosing(N::Int64,L::Int64)
    if iseven(N) && 2<N<L-1
        return Integer((L-2)/2)
    elseif isodd(N) && 2<N<L-1
        return Integer((L-N)/2)
    end
end
nosing(5,10)
function switch(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},N::Int64#=position of bipartition (ie next site)=#,L::Int64#=sites in chain=#)
    if !isodd(L)
        print("error, please enter an odd chain length")
        return 0
    end

    S = Integer((L-1)/2)
    if N == 2
        return RHOB(LeftEdge,RightEdge,S,0) #I know this looks messed up but I coded it backwards, sue me.
    elseif L == N
        return RHOB(RightEdge,LeftEdge,S,0)
    elseif 3<N<L-1 && iseven(N)
        return RHOA(LeftEdge,RightEdge,S,Integer((N-2)/2))
    elseif 3<N<L-1 && isodd(N)
        RHOA(RightEdge,LeftEdge,S,Integer((L-N)/2))
    elseif N == 3
        RHOC(LeftEdge,RightEdge,S)
    elseif N == L-1
        RHOC(RightEdge,LeftEdge,S)
    end
end
false & false
SvN((1. +0*im)switch(LET,RET,11,11),0.01)

SvN((1. +0*im)RHOA(LET,RET,20,1),0.01)
SvN((1. +0*im)RHOB(LET,RET,10,0),0.01)
SvN((1. +0*im)RHOC(LET,RET,10),0.01)

ChainLength =3

XAX = []
YAX = []
for i in 1:10000
    LET = normalize(rand(Float64,4))
    RET = normalize(rand(Float64,4))
    AGR = (dot(LET,RET))^2
    push!(XAX,AGR)
    DMAT = switch(LET,RET,2,ChainLength)
    push!(YAX,SvN((1. +0*im)DMAT,0.01))
end

filename = string("Analytic_EE_",ChainLength,"Site.txt")
touch(filename)
savearray = [XAX,YAX]
display(savearray)

open(filename,"w") do io
    writedlm(filename,savearray)
end



filename2 = "/Users/alexjacoby/Documents/xVBS_Code/xVBSFigs/EDGEFIG/Analytic_EE_3Site.txt"
dat = open(filename2) do io
    readdlm(filename2)
end
XAX = dat[1,:]
YAX= dat[2,:]


plot(XAX,YAX, xlabel = L"\mathcal{E}_{\mu\nu}\mathcal{E}^{\nu\mu}",ylabel = L"\mathbf{S}_{vN}" ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15,markercolor = :black,markersize =4, seriestype = :scatter, fontfamily = "Times", label = "ED",legend = :none,markeralpha = 0.5, markerstrokealpha = 1)

savefig("EDGEFIG.pdf")
