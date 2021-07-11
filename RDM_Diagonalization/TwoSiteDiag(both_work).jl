using LinearAlgebra
using Combinatorics
using Plots
#=REINDEX THE POSITIONS YOU FORGOT T SAVE A LOT OF THIS=#


#Ls = 10
#L = 2*Ls+1
BasisC =([[1,2],[1,2]],[[1,3],[1,2]],[[1,4],[1,2]],[[2,3],[1,2]],[[2,4],[1,2]],[[3,4],[1,2]],[[1,2],[1,3]],[[1,3],[1,3]],[[1,4],[1,3]],[[2,3],[1,3]],[[2,4],[1,3]],[[3,4],[1,3]],[[1,2],[1,4]],[[1,3],[1,4]],[[1,4],[1,4]],[[2,3],[1,4]],[[2,4],[1,4]],[[3,4],[1,4]],[[1,2],[2,3]],[[1,3],[2,3]],[[1,4],[2,3]],[[2,3],[2,3]],[[2,4],[2,3]],[[3,4],[2,3]],[[1,2],[2,4]],[[1,3],[2,4]],[[1,4],[2,4]],[[2,3],[2,4]],[[2,4],[2,4]],[[3,4],[2,4]],[[1,2],[3,4]],[[1,3],[3,4]],[[1,4],[3,4]],[[2,3],[3,4]],[[2,4],[3,4]],[[3,4],[3,4]])
function delta(a::Int64,b::Int64)
    if a == b
        return 1
    else
        return 0
    end
end

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

LET = [1.,0,0,0]
RET = [1.,0,0,0]

test =RHOC(LET,RET,10)
eigvals(test)
tr(test)
SvN((1. + 0*im)*test,0.01)
