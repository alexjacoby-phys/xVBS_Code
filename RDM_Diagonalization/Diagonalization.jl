using LinearAlgebra
using Combinatorics

#=REINDEX THE POSITIONS YOU FORGOT T SAVE A LOT OF THIS=#


#Ls = 10
#L = 2*Ls+1
Basis = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
function delta(a::Int64,b::Int64)
    if a == b
        return 1
    else
        return 0
    end
end

function x_n(L::Int64)
    coefficients = zeros(Float64,L)
    coefficients[1] = 3
    for i in 2:L
        coefficients[i] = 3*(3*coefficients[i-1]-2)
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


#=term 1 is just δ^{a[bc]d}_{e[fg]h}=#
function term1(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return 2*delta(a,e)*delta(d,h)*(delta(b,f)*delta(c,g)-delta(b,g)*delta(c,f))
end


#=term 2 is δ^{ab[cd]}_{e[fg]h}=#
function term2(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return delta(a,e)*(delta(b,f)*delta(c,g)*delta(d,h)+delta(b,g)*delta(d,f)*delta(c,h)-delta(b,g)*delta(c,f)*delta(d,h)-delta(b,f)*delta(d,g)*delta(c,h))
end
#=term 3 is δ^{[ab]cd}_{e[fg]h}=#
function term3(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return (delta(a,e)*delta(b,f)*delta(c,g)+delta(b,e)*delta(a,g)*delta(c,f)-delta(b,e)*delta(a,f)*delta(c,g)-delta(a,e)*delta(b,g)*delta(c,f))*delta(d,h)
end
#=term 4 is δ^{[ab]cd}_{ef[gh]}=#
function term4(a::Int64,b::Int64,c::Int64,d::Int64,e::Int64,f::Int64,g::Int64,h::Int64)
    return (delta(a,e)*delta(b,f)-delta(b,e)*delta(a,f))*(delta(c,g)*delta(d,h)-delta(c,h)*delta(d,g))
end


function RDM(a::Int64#=μ_{i}=#,b::Int64#=ν_{f}=#,c::Int64#=α_{i}=#,d::Int64#=β_{f}=#,α::Int64,β::Int64,μ::Int64,ν::Int64,N::Int64#=Number of singlets in chain=#,L::Int64#=Number of singlets into chain for bipartition=#)
    tempcof = x_n(N)
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
    RDMel = (c1-1)*(c2-1)*term1(a,α,β,b,c,μ,ν,d)+(c1-1)*term2(a,b,α,β,c,ν,μ,d)+(c2-1)*term3(α,β,a,b,c,ν,μ,d)+term4(α,β,a,b,c,d,μ,ν)
    return RDMel
end
#=make sure to be careful not clear what this is normalized with respect to =#
LeftEdge = normalize([1.,0,0,0])#normalize(2*rand(4)-ones(Float64,4))
RightEdge = normalize([0,1.,0,0])#normalize(2*rand(4)-ones(Float64,4))
AGR = (dot(LeftEdge,RightEdge))^2


function RHO(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},N::Int64,L::Int64)
    Basis = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    #=structured with the first two indices α and β and the next two μ and ν; final two indices are RDM indices with hilbert space of the boundary site between block and bulk=#
    rho = zeros(Float64,4,4,4,4,6,6)
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4, e in 1:6, f in 1:6
        μ = Basis[e][1]
        ν = Basis[e][2]
        α = Basis[f][1]
        β = Basis[f][2]
        rho[a,b,c,d,e,f] = LeftEdge[a]*RightEdge[b]*LeftEdge[c]*RightEdge[d]*RDM(a,b,c,d,α,β,μ,ν,N,L)
    end
    test = dropdims(sum(rho, dims = (1,2,3,4)),dims = (1,2,3,4))
    fnl = test*(1/tr(test))
    return fnl
end



print(SvN(RHO(LeftEdge,RightEdge,1,0)))

eigvals(RHO(LeftEdge,RightEdge,1,0))


L=10
YAX = zeros(L+1)
for i in 1:L+1
    YAX[i] = SvN(RHO(LeftEdge,RightEdge,L,i-1))
end
XAX = 2*Vector(0:L)+ones(Float64,L+1)
YAX
plot(XAX,YAX,title = string("SvN of xVBS Chain of Length ",2*L+1), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
