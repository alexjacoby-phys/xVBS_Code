using LinearAlgebra
using Combinatorics
using Plots
#=REINDEX THE POSITIONS YOU FORGOT T SAVE A LOT OF THIS=#


#Ls = 10
#L = 2*Ls+1
Basis = [1,2,3,4]
function delta(a::Int64,b::Int64)
    if a == b
        return 1
    else
        return 0
    end
end

function x_n(L::Int64)
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


#=term 1 is just δ^{abc]d}_{e[fg]h}=#
function term1(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return delta(a,c)*delta(β,ν)*delta(b,d)
end


#=term 2 is δ^{ab[cd]}_{e[fg]h}=#
function term2(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return delta(a,c)*(2*delta(β,d)*delta(ν,b)+delta(β,ν)*delta(b,d))
end
#=term 3 is δ^{[ab]cd}_{e[fg]h}=#
function term3(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return (2*delta(ν,a)*delta(β,c)+delta(a,c)*delta(β,ν))*delta(d,b)
end
#=term 4 is δ^{[ab]cd}_{ef[gh]}=#
function term4(a::Int64,b::Int64,c::Int64,d::Int64,β::Int64,ν::Int64)
    return (delta(a,c)*delta(b,ν)*delta(d,β)+delta(a,ν)*delta(β,c)*delta(b,d)-delta(a,ν)*delta(b,c)*delta(β,d)-delta(β,c)*delta(a,d)*delta(b,ν))
end


function RDM(a::Int64#=μ_{i}=#,b::Int64#=ν_{f}=#,c::Int64#=α_{i}=#,d::Int64#=β_{f}=#,β::Int64,ν::Int64,N::Int64#=Number of singlets in chain=#,L::Int64#=Number of singlets into chain for bipartition=#)
    c1 = x_n(L)
    c2 = x_n(N-L)
    RDMel = 6(c1-1)*(c2-1)*term1(a,b,c,d,β,ν)+(c1-1)*term2(a,b,c,d,β,ν)+(c2-1)*term3(a,b,c,d,β,ν)+term4(a,b,c,d,β,ν)
    return RDMel
end
#=make sure to be careful not clear what this is normalized with respect to =#
LeftEdge = normalize([1.,1.,0,0])#normalize(2*rand(4)-ones(Float64,4))
RightEdge = normalize([1.,1.0,0,0])#normalize(2*rand(4)-ones(Float64,4))
AGR = (dot(LeftEdge,RightEdge))^2


function RHO(LeftEdge::Vector{Float64},RightEdge::Vector{Float64},N::Int64,L::Int64)
    Basis = [1,2,3,4]
    #=structured with the first two indices α and β and the next two μ and ν; final two indices are RDM indices with hilbert space of the boundary site between block and bulk=#
    rho = zeros(Float64,4,4,4,4,4,4)
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4, β in 1:4, ν in 1:4
        rho[a,b,c,d,β,ν] = LeftEdge[a]*RightEdge[b]*LeftEdge[c]*RightEdge[d]*RDM(a,b,c,d,β,ν,N,L)
    end
    test = dropdims(sum(rho, dims = (1,2,3,4)),dims = (1,2,3,4))
    fnl = test*(1/tr(test))
    return fnl
end





RHO(LeftEdge,RightEdge,6,2)
L=6
YAX = zeros(L-1)
for i in 1:L-1
    YAX[i] = SvN((1. +0*im)RHO(LeftEdge,RightEdge,L,i),0.01)
end
XAX = 2*Vector(1:L-1)+ones(Float64,L-1)
prepend!(XAX,1.)
append!(XAX,13.)
prepend!(YAX,1.1)
append!(YAX,1.1)
YAX
plot(XAX,YAX,title = string("SvN of xVBS Chain of Length ",2L+1), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")


savefig("13sitepenandpaper")

N = 100
XAX = []
YAX = []
for i in 1:N
    LeftEdge = normalize(rand(Float64,4)-0.5*ones(Float64,4))
    RightEdge = normalize(rand(Float64,4)-0.5*ones(Float64,4))
    append!(XAX,(dot(LeftEdge,RightEdge))^2)
    append!(YAX,SvN((1. +0*im)*RHO(LeftEdge,RightEdge,6,3),0.01))
end

plot(XAX,YAX, seriestype = :scatter)
