using LinearAlgebra

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
asm[2]
for i in 1:6
    BCOps[i] = BCOps[i]-asm[i]
end

BCOps
