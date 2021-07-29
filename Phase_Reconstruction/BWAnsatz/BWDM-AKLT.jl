using LinearAlgebra
using Combinatorics
using SparseArrays
using Arpack
#using Plots
#using DelimitedFiles


#creates list of operators to pass to makeop#

X = [0. 1 0 ; 1 0 1 ; 0 1 0]*(1/sqrt(2))*(1+0*im)
Y = [0. -im 0 ; im 0 -im ; 0 im 0]*(1/sqrt(2))*(1+0*im)
Z = [1. 0 0 ; 0 0 0; 0 0 -1]*(1+0*im)
plus = [1. 0 0 ; 0 0 0; 0 0 0]*(1+0*im)
minus = [0. 0 0 ; 0 1 0; 0 0 0]*(1+0*im)
projplus = [1. 0 0 ; 0 0 0 ; 0 0 0]*(1+0*im)
projzero = [0. 0 0 ; 0 1 0 ; 0 0 0]*(1+0*im)
projminus = [0. 0 0 ; 0 0 0 ; 0 0 1]*(1+0*im)
Id = [1. 0 0 ; 0 1 0 ; 0 0 1]*(1+0*im)



L = 7
bcopsarray = []
for i in 1:L-1
    push!(bcopsarray,sparse(Id))
end
push!(bcopsarray,L*sparse(projplus+projzero))
bcopsarray
BC = kron(bcopsarray...)


a = 1.
b = Float64(1/3)
BL = kron(X,X)+kron(Y,Y)+kron(Z,Z)
CM = a*(BL)+b*(BL*BL)+Float64(2/3)*kron(Id,Id)
CMs = sparse(CM)
Ids = sparse(Id)
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
@time H = sum(Harray, dims = 1)[1]-100*BC
@time states = eigs(H,which = :SR, nev = 10 )
ev = 2π*(real(states[1])-real(states[1])[1]*ones(6))
diagm(ev)
notrho = 2.71^(-diagm(ev))













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


#plot(XAX,YAX,title = string("SvN of xVBS Chain of Length ",L, ", N =",n), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")
#savefig("5siteED")
