import Pkg; Pkg.add("QuadGK")
using LinearAlgebra, QuadGK
using Plots
gr()

grau = 2            #Digite o grau do polinômio
a,b = 0 , 1         #Digite o intervalo [a,b] 
f(x) = 4x^3         #Digite a f(x) que deseja



M = zeros(grau+1,grau+1)
B = zeros(grau+1)
for i = 1:grau+1  
    for j =i:grau+1
        h(x) = x^(i-1)*x^(j-1)
        M[i,j] = quadgk(h, a, b)[1]
        M[j,i] = M[i,j]
    end
    g(x) = f(x)*x^(i-1)
    B[i] = quadgk(g, a, b)[1]
end
α = M\B             #Resolve M α = B
P(x) = sum(α[k+1]*x^(k) for k = 0:grau) #Polinômio aproximador
m = range(0,1, length=100)
plot(m,f.(m),lab="f(x)")
plot!(m,P.(m),lab="Aproximação",legend=:bottomright)
