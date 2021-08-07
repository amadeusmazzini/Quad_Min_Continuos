import Pkg; Pkg.add("QuadGK")
using LinearAlgebra, QuadGK, Plots
gr()

f(x) = (1+sqrt(2+x))/(1+sqrt(1+x^2)) +  sin(5*π*x)*cos(5x)x^2/(x^2+1)
w(x) = 1
a , b = 0 , 1
grau = 5


Φ = Any[]
B = zeros(grau+1)
C = zeros(grau+1)
push!(Φ, x -> 1)
g(x) = x * w(x) * Φ[1](x)^2
g₁(x) = w(x) * Φ[1](x)^2
B[1] = quadgk(g, a, b, rtol=1e-8)[1] / quadgk(g₁, a, b, rtol=1e-8)[1]
push!(Φ, x-> x - B[1])

for k = 3:grau+1
    h(x)  = x * w(x) * Φ[k-1](x)^2
    h₁(x) = w(x) * Φ[k-1](x)^2
    c(x)  = x * w(x) * Φ[k-1](x) * Φ[k-2](x)
    c₁(x) = w(x) * Φ[k-2](x)^2
    B[k]  = quadgk(h, a, b, rtol=1e-8)[1] / quadgk(h₁, a, b, rtol=1e-8)[1]
    C[k]  = quadgk(c, a, b, rtol=1e-8)[1] / quadgk(c₁, a, b, rtol=1e-8)[1]
    push!(Φ, x-> (x - B[k])Φ[k-1](x)-C[k]Φ[k-2](x))
end

u(x) = f(x) * w(x)
α = zeros(grau+1)
for n = 1:grau+1
    aux1(x) = u(x) * Φ[n](x)
    aux2(x) = w(x) * Φ[n](x)^2
    α[n] = quadgk(aux1, a, b, rtol=1e-8)[1]/quadgk(aux2, a, b, rtol=1e-8)[1]
end
o(x)=sum(α[k]*Φ[k](x) for k=1:grau+1)

m = range(a,b,length=100)
plot(m,f.(m),lab=:"f(x)",)
plot!(m,o.(m),legend=:topleft,lab=:"Aproximação")