import Pkg; Pkg.add("QuadGK")
using QuadGK, Plots
gr()

# ALGORITMO PARA APROXIMAR FUNÇÕES UTILIZANDO POLINÔMIOS DE CHEBYSHEV

#f(x) = (1+sqrt(2+x))/(1+sqrt(1+x^2))+sin(5*π*x)*cos(5x)x^2/(x^2+1) #(SUGESTÃO)
f(x) = x*sin(x)        ## Digite a função f(x) que deseja aproximar   
a, b = 0, 2                             ## Digite o intervalo [a,b] que deseja
grau = 5                                ## Digite o Grau da aproximação



chebyshev(x,n) = cos(n*acos(x)) # Polinômio de Chebyshev 
z(x) = (a + b) / 2 + x * (b - a) / 2 #Transformação de variáveis
xtf(z) = (2z - a - b) / (b - a) #Inversa de z
g(x) = f(z(x))
w(x) = 1/(sqrt(1-x^2)) #Função Peso
h(x) = g(x) * w(x)
θ = zeros(grau+1) # θ será o vetor de coeficientes
θ[1] = quadgk(h, -1+1e-12, 1-1e-12, rtol=1e-8)[1]/π
for n = 1:grau
    Tn(x) = h(x) * chebyshev(x,n)
    θ[n+1]= 2 * quadgk(Tn, -1+1e-12, 1-1e-12, rtol=1e-8)[1]/π
end
o(x) = sum(θ[k]*chebyshev(x,k-1) for k=1:grau+1)
otil(z) = o(xtf(z))
zg = range(a, b, length=100)
Error(x) = w(x)*(f(x)-otil(x))^2 ## Função Erro
E = quadgk(Error, -1+1e-12, 1-1e-12, rtol=1e-8)[1] #Erro da Aproximação
plot(zg,f.(zg),lab=:"f(x)") #Plot de f
plot!(zg,otil.(zg),legend=:bottomright,lab=:"Aproximação Por Chebyshev") #Plot da Aproximação
