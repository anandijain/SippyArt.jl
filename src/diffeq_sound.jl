# used to create samples, used in https://8102.bandcamp.com/track/8-21-21

using WAV, DifferentialEquations, SippyArt
using ModelingToolkit
sr=44100
u0 = [1.0;1.0]
p = [1.5,1.0,3.0,1.0]
tspan = (0.0, 10000.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)
saveat = tspan[1]:0.1:tspan[2]
sol = solve(prob; saveat=saveat)

fn = "/Users/anand/Music/samples/lotka.wav"
sol_to_wav(sol, fn; sr=sr)

w, sr = wavread("/Users/anand/ableton_files/081821081821.wav")
w, sr = wavread("/Users/anand/ableton_files/081821 Project/081821.wav")

function lorenz!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end

p = [10.0,28.0,8 / 3]
u0 = [1.0;0.0;0.0]
tspan = (0.0, 1000.0)
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob; saveat=tspan[1]:0.001:tspan[2])
arr = Array(sol)
for i in 1:length(u0)
    arr[i, :] .= arr[i, :] ./ maximum(arr[i, :])
end
arr = 2*arr .- 1
wavwrite(arr[1:2, :]', 44100, "foo.wav")

using Catalyst, DifferentialEquations, WAV

repressilator = @reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ;

odesys = convert(ODESystem, repressilator)

# parameters [α,K,n,δ,γ,β,μ]
p = (.5, 40, 2, log(2)/120, 5e-3, 20*log(2)/120, log(2)/60)

# initial condition [m₁,m₂,m₃,P₁,P₂,P₃]
u0 = [0.,0.,0.,20.,0.,0.]

# time interval to solve on
tspan = (0., 100000.)

# create the ODEProblem we want to solve
prob = ODEProblem(repressilator, u0, tspan, p)
sol = solve(prob; saveat=tspan[1]:0.1:tspan[2])
arr = Array(sol)
for i in 1:length(u0)
    arr[i, :] .= arr[i, :] ./ maximum(arr[i, :])
end
arr = 2*arr .- 1
for i in 1:2:length(u0)-1
    wavwrite(arr[i:i+1, :]', 44100, "repressilator3_$i.wav")
end

# redefine the initial condition to be integer valued
u₀ = [0,0,0,20,0,0]

# next we create a discrete problem to encode that our species are integer valued:
dprob = DiscreteProblem(repressilator, u₀, tspan, p)

# now, we create a JumpProblem, and specify Gillespie's Direct Method as the solver:
jprob = JumpProblem(repressilator, dprob, Direct(), save_positions=(false,false))

sol = solve(jprob, SSAStepper(); saveat=tspan[1]:0.1:tspan[2])
arr = Float64.(Array(sol))
for i in 1:length(u0)
    arr[i, :] .= arr[i, :] ./ maximum(arr[i, :])
end
arr = 2*arr .- 1
for i in 1:2:length(u0)-1
    wavwrite(arr[i:i+1, :]', 44100, "data/repressilator_jump_3_$i.wav")
end

using ModelingToolkit, LinearAlgebra, DifferentialEquations, Plots, Test 

function Pendulum(;name, g=9.81, r=1.)
    gval, rval = g, r
    @parameters g r
    @variables x(t) y(t) T(t)
    D2 = Differential(t)^2
    eqs = [ D2(x) ~ T * x
            D2(y) ~ T * y - g ]
    ODESystem(eqs, t, [x, y, T], [g, r]; name=name, defaults=(g => gval, r => rval))
end

function connect_pendulums(pends)
    eqs = [0 ~ pends[1].x^2 + pends[1].y^2 - pends[1].r^2]
    for i in 2:length(pends)
        push!(eqs, connect_eq(pends[i - 1:i]...))
    end
    eqs
end

connect_eq(a, b) = 0 ~ (a.x - b.x)^2 + (a.y - b.y)^2 - b.r^2 

@parameters t
D = Differential(t)

N = 2 # change this to add or remove pendulumns
pends = []
for i in 1:N
    push!(pends, Pendulum(;name=Symbol("pendulum", i)))
end

npendeqs = connect_pendulums(pends)

@named npend = compose(ODESystem(npendeqs, t; name=:connections), pends...)
@named npend2 = ODESystem(npendeqs, t; systems=pends)

u0 = Pair[]
for (i, p) in enumerate(pends)
    pairs = [
        p.x => i
        p.y => 0
        D(p.x) => 0
        D(p.y) => 0
        p.T => 0
    ]
    append!(u0, pairs)
end

pendulum_sys = structural_simplify(ode_order_lowering(flatten(npend)))
prob = ODAEProblem(pendulum_sys, u0, (0, 1000.0))

sol = solve(prob; saveat=tspan[1]:0.1:tspan[2])



plot(sol)
sol = solve(prob)

#todo duffing and van der pol



p = [0.1, 0.1, 1.4]
u0 = [1.0; -1.0]
tspan = (0.0, 1e5)


prob = ODEProblem(duffing!, u0, tspan, p)
sys = modelingtoolkitize(prob)
sol = solve(prob; saveat=tspan[1]:0.01:tspan[2])
sol_to_wav(sol, "data/wav/duffing.wav")



p = [0.32899]
u0 = rand(3)
tspan = (0.0, 1e4)


prob = ODEProblem(cyclically_symmetric, u0, tspan, p)
sys = modelingtoolkitize(prob)
sol = solve(prob; saveat=tspan[1]:0.01:tspan[2])
sol_to_wav(sol, "data/wav/cyclically_symmetric.wav")

p = [0.2, 0.2, 5.7]
u0 = rand(3)
tspan = (0.0, 1e4)


prob = ODEProblem(rossler, u0, tspan, p)
sys = modelingtoolkitize(prob)
sol = solve(prob; saveat=tspan[1]:0.01:tspan[2])
sol_to_wav(sol, "data/wav/rossler.wav")


p = [1, 3, 5, 1, 5, 1, 4, -8/5]
u0 = rand(3)
tspan = (0.0, 1e4)


prob = ODEProblem(hindmarsh_rose, u0, tspan, p)
sys = modelingtoolkitize(prob)
sol = solve(prob; saveat=tspan[1]:0.01:tspan[2])
sol_to_wav(sol, "data/wav/hindmarsh_rose.wav")

p = (a = 0.95, b = 0.7, c = 0.6, d = 3.5, e = 0.25, f = 0.1)
u0 = rand(3)
tspan = (0.0, 1e4)


prob = ODEProblem(aizawa, u0, tspan, p)
sys = modelingtoolkitize(prob)
sol = solve(prob; saveat=tspan[1]:0.01:tspan[2])
sol_to_wav(sol, "data/wav/aizawa.wav")


# discrete 
using WAV, DifferentialEquations, SippyArt, Plots
using ModelingToolkit
sr = 88200
@parameters t a=440 b=439
# @parameters c nsteps δt
# Diff = Difference(t; dt=1//sr)
Diff = Difference(t; dt=1) 
D2 = Differential(t)^2
D = Differential(t)
# Diff2 = Difference(t; dt=1)^2
@variables x(t)=0.5 y(t)=0.5
tspan = (0.0,10*sr) 

eqs = [
    Diff(x) ~ 1-a*x^2+y,
    Diff(y) ~ b*x
]

eqs = [
    Diff(x) ~ cos(a*x),
    Diff(y) ~ sin(b*y)
]

@named sys = DiscreteSystem(eqs, t, [x, y], [a,b]) 
u0 = [x =>0.5, y=>0.5]
p = [a=>440,b=>439]
prob = DiscreteProblem(sys, u0, tspan, p;saveat=tspan[1]:tspan[2])
sol = solve(prob,FunctionMap())
arr= Array(sol)
arr2 = sol_to_wav(sol, "sincos.wav";sr=sr)

arr = sin.(1:10*sr)
wavwrite(arr, "idk.wav")



eq = D2(x) ~ -a*x
eq = D(x) ~ sin(x)
@named sys = ODESystem(eq)
sys = ode_order_lowering(sys)

prob = ODEProblem(sys, [D(x)=>0], (0, 10.); saveat=0.:1//sr:10)
prob = ODEProblem(sys, [], (0, 10.); saveat=0.:1//sr:10)
sol = solve(prob)
plot(sol; vars=(1,2))
plot(sol)
sol_to_wav(sol[x], "foodiffeq.wav";sr=sr)


using DiffEqProblemLibrary, SippyArt, ModelingToolkit, DifferentialEquations
DiffEqProblemLibrary.ODEProblemLibrary.importodeproblems()
@parameters t b=0.208186
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ sin(y) - b*x,
       D(y) ~ sin(z) - b*y,
       D(z) ~ sin(x) - b*z]

@named thomas = ODESystem(eqs)
tspan = (0.0,100.0)

@parameters t a=0.95 b=0.7 c=0.6 d=3.5 e=0.25 f=0.1
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ (z-b)*x-d*y,
       D(y) ~ d*x+(z-b)*y,
       D(z) ~ c+a*z-z^3/3-(x^2+y^2)*(1+e*z)+ f*z*x^3]

@named aizawa = ODESystem(eqs)

@parameters t a=3 b=2.7 c=1.7 d=2 e=9
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ y-a*x+b*y*z,
       D(y) ~ c*y-x*z+z,
       D(z) ~ d*x*y-e*z]

@named dadras = ODESystem(eqs)

@parameters t a=35 b=3 c=28
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ a*(y-x),
       D(y) ~ (c-a)*x-x*z+c*y,
       D(z) ~ x*y-b*z]

@named chen = ODESystem(eqs)


@parameters t a=0.14 b=0.10
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ y*(z-1+x^2) + b*x,
       D(y) ~ x*(3*z+1-x^2)+b*y,
       D(z) ~ -2*z*(a+x*y)]

@named rabinovich_fabrikant = ODESystem(eqs)

@parameters t a=2.07 b=1.79
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ y+a*x*y+x*z,
       D(y) ~ 1-b*x^2+y*z,
       D(z) ~ x-x^2-y^2]

@named sprott = ODESystem(eqs)

@parameters t a=1 b=3 c=1 d=5 r=1e-2 s=4 xr=-8/5 i=5
@variables x(t)=1 y(t)=0 z(t)=0
D = Differential(t)

eqs = [D(x) ~ y-a*x^3+b*x^2-z+i,
       D(y) ~ c-d*x^2-y,
       D(z) ~ r*(s*(x-xr)-z)]

@named hindmarsh_rose = ODESystem(eqs)

@parameters t p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12
@variables y1(t) y2(t) y3(t) y4(t) y5(t) y6(t) y7(t) y8(t)
D = Differential(t)
eqs = [D(y1) ~ -p1*y1 + p2*y2 + p3*y3 + p4,
       D(y2) ~ p1*y1 - p5*y2,
       D(y3) ~ -p6*y3 + p2*y4 + p7*y5,
       D(y4) ~ p3*y2 + p1*y3 - p8*y4,
       D(y5) ~ -p9*y5 + p2*y6 + p2*y7,
       D(y6) ~ -p10*y6*y8 + p11*y4 + p1*y5 -
                p2*y6 + p11*y7,
       D(y7) ~  p10*y6*y8 - p12*y7,
       D(y8) ~ -p10*y6*y8 + p12*y7]
de = ODESystem(eqs; name=:hires)

function sys_to_wav(sys, tspan; sr=88200)
    prob = ODEProblem(sys,[],tspan; saveat=tspan[1]:100//sr:tspan[2])
    sol = solve(prob)
    sol_to_wav(sol, "/Users/anand/Music/samples/$(nameof(sys)).wav";sr=sr)
    # plot(sol)
    # sys, prob, sol
    nothing
end

sys_to_wav(hindmarsh_rose, tspan)

syss= [thomas, aizawa, dadras, chen, rabinovich_fabrikant]
for sys in syss
    sys_to_wav(sys, tspan)
end
