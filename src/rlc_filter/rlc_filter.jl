
using DifferentialEquations, ModelingToolkit, Symbolics, ModelingToolkitStandardLibrary, LinearAlgebra
using ModelingToolkitStandardLibrary: Blocks, Blocks.Constant
import ModelingToolkitStandardLibrary: Blocks, Blocks.Constant
using ModelingToolkitStandardLibrary.Electrical

using WAV

using DataInterpolations
using Plots

function ConstantInterp(; name, interp)
    @named output = Blocks.RealOutput()
    eqs = [
        output.u ~ interp(t),
    ]
    compose(ODESystem(eqs; name=name), [output])
end


const MTK = ModelingToolkit
const MSL = ModelingToolkitStandardLibrary

@parameters t
D = Differential(t)
# @variables Lt(t) Ct(t) = 100 Vt(t) = 100
# s = wavread("square.wav")
ws, sr = wavread("/Users/anand/ableton_files/081821 Project/081821.wav")
w = ws[:, 1]
ws, sr = wavread("10secs.wav")
w = ws[:, 1]
tspan = (0, length(w) / sr-1)
ts = 0:1/sr:tspan[2]

interp = ConstantInterpolation(100w, ts)
# 2 cases, control voltage based on wav file or equation
# eq = sig ~ sin(sin())

@named r = Resistor(R=1)
@named v = Voltage()
@named c = Capacitor(C=1)
@named g = Ground()
@named l = Inductor(L=1)


interp_(t) = interp(t)
@register_symbolic interp_(t)
@named sig = ConstantInterp(interp=interp)
rlc_eqs = [
    # v.V ~ interp_(t)
    connect(sig.output, v.V)
    connect(v.p, r.p)
    connect(r.n, l.p)
    connect(l.n, c.p)
    connect(c.n, v.n, g.g)]

@named rlc = ODESystem(rlc_eqs, t; systems=[sig, v, r, l, c, g])
sys = structural_simplify(rlc)
Ls = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
for L in Ls
    prob = ODEProblem(sys, [l.L => L], tspan)
    sol = solve(prob;saveat=ts)
    gen_dir = "/Users/anand/Music/gen/"
    signal = Array(sol[1, :])
    sig2 = 2 * (((signal .- minimum(signal)) / (maximum(signal) - minimum(signal)) .- 0.5))
    wavwrite(sig2, joinpath(gen_dir, "RLC_$L.wav"); Fs=sr)

end