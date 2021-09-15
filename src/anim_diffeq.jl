using SippyArt
using DifferentialEquations
using ModelingToolkit
using FixedPointNumbers
using Colors, Images
using Distributions, QuasiMonteCarlo
using BenchmarkTools
using Test
using Images, Colors
using JLD2

vid_name = "lorenz_81"
# p = [10.0,28.0,8 / 3]
u0 = [1.0;0.0;0.0]
tspan = (0.0, 10.0)
# prob = ODEProblem(lorenz, u0, tspan, p)

# p = [0.32899]
# u0 = rand(3)
# p = collect((a = 0.95, b = 0.7, c = 0.6, d = 3.5, e = 0.25, f = 0.1))
# p = collect((a=5, b=-10, c=-0.38))
# p = collect((a=10, b=-28, c=8/3))


# tspan = (0.0, 1e4)

prob = ODEProblem(SippyArt.lorenz, u0, tspan, p)


# lb = zeros(3)
# ub = ones(3)
# N = 100
# n = N^2
# u0s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())

img_fn = "/Users/anand/Pictures/81.jpeg"
img_tmp = load(img_fn)
img_tmp = convert(Matrix{RGB{N0f8}}, img_tmp)
img = img_tmp[1:50:end-50, 1:50:end-50]
img_size = size(img)
pix_to_u0(pix) = convert(Array{Float64}, [pix.r, pix.g, pix.b])

# function prob_func(prob, i, repeat)
#     remake(prob, u0=@view(u0s[:, i]))
# end
function prob_func2(prob, i, repeat)
    remake(prob, u0=pix_to_u0(img[i]))
end

# ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func2)

# sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=n, saveat=tspan[1]:0.001:tspan[2])
sim = solve(ensemble_prob, EnsembleThreads(), trajectories=length(img), saveat=tspan[1]:0.001:tspan[2])

@save "data/sols/$(vid_name).jld2" sol=sim
# colorview(RGB, Array(sim))

pixs = map(y -> map(x -> RGB(x...), y), sim) # todo try colorview
frame = reshape(first.(pixs), size(img))
# frame = reshape(first.(pixs), (100, 100))
@test frame == img

frames = Matrix{RGB}[]
for i in 1:length(sim.u[1])
    # frame = reshape(getindex.(pixs, i), (N, N))
    frame = reshape(getindex.(pixs, i), img_size)
    frame = map(clamp01nan, frame)
    frame = convert.(RGB{N0f8}, frame)
    push!(frames, frame)
    save("data/out4/$i.jpeg", frame)
end

cmd = Cmd(["ffmpeg -framerate 60 -start_number 1 -i 'data/out4/%d.jpeg' -r 60", "-y", "data/$(vid_name).mp4"])
run(cmd)
