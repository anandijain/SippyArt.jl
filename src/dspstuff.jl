# some of thi scode was used in https://www.youtube.com/watch?v=2T7vzIOpS68
using SippyArt 
using WAV, DSP, Plots
using Images, FileIO
using BenchmarkTools
using ColorVectorSpace, ColorTypes, ColorBlendModes
using Dates
using ImageFiltering, ImageContrastAdjustment
using VideoIO

# using ImageSegmentation
# using DitherPunk

proj_dir = joinpath(datadir, string(today()))
proj_dir = "data/2021-09-05"
src_dir = joinpath(proj_dir, "sources")
out_dir = joinpath(proj_dir, "outputs")

fps = 60
# w, sr = wavread(joinpath(src_dir, "081821.wav"))
img = load("/Users/anand/Pictures/81.jpeg")
img2 = copy(img)
w, sr = wavread(joinpath(src_dir, "082121.wav"))
sr = Int(sr)
N = size(w, 1)
seconds = N / sr
nframes = round(Int, fps *  seconds, RoundUp)
s = w[:, 1] # just do one channel for now 
samples_per_frame = N / nframes
spf = round(Int, samples_per_frame)
resolution = (1000, 1000)
img_W, img_H = resolution


plot(s[1:sr])
plot(s[1:fps*spf])

ws = arraysplit(s, spf, 1)
wss = collect(ws)
avgs = mean.(wss)
plot(avgs)


center = rand(RGB{N0f8}, 500, 500) - rand(RGB{N0f8}, 500, 500)

frames_to_render = 10*fps # nseconds
center = rand(Gray, 500, 500) - rand(Gray, 500, 500)

function morpho_vid(s, frames_to_render)
    itval = round(Int, 1000000/frames_to_render) ÷ 2
    frames = [] 
    for i in 1:frames_to_render
        frame = Gray.(reshape(s[i*itval:1_000_000 + itval*i - 1], (1000, 1000)))
        # frame[250:749, 250:749] = center
        frame = map(clamp01nan, frame)
        frame = convert.(RGB{N0f8}, frame)
        # dil_ero!(center)
        center = morphogradient(center)
        push!(frames, frame)
    end
    frames
end

function morpho_vid_pixs(pixs, frames_to_render)
    itval = round(Int, 1000000/frames_to_render) ÷ 2
    frames = [] 
    for i in 1:frames_to_render
        frame = reshape(pixs[i*itval:1_000_000 + itval*i - 1], (1000, 1000))
        frame = map(clamp01nan, frame)
        frame = convert.(RGB{N0f8}, frame)
        push!(frames, frame)
    end
    frames
end

lows = Bandpass(1, 500; fs=sr)
mids = Bandpass(500, 2_000; fs=sr)
highs = Bandpass(2_000, 10_000; fs=sr)

designmethod = Butterworth(4)
filt_sigs = []
resps = [lows, mids, highs]
for rt in resps
    f = digitalfilter(rt, designmethod)
    filt_s = filt(f, s)
    push!(filt_sigs, filt_s)
end
l, m, h = filt_sigs

l,m,h = SippyArt.wav_to_img.([l,m,h])

function eqs_to_pixs(l, m, h)
    pixs = RGB{N0f8}[]
    for i in eachindex(l)
        pix = RGB(l[i], m[i], h[i])
        push!(pixs, convert(RGB{N0f8}, clamp01nan(pix)))
    end
    pixs
end
Vector{Matrix{RGB{N0f8}}}

# RGB(l[10000], m[10000], h[10000])
i = 10000
pix = RGB(l[i], m[i], h[i])

pixs = eqs_to_pixs(filt_sigs...)
frame = reshape(pixs[1:1_000_000], (1000, 1000))
frames = morpho_vid_pixs(pixs, 240)

alg = Equalization(nbins = 7)
img_adjusted = adjust_histogram(img, alg)


new_frames = []
for i in 1:240
    # alg = Equalization(nbins = i % 10 + 2)
    alg = Equalization(nbins = i ÷ 10 + 2)
    push!(new_frames, adjust_histogram(frames[i], alg))
end
# my_videowrite("data/2021-09-05/outputs/contrast_adj_81_3.mp4", "data/2021-09-05/outputs/tmp3/", new_frames)

img = rand(RGB{N0f8}, 1000, 1000)
img2 = imrotate(img, π/2)
img2 = img2[1:1000, 1:1000]
img3 = similar(img)

@benchmark blend.(img, img2, mode=BlendDarken)
@benchmark rand(100)

p = img[1]
p2 = img2[2]

@benchmark blend(p, p2, mode=BlendDarken)

v = VideoIO.load("data/2021-09-05/outputs/contrast_adj.mp4")
v2 = VideoIO.load("data/2021-09-05/outputs/contrast_adj2.mp4")
v2 = map(x->imrotate(x, π/2)[1:1000, 1:1000], v2)
v = Vector{eltype(v2)}(v)
typeof(v)
typeof(v2)
size(v2[1]) == size(v[1])

fs = SippyArt.blend_vids(v, v2, BlendDifference)



bw = Gray.(rand(Bool, 1000, 1000))
bw = tophat(bw)

imgg = imfilter(bw, Kernel.gaussian(3))
img = imfilter(img, Kernel.Laplacian())


d = dither(d, BalancedCenteredPoint())
d = dither(d, Bayer())
d = dither(bw, Bayer())
d = dither(bw, SimpleErrorDiffusion())


img = rand(RGB{N0f8}, 1000, 1000) - rand(RGB{N0f8}, 1000, 1000)
white = RGB{Float32}(1, 1, 1)
yellow = RGB{Float32}(1, 1, 0)
green = RGB{Float32}(0, 0.5, 0)
orange = RGB{Float32}(1, 0.5, 0)
red = RGB{Float32}(1, 0, 0)
blue = RGB{Float32}(0, 0, 1)

rubiks_colors = [white, yellow, green, orange, red, blue]
d = dither(cimg, FloydSteinberg(), rubiks_colors)
d = dither(d, WhiteNoiseThreshold(), rubiks_colors)


f(img) = imfilter(img, Kernel.Laplacian())
f(img) = imfilter(img, Kernel.gaussian(5))


img = rand(RGB, 1000, 1000) # - rand(RGB, 1000, 1000)

# VideoIO.save("video.mp4", fs, framerate=30)

img = map(clamp01nan, img)
img = convert.(RGB{N0f8}, img)


img = rand(RGB{N0f8}, 1000, 1000) # - rand(RGB{N0f8}, 1000, 1000)
# center = morphogradient(center)

using Images, FileIO, SippyArt
img = Gray.(rand(Bool, 1000, 1000)) # - rand(RGB{N0f8}, 1000, 1000)
img = rand(Gray, 1000, 1000) # - rand(RGB{N0f8}, 1000, 1000)
fs = my_apply(img -> imfilter(morphogradient(img), Kernel.gaussian(5)) , img; N=600)
fs = my_apply(img -> imfilter(img, Kernel.gaussian(5)) , img; N=600)
my_videowrite("data/2021-09-05/outputs/grad_rand3.mp4", "data/2021-09-05/outputs/tmp2/", fs)