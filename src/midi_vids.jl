#=
Goal here is maybe eventually to have a more MIDI like video editor

Starting by having each note in a midi clip trigger a video.
=#

using MIDI, DataFrames
using VideoIO, Colors
using Images, FileIO

"assuming 4/4"
function set_qpm!(m, qpm)
    microseconds_per_quarter_note = round(Int, inv(qpm//(60*1_000_000)))
    e = MIDI.SetTempoEvent(0, 0x51, microseconds_per_quarter_note)
    map(t->addevent!(t, 0, e), m.tracks)
    nothing
end

m = readMIDIFile("data/2021-09-05/sources/basic_drums.mid")
m = readMIDIFile("data/2021-09-05/sources/basic_melody.mid")

set_qpm!(m, 90)


ps = unique(df.pitch) # number of clips to use 

vidpaths = joinpath.("/Users/anand/Movies/", ["81.mp4", "green.mp4", "fluid.mp4"])
vids = collect.(VideoIO.load.(vidpaths))
firsts = first.(vids)
imgs = Matrix{RGBA{N0f8}}.(firsts)

# want to find out for each frame, which clips are on
ts = collect(seconds_per_frame .* (0:240)) # in seconds
tms = 1000*ts 

n_frames = 400
# to start it might be easier to make framerate = tickrate
ticks_per_frame = 1
# ticks_per_frame = 24 # ? 

pitch_to_clip = Dict(ps .=> vids)

ms_t = ms_per_tick(m)
ts = collect(ms_t .* (0:n_frames)) # in seconds


function build_frame!(canvas, pitchmap, ps)
    for p in ps
        f = pitchmap[p]
        h, w = size(f)
        canvas[1:h, 1:w] = f
    end
    canvas
end

"df = DataFrame(getnotes(midi))"
function drums_to_clip(df, pitchmap, outdir)
    canvas = Matrix{RGBA{N0f8}}(undef, 1000, 1000)
    df[!, :note_end] .=  df.position .+ df.duration
    n_frames = maximum(df.note_end)
    for i in 1:n_frames
        sdf = @view df[df.position .<= i .<= (df.position .+ df.duration), :]
        ps = sort(unique(sdf.pitch)) # frames to show 
        fill!(canvas, zero(RGBA))
        build_frame!(canvas, pitchmap, ps)
        save(joinpath(outdir, "$i.jpeg"), canvas)
    end
end

function get_pitchmap(pitches)
    # red = fill(RGBA(1, 0, 0, 1), (1000,1000))
    red = imgs[1][1:end, 1:1000]
    blue = fill(RGBA(0, 1, 0, 1), (500, 500))
    green = fill(RGBA(0, 0, 1, 1), (250, 250))
    xx = fill(RGBA(1, 0, 1, 1), (150, 150))
    xxx = fill(RGBA(0, 1, 1, 1), (1000, 100))
    aa = imgs[2][1:end, 1:1000]
    bb = imgs[3][1:1000, 1:1000]
    Dict(pitches .=> [red, blue, green, xx, xxx, aa, bb])
end


function doit(fn, qpm, outdir)
    m = readMIDIFile(fn)
    set_qpm!(m, qpm)
    ns = getnotes(m).notes
    df = DataFrame(ns)
    pitchmap = get_pitchmap(unique(df.pitch))
    drums_to_clip(df, pitchmap, outdir)
end


audiopath = "data/2021-09-05/sources/basic_melody.wav"
outdir = "data/2021-09-05/outputs/tmp4/"
vidpath = "data/2021-09-05/outputs/midi6.mp4"
mymidiwrite(vidpath, outdir, audiopath; fps=144)

fn = "data/2021-09-05/sources/basic_melody.mid"
myqpm = 90
doit(fn, myqpm, outdir)
run(Cmd(["ffmpeg", "-framerate", "144", "-start_number", "1", "-i", "$(outdir)/%d.jpeg", "-i", audiopath, "-r", "144", "-y", "$(vidpath)"]))
run(Cmd(["ffmpeg", "-framerate", "60", "-start_number", "1", "-i", "$(outdir)/%d.jpeg", "-r", "60", "-y", "$(vidpath)"]))
using WAV

w, sr = wavread(audiopath)
seconds = length(w) / sr
n_frames = Int(maximum(df.note_end))
n_frames/seconds