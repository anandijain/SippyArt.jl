module SippyArt
using WAV, Images
using ColorBlendModes, ColorTypes
using Dates

datadir = joinpath(@__DIR__, "../data/")
outdir = joinpath(datadir, "outs")

function new_project()
    p = joinpath(datadir, string(today()))

    mkpath.([
        joinpath(p, "sources"),
        joinpath(p, "outputs")
    ])
end

function my_videowrite(vidname, imgs_outdir, imgs)
    for (i, img) in enumerate(imgs)
        save(joinpath(imgs_outdir, "$i.jpeg"), img)
    end
    run(Cmd(["ffmpeg", "-framerate", "60", "-start_number", "1", "-i", "$(imgs_outdir)/%d.jpeg", "-r", "60", "-y", "$(vidname)"]))
end

function my_apply!(frames, f, img)
    frames[1] = img
    for i in 2:length(frames)
        frames[i] = f(frames[i-1])
    end
end
function my_apply(f, img; N=50)
    frames = Vector{typeof(img)}(undef, N)
    my_apply!(frames, f, img)
    frames
end

"WAV encodes in bounds (-1, 1)"
function wav_to_img(s)
    (s .+ 1) ./ 2
end

function clamp_convert(img)
    convert.(RGB{N0f8}, clamp01nan!(img))
end

function clamp_convert_pix(pix)
    convert(RGB{N0f8}, clamp01nan(pix))
end

function blend_vids!(out, v, v2, mode=BlendHardLight)
    for i in eachindex(v)
        out[i] = blend.(v[i], v2[i], mode=mode)
    end
    nothing
end

function blend_vids(v, v2, mode=BlendHardLight)
    out = similar(v)
    blend_vids!(out, v, v2, mode)
    out
end

export datadir
export new_project, my_videowrite
export my_apply, my_apply!

end # module
