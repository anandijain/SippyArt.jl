using WAV, DataInterpolations 
using BenchmarkTools

w, sr = wavread("/Users/anand/.julia/dev/SippyArt/data/2021-09-05/sources/082121.wav")


function _slow_nx_2(w, n)
    arr = Matrix{eltype(w)}(undef, size(w, 1) * n, size(w, 2))
    for (i, r) in enumerate(eachrow(w))
        for j in 1:n
            arr[n * i + j - n, :] .= r
        end 
    end
    arr
end



function resample(itp, ts)
    permutedims(reduce(hcat, itp.(ts)))
end

function resample2!(w, itp, ts)

    w .= permutedims(reduce(hcat, itp.(ts)))
    nothing
end


function resample3!(w, itp, ts)
    rs = itp.(ts)
    for (i,r) in enumerate(rs)
        w[i, :] .= r
    end
    nothing
end

function resample2(itp, ts)
    w = Matrix{eltype(first(itp.u))}(undef, length(ts), 2)
    resample2!(w, itp, ts)
    w
end


rows = map(Array, eachrow(w))
N = length(rows)
ratio = 2//3
ts = 1:ratio:N

wnew = Matrix{eltype(first(itp.u))}(undef, length(ts), 2)

Ts = [LinearInterpolation, QuadraticInterpolation, QuadraticSpline, CubicSpline]
for T in Ts[1:2]
    fn = "/Users/anand/.julia/dev/SippyArt/data/2021-09-05/sources/itps/itp_$(string(T))_$(ratio.num)_$(ratio.den).wav"
    itp = T(rows, 1:N)
    resample3!(wnew, itp, ts)
    wavwrite(wnew, fn; Fs=sr)
end


itp = LagrangeInterpolation(rows, 1:N)