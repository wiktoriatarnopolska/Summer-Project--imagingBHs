using Gradus
import Gradus.propernorm

@inline function tetradtest(g::AbstractMatrix{T}, v) where {T}
    # normalise the vector to unit length
    vt = v ./ √abs(propernorm(g, v))

    if v[1] > 0 && v[2] == 0 && v[3] == 0 && v[4] > 0
        vϕ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 0, 0, 1), (vt,), g)
        vr = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 0, 0), (vt,vϕ), g)
        vθ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 0, 1, 0), (vt,vϕ,vr), g)
        (vt, vϕ, vr, vθ)
    elseif v[1] > 0 && v[2] > 0 && v[3] == 0 && v[4] > 0
        vϕ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 0, 1), (vt,), g)
        vr = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 0, 0, 1), (vt, vϕ), g)
        vθ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 0, 1, 0), (vt, vϕ, vr), g)
        (vt, vϕ, vr, vθ)
    elseif v[1] > 0 && v[2] > 0 && v[3] > 0 && v[4] > 0
        vϕ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 0, 0), (vt,), g)
        vr = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 1, 0), (vt, vϕ), g)
        vθ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 1, 1), (vt, vϕ, vr), g)
        (vt, vϕ, vr, vθ)
    elseif v[1] > 0 && v[2] > 0 && v[3] == 0 && v[4] == 0
        vϕ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 0, 0), (vt,), g)
        vr = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 1, 0), (vt, vϕ), g)
        vθ = Gradus.GradusBase.gramschmidt(SVector{4,T}(1, 1, 1, 1), (vt, vϕ, vr), g)
        (vt, vϕ, vr, vθ)
    else
        print("error")
    end
end

tetradtest(m::AbstractMetric, x, v) = tetradtest(Gradus.metric(m, x), v)

export tetradtest