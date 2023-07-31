using Gradus
using Plots

# defining new type -- coronal spectrum (photon index (power of g))
abstract type AbstractCoronalSpectrum end

coronal_spectrum(spectrum::AbstractCoronalSpectrum, g) = error("not implemented for $(typeof(spectrum))")

# creating the power law -> g^Γ
struct PowerLawSpectrum{T} <: AbstractCoronalSpectrum
    Γ::T
end

function coronal_spectrum(spectrum::PowerLawSpectrum, g)
    g^(spectrum.Γ)
end

# now replacing g^2 with function returning g^Γ
function source_to_disc_emissivity(m::AbstractStaticAxisSymmetric, N, A, x, g, spectrum::AbstractCoronalSpectrum)
    v = CircularOrbits.fourvelocity(m, SVector(x[2], x[3]))
    # account for relativistic effects in area due to lorentz shift
    γ = lorentz_factor(m, x, v)
    # divide by area to get number density
    N / (coronal_spectrum(spectrum, g) * A * γ)
end

point_source_equitorial_disc_emissivity(θ, g, A, γ, spectrum) = sin(θ) / (coronal_spectrum(spectrum, g) * A * γ)

# adding argument to the emissivity profile function + radial disc  profile
function Gradus.emissivity_profile(
    sampler::AbstractDirectionSampler,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::Gradus.AbstractCoronaModel,
    spectrum::AbstractCoronalSpectrum;
    grid = GeometricGrid(),
    N = 100,
    kwargs...,
)

    RadialDiscProfile(
        tracecorona(m, d, model, spectrum; sampler = sampler, kwargs...);
        grid = grid,
        N = N,
    )
end

Gradus.emissivity_profile(
    ::Nothing,
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    spectrum::AbstractCoronalSpectrum,
    model::LampPostModel;
    kwargs...,
) = Gradus._point_source_symmetric_emissivity_profile(m, d, model; kwargs...)



m = KerrMetric(1.0, 0.998)
d = GeometricThinDisc(0.0, 100.0, π/2)

model1 = LampPostModel(h = 10.0)

# first case -- g^2 spectrum
em_prof1 = Gradus.emissivity_profile(
m, 
d, 
model1, 
n_samples = 100_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile1 = Gradus.RadialDiscProfile(em_prof1)

plot(profile1, linecolor = "orchid3", title = "Disc profile Γ = 2.0")

# introducing Γ = 3.0
spectrum = PowerLawSpectrum(3.0)

model2 = LampPostModel(h = 10.0)

# now with spectrum as an argument -- returns error
em_prof2 = Gradus.emissivity_profile(
spectrum,
m, 
d, 
model2, 
n_samples = 100_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

profile2 = Gradus.RadialDiscProfile(em_prof2)
plot!(profile2, linecolor = "royalblue4", title = "Disc profile Γ = 3.0")