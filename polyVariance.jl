using TypedPolynomials
using SphericalHarmonicExpansions

TypedPolynomials.@polyvar x y z

l = 5
m = 2
Zₗᵐ = SphericalHarmonicExpansions.zlm(l, m, x, y, z)

println("Zₗᵐ: $Zₗᵐ")

# Variances
σ²ᵣ = [0.01, 0.01, 0.01]  # [σ²ₓ, σ²ᵧ, σ²𝓏]

cov_xy = 0.0
cov_xz = 0.0
cov_yz = 0.0

# Derivatives
δZₗᵐ = [differentiate(Zₗᵐ, var) for var in (x, y, z)]

# Variance contribution
σ² = sum(δZₗᵐ[i]^2 * σ²ᵣ[i] for i in 1:3)

# Covariance contribution
σ² += 2 * (δZₗᵐ[1]*δZₗᵐ[2]*cov_xy + δZₗᵐ[1]*δZₗᵐ[3]*cov_xz + δZₗᵐ[2]*δZₗᵐ[3]*cov_yz)

println("Variance: $σ²")