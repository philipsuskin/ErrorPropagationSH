using TypedPolynomials
using SphericalHarmonicExpansions
include("plotMagneticField.jl")

TypedPolynomials.@polyvar x y z

varianceCoeff = 1

L = 6

cx = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cy = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cz = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)

for l in 0:L
  for m in -l:l
    for cs in [cx, cy, cz]
      cs[l, m] = varianceCoeff
    end
  end
end

using MPISphericalHarmonics, MPIUI
coeffs = MagneticFieldCoefficients([cx; cy; cz;;], 0.045)
MagneticFieldViewer(coeffs)