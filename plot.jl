using TypedPolynomials
using SphericalHarmonicExpansions
include("plotMagneticField.jl")

TypedPolynomials.@polyvar x y z
TypedPolynomials.@polyvar c[1:49]

varianceCoeff = 1

L = 6

cx = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cy = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cz = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)

for l in 0:L
  for m in -l:l
    cx[l, m] = varianceCoeff
    cy[l, m] = varianceCoeff
    cz[l, m] = varianceCoeff
  end
end

if false
  Bx = sphericalHarmonicsExpansion(cx, x, y, z)
  By = sphericalHarmonicsExpansion(cx, x, y, z)
  Bz = sphericalHarmonicsExpansion(cx, x, y, z)

  fx(ix, iy, iz) = Bx(x => ix, y => iy, z => iz)
  fy(ix, iy, iz) = By(x => ix, y => iy, z => iz)
  fz(ix, iy, iz) = Bz(x => ix, y => iy, z => iz)

  fignumber = plotMagneticField([fx, fy, fz], 0.045, center=[0.0, 0.0, 0.0])
  savefig("varianceCoeff-$(varianceCoeff).png")

  fieldError = (ix, iy, iz) -> (fx(ix, iy, iz), fy(ix, iy, iz), fz(ix, iy, iz))
  println("Field error: $(fieldError(0.0, 0.0, 0.0))")
else
  using MPISphericalHarmonics, MPIUI
  coeffs = MagneticFieldCoefficients([cx; cy; cz;;], 0.045)
  MagneticFieldViewer(coeffs)
end