using DynamicPolynomials
using SphericalHarmonicExpansions
include("plotMagneticField.jl")

@polyvar x y z
@polyvar c[1:49]
polynomialGlobal = sphericalHarmonicsExpansion(c, x, y, z, true)

L = 6

cx = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cy = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)
cz = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, 0.045, true)

for l in 0:L
  for m in -l:l
    coeff = SphericalHarmonicExpansions.zlm(l, m, x, y, z)
    println("l = $l, m = $m")
    println("coeff = $coeff")
    cx[l, m] = coeff
    cy[l, m] = coeff
    cz[l, m] = coeff
  end
end

# using MPISphericalHarmonics, MPIUI
# coeffs = MagneticFieldCoefficients([cx; cy; cz;;], 0.045)
# MagneticFieldViewer(coeffs)

Bx = polynomialGlobal(c => vec(cx.coefficients))
By = polynomialGlobal(c => vec(cy.coefficients))
Bz = polynomialGlobal(c => vec(cz.coefficients))

fx(ix, iy, iz) = Bx(x => ix, y => iy, z => iz)
fy(ix, iy, iz) = By(x => ix, y => iy, z => iz)
fz(ix, iy, iz) = Bz(x => ix, y => iy, z => iz)

fignumber = plotMagneticField([fx, fy, fz], 0.045, center=[0.0, 0.0, 0.0])
savefig("magneticField.png")