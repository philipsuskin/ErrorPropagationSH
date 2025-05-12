using TypedPolynomials
using SphericalHarmonicExpansions
include("plotMagneticField.jl")

using MonteCarloMeasurements

L = 6
R = 0.045
center = [0.0, 0.0, 0.0]
N = 1000

cx = cy = cz = SphericalHarmonicExpansions.SphericalHarmonicCoefficients(L, R, true)

# Assuming variance of coefficients (σ^2_{γ})
function case1(variance = 1)
  for l in 0:L
    for m in -l:l
      varianceCoeff = variance

      for cs in [cx, cy, cz]
        cs[l, m] = varianceCoeff
      end
    end
  end

  return "varianceCoeff-$(variance).png"
end

# Assuming variance of solid harmonics (Z_l^m)
function case2(variance = 1, particleCount = 100)
  Bik = 1
  Zlm = Particles(particleCount, Gamma(variance))
  for l in 0:L
    for m in -l:l
      pCoeff = (2l + 1) / (R^l * N) * sum(Bik .* Zlm)

      for cs in [cx, cy, cz]
        cs[l, m] = rand(pCoeff.particles)
      end
    end
  end

  return "varianceSolidHarmonic-$(variance).png"
end

# figname = case1()
figname = case2()

TypedPolynomials.@polyvar x y z

Bx = sphericalHarmonicsExpansion(cx, x, y, z)
By = sphericalHarmonicsExpansion(cy, x, y, z)
Bz = sphericalHarmonicsExpansion(cz, x, y, z)

fx(ix, iy, iz) = Bx(x => ix, y => iy, z => iz)
fy(ix, iy, iz) = By(x => ix, y => iy, z => iz)
fz(ix, iy, iz) = Bz(x => ix, y => iy, z => iz)

fignumber = plotMagneticField([fx, fy, fz], R, center=center)
savefig(figname)

fieldError = (ix, iy, iz) -> (fx(ix, iy, iz), fy(ix, iy, iz), fz(ix, iy, iz))
println("Field error: $(fieldError(center...))")