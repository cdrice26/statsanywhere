"""
Generate reference_values.txt using SciPy. Exits with non-zero if SciPy is not available.

Usage: python3 generate_reference.py
Writes tests/reference_values.txt containing lines of the form:
<key> <args...> <value>

This must be run on a machine with SciPy installed to produce verified reference values.
"""
import sys
try:
    import scipy.stats as stats
    import scipy.special as spec
except Exception as e:
    sys.stderr.write('SciPy is required to generate reference values: ' + str(e) + '\n')
    sys.exit(2)

out_path = 'tests/reference_values.txt'
with open(out_path, 'w') as out:
    # binomial
    out.write(f"binompmf 5 3 0.7 {stats.binom.pmf(3,5,0.7):.12g}\n")
    out.write(f"binomcdf 6 2 0.7 {stats.binom.cdf(2,6,0.7):.12g}\n")
    out.write(f"binompmf 100 50 0.7 {stats.binom.pmf(50,100,0.7):.12g}\n")
    out.write(f"binomcdf 100 50 0.7 {stats.binom.cdf(50,100,0.7):.12g}\n")

    # chi-squared
    out.write(f"X2pdf 5 3 {stats.chi2.pdf(5,3):.12g}\n")
    out.write(f"X2cdf 5 3 {stats.chi2.cdf(5,3):.12g}\n")

    # F distribution
    out.write(f"fpdf 3 5 4 {stats.f.pdf(3,5,4):.12g}\n")
    out.write(f"fcdf 3 5 4 {stats.f.cdf(3,5,4):.12g}\n")

    # normal
    out.write(f"normalcdf_std -1 1 {(stats.norm.cdf(1) - stats.norm.cdf(-1)):.12g}\n")
    out.write(f"normalcdf -1 1 1 2 {(stats.norm.cdf((1-1)/2) - stats.norm.cdf((-1-1)/2)):.12g}\n")

    # poisson
    out.write(f"poissonpmf 3 5 {stats.poisson.pmf(3,5):.12g}\n")
    out.write(f"poissoncdf 3 5 {stats.poisson.cdf(3,5):.12g}\n")

    # t distribution
    out.write(f"tpdf 3 5 {stats.t.pdf(3,5):.12g}\n")
    out.write(f"tcdf 3 5 {stats.t.cdf(3,5):.12g}\n")
    out.write(f"tcdf 0.5 5 {stats.t.cdf(0.5,5):.12g}\n")

    # operations
    out.write(f"factorial 9 {float(__import__('math').factorial(9)):.12g}\n")
    out.write(f"combination 6 4 {float(__import__('math').comb(6,4)):.12g}\n")

    # special functions via scipy.special
    out.write(f"gamma 5.0 {spec.gamma(5.0):.12g}\n")
    out.write(f"gamma 6.7 {spec.gamma(6.7):.12g}\n")

    # SciPy's gammainc returns the regularized lower incomplete gamma P(a,x).
    # Unregularized lower incomplete gamma = P(a,x) * Gamma(a)
    out.write(f"lower_incomplete_gamma 1 5 {(spec.gammainc(1,5) * spec.gamma(1)):.12g}\n")
    out.write(f"lower_incomplete_gamma 3 5.5 {(spec.gammainc(3,5.5) * spec.gamma(3)):.12g}\n")

    out.write(f"beta 3 2 {spec.beta(3,2):.12g}\n")
    out.write(f"beta 4 10 {spec.beta(4,10):.12g}\n")

    # incomplete beta (unregularized): B_x(a,b) = betainc(a,b,x) * beta(a,b)
    bx_1 = spec.betainc(1, 3, 0.4) * spec.beta(1, 3)
    out.write(f"incomplete_beta 0.4 1 3 {bx_1:.12g}\n")
    bx_2 = spec.betainc(2, 3, 1.0) * spec.beta(2, 3)
    out.write(f"incomplete_beta 1 2 3 {bx_2:.12g}\n")

    # regularized incomplete beta (I_x)
    out.write(f"regularized_incomplete_beta 0.4 1 3 {spec.betainc(1,3,0.4):.12g}\n")
    out.write(f"regularized_incomplete_beta 1 2 3 {spec.betainc(2,3,1.0):.12g}\n")

    out.write(f"erf 1 {spec.erf(1):.12g}\n")
    out.write(f"erf 0.5 {spec.erf(0.5):.12g}\n")

print('Wrote', out_path)
