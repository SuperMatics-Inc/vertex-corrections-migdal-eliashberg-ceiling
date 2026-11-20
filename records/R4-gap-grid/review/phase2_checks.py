"""Phase 2 review: execute local source mutants and audit unpadded quantization."""
from fractions import Fraction as F
from math import ceil, nextafter, sqrt
from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / 'G1-GAP'
sys.path.insert(0, str(SRC))

from finite import PUBLISHED, Tube
import finite_check as folded
from replay import C_BOUND, ETA, LAM, box, iv


def check(condition, what):
    if not condition:
        raise ArithmeticError(what)


def patch_check(old, new):
    path = SRC / 'finite_check.py'
    source = path.read_text(encoding='utf-8')
    check(source.count(old) == 1, f'mutation anchor not unique: {old!r}')
    namespace = {}
    exec(compile(source.replace(old, new, 1), str(path), 'exec'), namespace)
    return namespace['Check']


def interpolation_mutant():
    bad = patch_check('error = enclosure(dot(self.fold_hi,modulus)/8)', 'error = 0.')
    lam, eta, x = F(2), F(1), F(1,4)
    tube = Tube(lam, eta, x)
    d0, d1 = (folded.np.pad(d[:tube.n], (0, 9+tube.max_m-tube.n)) for d in tube.centers(1))
    args = (0., 1., d0, d1, 2*folded.np.arange(len(d0))+1)
    good_e = folded.Check(lam,eta,x,96).rainbow_error(*args)[0]
    bad_e = bad(lam,eta,x,96).rainbow_error(*args)[0]
    # Exact prescribed midpoint widths with an interval recomputation of the signed row.
    a = 2*box(lam/x)/iv.pi
    b = (2*iv.pi*C_BOUND)**2*box(lam*eta**2)
    k = lambda j: a/(1+b*j*j)
    t = box(F(1,2))
    X = box(x)
    partial = sum((k(m)-k(m+1))*X*iv.atan2(1, X*(2*m+1+t*box(F(d1[m]))))
                  for m in range(len(d1)))
    midpoint = t*(partial+iv.mpf([0, 1])*(2*a/(b*(len(d1)-1)*(2*len(d1)+1)))) \
               - t*box(F(d1[0]))
    check(abs(midpoint) < box(F(good_e)) and abs(midpoint) > box(F(bad_e)),
          'curvature mutant not distinguished by exact midpoint')
    print(f'killed zero-curvature mutant: true midpoint |defect|={abs(midpoint)}; '
          f'good bound={good_e:.9g}, mutant bound={bad_e:.9g}')


def disconnected_tube_mutant():
    bad = patch_check('require(np.array_equal(prior,d0) and d1.shape==(self.n,)',
                      'require(d1.shape==(self.n,)')
    lam, eta, x = F(1,9), F(1), F(1,2)
    result, why = Tube(lam,eta,x).verify(48)
    check(result is not None, why)
    witness = result[-1]
    verifier = folded.Check(lam,eta,x,96)
    check(verifier.validate(witness)[0] is not None, 'unmutated witness')
    broken = list(witness)
    a,b,r,p = broken[10]
    a = a.copy()
    a[0] += 2**-20
    broken[10] = (a,b,r,p)
    try:
        verifier.validate(broken)
    except ArithmeticError as e:
        check(str(e)=='inhomogeneous endpoints or invalid radius', 'unexpected failure')
    else:
        raise ArithmeticError('bad junction survived original verifier')
    check(bad(lam,eta,x,96).validate(broken)[0] is not None,
          'broken witness failed for a reason other than junction check')
    print('killed missing-junction mutant: discontinuous centers are rejected only with the check')


def quantization_audit():
    den = folded.DEN
    # The quotient at line 202 is rounded to nearest before the integer ceiling.
    # A representable center can lie below the exact root of a quantization boundary.
    w, x = 1, F(1,8)
    for k in range(57457, 64527):
        target = den/k - float(1/x**2)
        float_v = sqrt(target)
        for v in (float_v,nextafter(float_v,0.)):
            exact = 1/(F(v)**2+1/x**2)
            used = float(folded.quantize(folded.np.array([1/(v*v+64)]),True)[0])/den
            if F(used) < exact:
                check(1 < v < 3, 'adversarial center outside test range')
                vp = folded.enclosure(v,-1)
                padded = folded.enclosure(1/(vp*vp+64))
                check(F(padded) >= exact, 'pad did not enclose inverse moment')
                print(f'unguarded quotient rounds inward: w={w}, x={x}, d={v-w:.17g}, '
                      f'k={k}, exact*DEN-k={float(exact*den-k):.9g}')
                break
        else:
            continue
        break
    else:
        raise ArithmeticError('expected inward-quantization adversary not found')

    repaired = patch_check('j = self.weight/(v*v+1/self.x**2)',
                           'v_lower = enclosure(v,-1)\n'
                           '                j = enclosure(self.weight/(v_lower*v_lower+1/self.x**2))')
    count = 0
    worst = 0.
    for eta, rows in zip(ETA[:2],PUBLISHED):
        for lam, xs in zip(LAM,rows):
            for x in map(F,xs):
                result, why = Tube(lam,eta,x).verify(48)
                check(result is not None, why)
                report, why = folded.Check(lam,eta,x,96).validate(result[-1])
                check(report is not None, why)
                ref, why = repaired(lam,eta,x,96).validate(result[-1])
                check(ref is not None, why)
                worst = max(worst, result[0], ref[0])
                for d0,d1,_,pre in result[-1]:
                    if not pre:
                        continue
                    for i,d in enumerate(folded.np.minimum(d0,d1)):
                        w = 2*i+1
                        weight = 1 if w < 17 else 4
                        v = F(w)+F(d)
                        exact = weight/(v*v+1/x**2)
                        vf = w+d
                        j = weight/(vf*vf+1/float(x)**2)
                        check(ceil(j*den) >= ceil(exact*den), 'published quotient quantized inward')
                        count += 1
                # The zero-center extra frequencies have exact integer widths.
                for i in range(96,297):
                    w = 2*i+1
                    exact = 4/(F(w*w)+1/x**2)
                    j = 4/(float(w*w)+1/float(x)**2)
                    check(ceil(j*den) >= ceil(exact*den), 'published tail quotient quantized inward')
                    count += 1
    check(F(worst) < F('0.50417617'), 'padding broke published pairing target')
    print(f'published inverse moments: {count} exact ceiling comparisons passed; '
          'padded verifier accepted all 15 with bound < 0.50417617')


if __name__ == '__main__':
    interpolation_mutant()
    disconnected_tube_mutant()
    quantization_audit()
