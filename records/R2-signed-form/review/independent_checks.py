"""Independent finite algebra and rational-guard review of G1-SIGNED.

No imports from the author replays. The continuum estimates are checked on
paper in REVIEW.md; symbolic polynomials do not certify the continuum proof.
"""
from fractions import Fraction as F
from sympy import I, Matrix, expand, symbols


def audit_crossed_identity():
    r, q, y = [symbols(f'{prefix}:4', real=True) for prefix in 'rqy']
    e = symbols('e', real=True)
    s3 = Matrix.diag(1, -1)
    g = [Matrix([[r[i]-I*q[i], -e*q[i]*y[i]],
                 [-e*q[i]*y[i], -r[i]-I*q[i]]]) for i in range(4)]

    def anchor_rows(middle_no_hole=False):
        # Corners are indexed around the parallelogram; using the opposite
        # traversal of both transfers takes the real part automatically.
        normal_rows = []
        anomalous_rows = []
        for p in range(4):
            a, b, c = ((p+i) % 4 for i in (1, 2, 3))
            middle = g[b] if not middle_no_hole else g[b].conjugate()
            cross = s3*g[a]*s3*middle*s3*g[c]*s3
            normal = -expand(cross[0, 0].subs(e, 0)).as_real_imag()[1]
            linear = expand(cross[0, 1]).coeff(e, 1).as_real_imag()[0]
            normal_rows.append(q[p]*y[p]**2*normal)
            anomalous_rows.append(q[p]*y[p]*linear)
        return sum(normal_rows)/4, sum(anomalous_rows)/4

    normal, anomalous = anchor_rows()
    polynomial = expand(normal-anomalous)
    alt = sum((-1)**i*y[i] for i in range(4))
    # Each undirected pair occurs at two of the four anchors.
    pair = sum(q[i]*q[j]*(y[i]-y[j])**2 * r[k]*r[l]/4
               for i in range(4) for j in range(i+1, 4)
               for k, l in [[k for k in range(4) if k not in (i, j)]])
    expected = -q[0]*q[1]*q[2]*q[3]*alt**2/4 + pair
    assert polynomial == expand(expected)
    assert expand(polynomial + anomalous) != 0  # omitted normal subtraction
    _, mutated = anchor_rows(middle_no_hole=True)
    assert polynomial != expand(normal-mutated)
    assert polynomial != expand(expected-pair)
    # One static transfer and a momentum-independent trial: paired frequency
    # corners are p=a and b=c. The result is a nonzero real-line pair term.
    duplicate = {r[1]:r[0], q[1]:q[0], y[1]:y[0],
                 r[3]:r[2], q[3]:q[2], y[3]:y[2]}
    one_static = expand(polynomial.subs(duplicate))
    assert one_static == expand(q[0]*q[2]*r[0]*r[2]*(y[0]-y[2])**2)
    assert one_static != 0
    # A coincident-frequency corner is still counted with multiplicity.
    assert (2 - (-1) + 2 - (-1))**2 > 4*(2**2+(-1)**2+(-1)**2)


def audit_rational_guards():
    delta = F(1, 65536)
    # Direct witness domain: two independent uniform transfers, both in
    # [1/2,3/4]. The other vertices miss [-1/16,1/16] by interval bounds.
    edge = F(1, 16)
    assert edge-F(1, 2) < -edge and edge-1 < -edge
    reduction = (F(3, 4)-F(1, 2))**2 / 4
    negative = reduction / 64   # quarter times four lower |Q| bounds
    real_error = 6*16*(16*delta)**2
    tail_bound = F(1, 1 << 24)
    # Compare to the stated 3s bound, and separately use the exact 2s+s².
    margin = negative-real_error-160*3*tail_bound
    assert margin == F(55, 262144) > F(1, 8192)
    assert negative-real_error-160*(2*tail_bound+tail_bound**2) > margin
    assert negative-real_error-160*3*tail_bound*8 < F(1, 8192)
    assert negative/16-real_error-160*3*tail_bound < 0

    # From the ordered one-static sum with I(1,1)/π=1/16,
    # I(1,3)/π=1/128, and the two choices of static transfer.
    one_static = 2*(-8*F(1,16)+4*F(1,128))
    assert one_static == -F(15,16)
    assert one_static + F(10*25,1000) < -F(1,4)
    assert (-one_static) + F(10*25,1000) > 0

    # Monotonic scalar polynomial guards, evaluated at their upper edges.
    r = F(1,16)
    assert 2*r+8*r*r == F(5,32) < F(1,2)
    assert 4*r+48*r*r == F(7,16) < 1
    R = F(1,128)
    assert F(1, 1<<20)+4*F(1, 1<<10) < R
    assert F(2*4096, 3**3) < 304
    L = 304*R*R
    assert 3*R+L == F(43,1024) < F(1,2)
    assert 4*R+6*L == F(73,512) < 1
    # Exact comparisons to already-certified decimal bounds on C_infty.
    lo, hi, printed = map(F, ('0.1827262477723295163523252',
                             '0.1827262477723295176168720',
                             '0.18272624777232952'))
    assert F(1,8) < lo < hi < printed
    assert lo*(1 << 16) > (1 << 12)


if __name__ == '__main__':
    audit_crossed_identity()
    audit_rational_guards()
    print('Independent exact polynomial, one-static limit, and rational guards agree.')
