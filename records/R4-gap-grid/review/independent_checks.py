"""Independent rational checks; imports neither of the submitted verifiers."""
from fractions import Fraction as F

C_MIN = F('0.18272624777232950')  # GAMMA certified interval, rounded down
LAM = (F(1, 9), F(1, 4), F(1, 2), F(1), F(2))
XS = tuple(F(2) ** j for j in range(-4, 5))
SCALE = 10**26


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def round_up(t):
    return F(-(-t.numerator * SCALE // t.denominator), SCALE)


def decimal_bound(t, places, upper):
    scale = 10**places
    num = t.numerator * scale
    units = -(-num // t.denominator) if upper else num // t.denominator
    return f'{units // scale}.{units % scale:0{places}d}'


def atan_pair(t, terms=24):
    # Even alternating partial sum is below atan, next partial sum above.
    check(0 <= t <= 1 and terms % 2 == 0, 'atan domain')
    low = sum((-1)**j * t**(2*j+1) / (2*j+1) for j in range(terms))
    return low, low + t**(2*terms+1) / (2*terms+1)


A5, B5 = atan_pair(F(1, 5))
A239, B239 = atan_pair(F(1, 239))
PI_MIN = 16*A5 - 4*B239  # Machin's identity
PI_MAX = 16*B5 - 4*A239


def log_upper(a):
    t = (a-1) / (a+1)
    n = 48
    return 2 * (sum(t**(2*j+1) / (2*j+1) for j in range(n))
                + t**(2*n+1) / ((2*n+1)*(1-t*t)))


def small_x():
    r, gamma, ncut = F(23, 5), F(13, 5), 128
    check(2 + 1/(12*C_MIN**2) < r, 'uniform transfer sum')
    check(2*(2*PI_MIN*C_MIN)**2 > gamma, 'transfer envelope')
    k = lambda j: 2/(1+gamma*j*j)
    tail = F(25, 9)/(2*gamma*(ncut-1)**2)
    rows = [sum(round_up((k(n+m+1) + (k(n-m) if n != m else 0))/(2*m+1))
                for m in range(ncut)) + tail for n in range(16)]
    check(max(*rows, F(22)/(16*gamma)) < F(19, 20), 'canceled Schur rows')
    x = F(1, 2**40)
    logs = log_upper(F(3)) + 40*log_upper(F(2))
    l = 2*r*r*x*(3+4*logs)**3/(3*PI_MIN**2)
    eps = 2*(1+8*r)*l
    check(max(l, (1+8*r)*l) < F(1, 8), 'line box')
    check(r*x/PI_MIN + 4*r*l < F(1, 4) and 4*r*x/PI_MIN < F(1, 2), 'image')
    check(4*r*x/PI_MIN + 32*r*l/PI_MIN + F(1, 4) < 1, 'rainbow contraction')
    check(6*l*(1+16*r) < 1 and 1+r+(1+8*r)*l < 6, 'crossed contraction and width')
    gap = (F(1, 20) - 2*r*eps - 2*x*(r+r*r)/PI_MIN - 10*l)/6
    check(gap > F('0.00807313'), 'small-x published gap')
    loose_l = 2*r*r*x*(3+4*(F(11, 10)+40*F(7, 10)))**3/(3*F('3.14159')**2)
    loose_gap = (F(1, 20) - 4*r*(1+8*r)*loose_l
                 - 2*x*(r+r*r)/F('3.14159') - 10*loose_l)/6
    print(f'small x: Schur row <= {decimal_bound(max(rows), 9, True)}; '
          f'gap >= {decimal_bound(gap, 12, False)}; '
          f'submitted rational gap >= {decimal_bound(loose_gap, 15, False)}')
    check(l < F('2.21278157e-06') < loose_l, 'rounded-down proxy still bounds actual L')
    check(F('0.00807313972') <= loose_gap, 'submitted lower print direction')
    print('small x: printed L rounds its rational proxy down, but bounds L independently')


def volume():
    # Interchange the separation d and axial cap coordinate y:
    # V = (3/2) integral_0^1 (1-y^2) min(1,2y)^3 dy.
    v = F(3, 2)*(8*(F(1, 2)**4/4-F(1, 2)**6/6) + F(2, 3) - F(1, 2) + F(1, 2)**3/3)
    check(v == F(15,32), 'unit-ball triple overlap')
    return v


def peak(lam, eta, x, h, ncut=80):
    alpha = (2*PI_MIN*C_MIN)**2*lam*eta*eta
    k = lambda j: 1/(1+alpha*j*j)
    a = 2*lam/(PI_MIN*x)
    # Direct unfolded rows; omitted pairs are bounded by m^(-p-2).
    s, t = (sum(round_up(k(j)/(abs(2*j+1)-h)**p) for j in range(-ncut,ncut))
            + 2/(alpha*2**p*(p+1)*(ncut-1)**(p+1)) for p in (1,2))
    z = 1/(1-h)
    cs = a*s
    ct = a*t
    return cs+V*cs**2*z, ct+V*(2*cs*ct*z+(cs*z)**2)


def witnesses(test, lam, eta, x):
    for n in range(1,26):
        h = F(n,50)
        f, b = test(lam,eta,x,h)
        if f < h and b < 1:
            yield b, h, f


def high_x_and_grid():
    high, h, f = min(witnesses(peak,F(2),F(1),F(19,2)))
    check(f < h and high < F('0.70468692'), 'published large-x bound')
    print(f'large x: h={h}, image <= {decimal_bound(f, 12, True)}, '
          f'row <= {decimal_bound(high, 12, True)}; V={V}')
    for eta, thresholds in ((F(1),(1,2,4,8,16)),(F(3,2),(1,2,4,8,8))):
        cnt = 0
        for lam, x0 in zip(LAM, thresholds):
            x = F(x0)
            witness = min(witnesses(peak,lam,eta,x),default=None)
            check(witness is not None, 'finite Einstein peak anchor')
            cnt += sum(y >= x for y in XS)
            print(f'anchor eta={eta}, lambda={lam}, x={x}: peak, h={witness[1]}, '
                  f'row <= {decimal_bound(witness[0], 9, True)}')
        check(cnt == (15 if eta == 1 else 16), 'finite Einstein count')


def near_static_budget():
    errs = []
    for lam in LAM:
        for x in XS:
            delta = 1/(12*C_MIN**2*lam*(2**18)**2)
            g = 2*lam*x/PI_MIN
            zmax = 2/x
            omitted = delta*(2+delta)
            map_error = (g*delta*zmax + V*g*g*omitted*zmax**3)/x
            jac_error = g*delta*zmax**2 + 3*V*g*g*omitted*zmax**4
            check(map_error < F(1,2**17), 'near-static normal remainder')
            errs.append((map_error,jac_error))
    res, jac = (max(col) for col in zip(*errs))
    print(f'near-static: max residual/x <= {decimal_bound(res, 15, True)}; '
          f'max derivative remainder <= {decimal_bound(jac, 15, True)}')


if __name__ == '__main__':
    V = volume()
    small_x()
    high_x_and_grid()
    near_static_budget()
