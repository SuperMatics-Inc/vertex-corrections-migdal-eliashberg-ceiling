"""Reproduce the small-x and high-x certificates. Intermediate misses stay unresolved."""
from fractions import Fraction as Q
from functools import lru_cache

from math import isqrt
from mpmath import iv

iv.dps = 43
LAM = tuple(map(Q, ('1/9', '1/4', '1/2', '1', '2')))
ETA = (Q(1), Q(3, 2), None)
GRID = tuple(map(Q, ('1/16', '1/8', '1/4', '1/2', '1', '2', '4', '8', '16')))
H = tuple(Q(j, 100) for j in range(2, 51, 2))
PI_LOWER = Q(314159, 100000)
C_LOWER = Q('0.18272624777232950')
C_UPPER = Q('0.18272624777232952')
C_BOUND = iv.mpf(['0.18272624777232950', '0.18272624777232952'])
VOLUME = Q(15, 32)


def box(q):
    q = Q(q)
    return iv.mpf(q.numerator) / q.denominator


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def decimal(q, places, upper=True):
    """Format exact rationals or binary64 values with directed decimal rounding."""
    scale = 10**places
    units = -((-Q(q)*scale)//1) if upper else (Q(q)*scale)//1
    return f'{units//scale}.{units%scale:0{places}d}'


def small_guard(cutoff=40, mutant=None):
    # Exact rational Schur check, plus an interval evaluation of the normal and form guards.
    gamma = Q(13, 5)
    require((2 * PI_LOWER * C_LOWER)**2 * 2 > gamma, 'Einstein envelope')
    def major(j):
        return 2 / (1 + gamma*j*j)
    rows = [sum((major(n+m+1) + (major(n-m) if m != n else 0)) / (2*m+1)
                for m in range(64)) for n in range(16)]
    tail = Q(1, 63**2) if mutant != 'missing-small-tail' else 0
    fixed_rows = max(rows) + tail
    require(max(fixed_rows, 22 / (16 * gamma)) < Q(19, 20), 'small-x Schur')
    # Independent interval construction, with a longer cutoff and the Hankel/Toeplitz sums split.
    rating = []
    for n in range(16):
        k = lambda j: 2 / (1 + box(gamma)*j*j)
        toeplitz = sum(k(n-m) / (2*m+1) for m in range(96) if m != n)
        hankel = sum(k(n+m+1) / (2*m+1) for m in range(96))
        rating.append(toeplitz + hankel + box(Q(1, 95**2)))
    require(all(v < iv.mpf('0.95') for v in rating), 'independent Schur')
    require(all(rating[n] < box(rows[n] + tail) for n in range(16)), 'small cutoff disagreement')

    x = box(Q(1, 2**cutoff))
    r = box(Q(23, 5))
    require(2 + 1 / (12 * C_BOUND**2) < r, 'transfer total')
    j = 3 + 4 * iv.ln(3/x)
    l = 2 * r*r * x * j**3 / (3 * iv.pi**2)
    eps = 2 * (1 + 8*r) * l
    contraction1 = 4*r*x/iv.pi + 32*r*l/iv.pi + box(Q(1, 4))
    contraction2 = 6*l*(1 + 16*r)
    image = r*x/iv.pi + 4*r*l
    margin = 1 - box(Q(19, 20)) - 2*r*eps - 2*x*(r+r*r)/iv.pi - 10*l
    width = 1 + r + (1 + 8*r)*l
    require(j > 12 and l < box(Q(1, 8)) and (1 + 8*r)*l < box(Q(1, 8)), 'normal box')
    require(image < box(Q(1, 4)) and 4*r*x/iv.pi < box(Q(1, 2)), 'rainbow image')
    require(contraction1 < 1 and contraction2 < 1 and margin > 0, 'small-x guards')
    # Independent exact rational upper bounds for logs and pi, and a conservative lower margin.
    q = Q(1, 2**cutoff)
    ll = 2*Q(23, 5)**2*q*(3+4*(Q(11, 10)+cutoff*Q(7, 10)))**3/(3*PI_LOWER**2)
    require(l < box(ll) and iv.ln(3) < box(Q(11,10)) and iv.ln(2) < box(Q(7,10)), 'independent logs')
    require(6*ll*(1+16*Q(23,5)) < 1, 'independent contraction')
    rr = Q(23,5)
    lower = (1-Q(19,20)-4*rr*(1+8*rr)*ll-2*q*(rr+rr*rr)/PI_LOWER-10*ll)/6
    require(lower > Q('0.00807313') and width < 6, 'published small-x spectral gap')
    print(f'SMALL x<=2^-{cutoff}: L<={decimal(ll,14)}, raw spectral gap>={decimal(lower,11,False)}')
    return lower


@lru_cache(None)
def exact_center(lam, eta, h, p, count):
    if eta is None:
        return (1-h)**(-p)
    a = (2 * PI_LOWER * C_LOWER)**2 * lam * eta**2
    def term(j):
        return 1 / (1 + a*j*j)
    finite = sum((term(m) + term(m+1))/(2*m+1-h)**p for m in range(count))
    extra = 2 / (a * 2**p * (p+1) * (count-1)**(p+1))
    return finite + extra


@lru_cache(None)
def independent_center(lam, eta, h, p, count):
    if eta is None:
        return 1 / (1 - box(h))**p
    s = (2 * iv.pi * C_BOUND)**2 * box(lam * eta**2)
    total = iv.mpf(0)
    for k in range(-count, count):
        w = abs(2*k+1)
        total += 1 / ((1+s*k*k) * (w-box(h))**p)
    # One independent one-sided transfer-tail bound for each half-line.
    tpos = 1 / (s * 2**p * (p+1) * (count-1)**(p+1))
    tneg = 1 / (s * 2**p * (p+1) * (count-1)**(p+1))
    return total + tpos + tneg


def guards(high, low, r, h, mass=VOLUME):
    # The same derivative majorant bounds the full-momentum absolute pairing Schur row.
    z = 1/(1-h)
    return r*high + mass*r*r*high*high*z, r*low + mass*r*r*(2*high*low*z + high*high*z*z)


def independent_guards(lam, eta, x, h):
    j1 = independent_center(lam, eta, h, 1, 43)
    j2 = independent_center(lam, eta, h, 2, 43)
    a = 2*box(lam/x) / iv.pi
    u = 1-box(h)
    f = a*j1 + (15*box(Q(1,32))) * (a*j1)**2/u
    cross = (a*j1)**2/(u*u) + 2*(a*j1/u)*(a*j2)
    b = a*j2 + 15*box(Q(1,32))*cross
    return f, b


def dos_guards(lam, eta, x, h, independent):
    """A second box using the bare complex mean and spatial L2 rather than peak lines."""
    n = 16 if independent else 12
    if independent:
        X,H0 = box(x),box(h)
        alpha = (2*iv.pi*C_BOUND)**2 * box(lam*(eta or 1)**2)
        factor = 2*box(lam/x)/iv.pi
        root = iv.sqrt
        angle = lambda u: iv.atan2(1,u)
        weights = range(-n,n)
        def k(j):
            return 1/(1+alpha*j*j)
    else:
        X,H0 = x,h
        alpha = (2*PI_LOWER*C_LOWER)**2 * lam*(eta or 1)**2
        factor = 2*lam/(PI_LOWER*x)
        scale = 2**28
        def root(t):
            return Q(1+isqrt(t.numerator*scale**2//t.denominator),scale)
        def angle(u):
            return min(1/u,Q(3927,2500),Q(3927,2500)-u+u**3/3)
        weights = range(n)
        def k(m):
            return 1/(1+alpha*m*m) + 1/(1+alpha*(m+1)**2)
    def line(w):
        z = X*angle(X*w)/(w-H0)
        return w*z, root(w*z/(w-H0)), w*z/(w-H0)
    z = line(1)
    if eta is None:
        s1,s2,s3 = z
    else:
        rows = [tuple(k(j)*v for v in line(abs(2*j+1))) for j in weights]
        s1,s2,s3 = (sum(v) for v in zip(*rows))
        t = 1/(2*alpha*(n-1)**2)
        s1,s2,s3 = s1+t,s2+t,s3+1/(6*alpha*(n-1)**3)
    f = factor*s1 + (factor*s2)**2*z[1]
    b = factor*s3 + factor**2*(2*s3*s2*z[1]+s2**2*z[2])
    return f,b


def choose_dos(lam, eta, x):
    winner = None
    for h in H:
        f,b = dos_guards(lam,eta,x,h,False)
        if f < h and b < 1 and (winner is None or b < winner[0]):
            winner = b,h
    if winner is None:
        return None
    b,h = winner
    f2,b2 = dos_guards(lam,eta,x,h,True)
    require(f2 < box(h) and b2 < box(b), 'independent DOS box')
    return b,h


def choose(lam, eta, x):
    r = 2*lam/(PI_LOWER*x)
    winner = None
    for h in H:
        t1,t2 = (exact_center(lam, eta, h, p, 34) for p in (1,2))
        f, b = guards(t1, t2, r, h)
        if f < h and b < 1 and (winner is None or b < winner[0]):
            winner = b, h
    if winner is None:
        return None
    b, h = winner
    f2, b2 = independent_guards(lam, eta, x, h)
    require(f2 < box(h) and b2 < box(b), 'independent high-x bound')
    return b,h


def mutants():
    # Volume: independent slice-by-slice cap integral, versus the translated lens polynomial.
    # Two equal unit caps at distance d have height 1-d/2 and volume ratio (3 h^2-h^3)/2.
    coeff = (Q(1), -Q(3,4), Q(0), Q(1,16))
    radial = 3 * sum(a / (i+3) for i,a in enumerate(coeff))
    import sympy as sp
    d = sp.Symbol('d')
    h = 1-d/2
    slice_integral = sp.integrate(3*d*d*(3*h*h-h**3)/2, (d,0,1))
    require(radial == slice_integral == VOLUME, 'overlap geometry')
    rejected = 0
    reflected = independent_center(Q(2), Q(1), Q(1,4), 1, 43)
    folded = exact_center(Q(2), Q(1), Q(1,4), 1, 34)
    require(reflected < box(folded), 'reflection fixture agreement')
    r, h0, z = Q(1,5), Q(1,4), Q(4,3)
    # Resolvent derivative dc/dSigma=c^2 at c=4/3; extract the exact linear coefficient.
    t = sp.Symbol('t')
    derivative = (r*(z+t*z*z)+VOLUME*r*r*(z+t*z*z)**3).expand().coeff(t)
    full = guards(z,z*z,r,h0)[1]
    require(full == derivative, 'three insertion derivative fixture')
    fixtures = (
        ('lost overlap', Q(1) != slice_integral),
        ('wrong reflected tail', box(folded - Q(1,100)) < reflected),
        ('missing middle insertion', full - VOLUME*r*r*z**4 < derivative),
        ('wrong Einstein normalization', (2*PI_LOWER*C_LOWER)**2*2/4 < Q(13,5)),
    )
    for name, caught in fixtures:
        require(caught, name + ' mutant survived')
        rejected += 1
    try:
        small_guard(mutant='missing-small-tail')
    except ArithmeticError as error:
        require(str(error) == 'small cutoff disagreement', 'wrong mutant failure')
        rejected += 1
    else:
        raise ArithmeticError('small-tail mutant survived')
    print(f'MUTANTS rejected={rejected}')


def weak_guard():
    r, h = Q(3,25), Q(3,10)
    j = 6-Q(22,15)*(1-Q(3,5))
    require(Q(3,5)**4 > Q(1,8) and j*j < 3*PI_LOWER**2, 'Hölder DOS constant')
    image = r/(1-h)+3*r*r/(1-h)**3
    derivative = r/(1-h)**2+9*r*r/(1-h)**4
    require(image < h and derivative < Q(157,200), 'weak Young ball')
    require(Q(7,60)*(1+Q(10,768)) < r and Q(7,60)+1/(12*C_LOWER**2*30**2)<r,
            'weak finite transfer bounds')
    # Independent interval endpoints keep the secant slope and pi rather than rounding to 3.
    s = iv.mpf(2)**(-iv.mpf(3)/4)
    c = ((6-box(Q(22,15))*(1-s))/iv.pi)**2
    a = 1-box(h)
    ri = box(r)
    require(c < 3 and ri/a+c*ri**2/a**3 < box(h) and ri/a**2+3*c*ri**2/a**4 < box(derivative),
            'weak interval cross-check')
    print(f'WEAK Young R<=3/25 image<={decimal(image,8)} derivative<={decimal(derivative,8)}')


def main():
    improved = 0
    small_guard()
    weak_guard()
    counts = dict.fromkeys(ETA, 0)
    anchors = {Q(1): tuple(map(Q, (1,2,4,8,16))), Q(3,2): tuple(map(Q, (1,2,4,8,8)))}
    print('lambda,eta,x,guard,h,interval_upper')
    for i,lam in enumerate(LAM):
        for eta in ETA:
            for x in GRID:
                witness = choose(lam, eta, x)
                if eta in anchors and x == anchors[eta][i]:
                    require(witness is not None, f'peak rectangle anchor {lam},{eta},{x}')
                kind = 'PEAK'
                if witness is None:
                    witness = choose_dos(lam,eta,x)
                    kind = 'DOS'
                    improved += witness is not None
                row = f'{lam},{eta if eta is not None else "static"},{x}'
                if witness is None:
                    print(f'{row},UNRESOLVED,,')
                    continue
                b,h = witness
                require(h in H and 0 < b < 1, 'invalid box witness')
                print(f'{row},{kind},{h},{decimal(b,8)}')
                counts[eta] += 1
    # The high-x criterion only improves with increasing x and increasing eta.
    # For continuous lambda <=2 both a K factors and a^2 K K increase with lambda.
    uniform = choose(Q(2), Q(1), Q(19,2))
    require(uniform is not None and uniform[1] == Q(9,25) and uniform[0] < Q('0.70468692'),
            'published uniform large-x witness')
    require(counts == {Q(1):15,Q(3,2):16,None:23} and improved == 3, 'published slice counts')
    print(f'GRID certified={sum(counts.values())}/135, slices={list(counts.values())}; '
          f'uniform x>=19/2 PEAK h={uniform[1]} upper<={decimal(uniform[0],8)}')
    mutants()


if __name__ == '__main__':
    main()
