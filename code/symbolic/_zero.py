"""Helper for the symbolic checks of Paper C: decide whether a sympy expression vanishes identically.
First tries sympy simplification; if that is inconclusive the expression is evaluated with mpmath at 50 digits at several
random points (reported explicitly as 'numerical-zero at N random points', never silently)."""
import random
import sympy as sp

def is_zero(expr, symbols=(), positive=True, seed=1, npts=6):
    e = sp.simplify(expr)
    if e == 0:
        return True, "exactly 0 (sympy.simplify)"
    rnd = random.Random(seed)
    worst = 0
    for _ in range(npts):
        sub = {s: sp.Rational(rnd.randint(11, 97), rnd.randint(7, 31)) for s in symbols}
        v = sp.N(e.subs(sub), 50)
        worst = max(worst, abs(v))
    return worst < sp.Float('1e-40'), f"simplify inconclusive; |value| <= {sp.N(worst, 3)} at {npts} random rational points (50 digits)"

def report(tag, ok, how):
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {how}", flush=True)
    return ok
