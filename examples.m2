load "MultigradedImplicitization.m2"

-- Example 1
R = QQ[x,y,z]
S = QQ[s,t]

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im)

D = maxGrading(f)

G = {}
L = for i in 0..60 list (
    G = G | componentOfIdeal({i}, G, f, D);
    {i, G}
);
netList L


--Example 2

load "~/Documents/my_scripts/sunlets/sunletQuadGens.m2"

n = 6;
R = qRing n;
f = sunletParam n;
S = ring(f#0);
f = map(S,R,f);

D = matrix{toList(2^(n-1):1)};

G = {};
L = for i in 0..3 list (
    G = G | componentOfIdeal({i}, G, f, D);
    {i, #G}
);
netList L
