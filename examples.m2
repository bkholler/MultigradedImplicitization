load "MultigradedImplicitization.m2"

-- Example 1
R = QQ[x,y,z]
S = QQ[s,t]

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im)

D = maxGrading(f)

G = {}
L = for i in 0..100 list (
    G = G | componentOfIdeal({i}, G, f, D);
    {i, G}
);
netList L


--Example 2

load "~/Documents/my_scripts/sunlets/sunletQuadGens.m2"

R = qRing(5)
f = sunletParam 5
S = ring(f#0)
f = map(S,R,f)

D = maxGrading f
deg = D * (transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0}})
deg = flatten entries deg
G = {}

componentOfIdeal(deg, G, f, D)