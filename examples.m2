
-- Example 1
-- this example is broken at the moment. i think it's because of the non-standard Z-grading
restart
load "MultigradedImplicitization.m2"
R = QQ[x,y,z];
S = QQ[s,t];

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im);

D = maxGrading(f);


G = {}

for i in 0..100 do (
    G = G | componentOfIdeal(i, G, f, R, new MutableHashTable)
)



--Example 2
--still fine
restart
load "MultigradedImplicitization.m2"
load "~/Documents/my_scripts/sunlets/sunletQuadGens.m2"

n = 5;
R = qRing n;
f = sunletParam n;
S = ring(f#0);
f = map(S,R,f);

componentsOfIdeal(f,3)

D = matrix{toList(2^(n-1):1)};







--Example 3
restart 
load "MultigradedImplicitization.m2"
needsPackage "Polyhedra"

R = QQ[x..z]
S = QQ[u,v]

f1 = map(S,R,{u^2,u*v,v^2})
D1 = maxGrading f1
I1 = ker f1

f2 = map(S,R,{(u+v)^2,(u+v)*(u-v),(u-v)^2})
D2 = maxGrading f2
I2 = ker f2

I1 == I2
D1 != D2

isSubset(I1, toricMarkov(D2,R))




