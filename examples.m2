load "MultigradedImplicitization.m2"

-- Example 1
R = QQ[x,y,z]
S = QQ[s,t]

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im)

D = maxGrading(f)
L = for i in 30..100 list (
    {i, componentOfIdeal({i}, f, MaxGrading => D)}
);
netList L

componentOfIdeal({30},f, MaxGrading => D)

restrictedComponentOfIdeal({60}, {z^5 - y^3 - x^2}, f, D)