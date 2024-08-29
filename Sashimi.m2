needsPackage "MultigradedImplicitization"

-- This is the parameterization of the d-th secant power of the Segre embedding of P^(a-1) x P^(b-1) x P^(c-1)
-- The case when d = a = b = c = 4 corresponds to the Salmon problem
secantSegre = (d, a, b, c) -> (
    S := QQ[x_(0,0) .. x_(d-1,a-1), y_(0,0) .. y_(d-1,b-1), z_(0,0) .. z_(d-1,c-1)];
    R := QQ[p_(0,0,0) .. p_(a-1, b-1, c-1)];
    images := toList for l in (0,0,0)..(a-1, b-1, c-1) list sum apply(d, i -> x_(i, l_0) * y_(i, l_1) * z_(i, l_2));
    map(S, R, images))

-- An easier alternate example is the case when d = a = b = c = 3 which corresponds to the 3 state general Markov model
-- Running the code below takes about 3 seconds to compute the degree 4 phylogenetic invariants of this model
d = 3;
a = b = c = 3;

F = secantSegre(d, a, b, c);
S = source F
R = target F

end--
restart
needs "Sashimi.m2"

G = time componentsOfKernel(4, F);
