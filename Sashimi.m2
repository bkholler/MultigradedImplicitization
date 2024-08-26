restart
needsPackage "MultigradedImplicitization"

-- This is the parameterization of the d-th secant power of the Segre embedding of P^(a-1) x P^(b-1) x P^(c-1)
-- The case when d = a = b = c = 4 corresponds to the Salmon problem
secantSegre = (d, a, b, c) -> (

	S := QQ[toList(for i from 0 to d-1 list x_(i, 0)..x_(i, a-1)) | toList(for i from 0 to d-1 list y_(i, 0)..y_(i, b-1)) | toList(for i from 0 to d-1 list z_(i, 0)..z_(i, c-1))];
	R := QQ[p_(0,0,0)..p_(a-1, b-1, c-1)];

	images := toList for l in (0,0,0)..(a-1, b-1, c-1) list sum apply(d, i -> x_(i, l_0)*y_(i, l_1)*z_(i, l_2));

	map(S, R, images)
	)


-- An easier alternate example is the case when d = a = b = c = 3 which corresponds to the 3 state general Markov model
-- Running the code below takes about 3 seconds to compute the degree 4 phylogenetic invariants of this model
l = 3;
a = 3;
b = 3;
c = 3;

F = secantSegre(l, a, b, c)
G = time componentsOfKernel(4, F)