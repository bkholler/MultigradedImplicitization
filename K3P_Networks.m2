restart

-- This loads our new package which contains our implementation of Algorithm 1
needsPackage "MultigradedImplicitization"

-- This loads the PhylogeneticTrees package for macaulay2
needsPackage "PhylogeneticTrees"


-- n, the number of leaves
-- M, a model, see PhylogeneticTrees.m2
-- Outputs the parameterization of a cyclically labeled n-leaf sunlet network with the leaf labeled 1 as the reticulation vertex
-- The output is given as a list whose entries correspond to the fourier coordinates in lexicographic order on the label sequence
sunletParam = (n, M) -> (

	-- This makes 
	indR := leafColorings(n,M);
	GH := new MutableHashTable from apply(group(M), toList(0..#group(M)-1), (i, j) -> i => j);

	-- This makes the ring of parameters which is the codomain of phi
	-- the "a" parameters correspond to the leaves of the network while the "b" parameters correspond to the internal edges
	indS := flatten for i from 1 to n list apply(group(M), j -> {i, GH#j});
	a := symbol a;
	b := symbol b;
	S := QQ[apply(indS, k -> a_k) | apply(indS, k -> b_k)];
	

	images := for g in indR list(

		aProd := product(apply(n, i -> a_{i+1, GH#(g_i)}));
		bProd1 := product(for j from 0 to #g - 2 list b_{j+1, GH#(sum(g_(toList(0..j))))});
    	bProd2 := product(for j from 1 to #g- 1 list b_{j+1, GH#(sum(g_(toList(1..j))))});
    	aProd*(bProd1 + bProd2)
		);

	return images
	)

-- Running this block of code below takes about 8 minutes
-- n is the number of leaves on the network and d is the total degree up to which one wants to compute the kernel
n = 4;
d = 3;
M = K3Pmodel;
R = qRing(n, M);
images = sunletParam(n, M);
phi = map(ring images_0, R, images);

G = time componentsOfKernel(d, phi)

