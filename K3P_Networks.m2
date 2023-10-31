restart
installPackage "MultigradedImplicitization"
loadPackage("PhylogeneticTrees", FileName => "/Users/bkholler/My Drive/MultigradedImplicitization/PhylogeneticTrees.m2")

sunletParam = (n, M) -> (

	indR := leafColorings(n,M);
	GH := groupHash(M);
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



M = K3Pmodel
R = qRing(4, M)
images = sunletParam(4, M);
phi = map(ring images_0, R, images);

G = time componentsOfKernel(3, phi)

dom = newRing(source phi, Degrees => maxGrading phi);
B = sub(basis(3, source phi), dom);
lats = unique apply(flatten entries B, m -> degree m);