restart
needsPackage "Graphs";
needsPackage "Polyhedra";
load "MultigradedImplicitization.m2"

-- Creates the ring of parameters whose variables are the root parameters and entries of the transition matrices
-- N, a phylogenetic network
-- sym1 is used for root parameters
-- sym2 is used for transition matrix entries
gmmParamRing = (N, sym1, sym2) -> (

	indSet := flatten for e in edges(N) list(


		flatten for i from 0 to 3 list for j from 0 to 3 list (toList(e), i, j)
		);

	r := getSymbol sym1; 
	m := getSymbol sym2;

	return R := QQ[apply(4, i -> r_i) | apply(indSet, i -> m_i)];
)



intVertices = T -> delete(null, vertices(T) / (i -> if degree(T, i) > 1 then i))


gmmTreeParam = T -> (

	int := sort intVertices(T);
	L := sort leaves(graph(graph(T)));
	n := #L;
	S := gmmParamRing(T, "r", "m");

	phi := for leafState in toList((n:0)..(n:3)) list(


			sum for intState in toList(toList(#int:0)..toList(#int:3)) list(


				states := hashTable(apply(L, i -> i => leafState_(i-1)) | apply(int, i -> i => intState_(i - #L - 1)));

				(r_(intState_0))*product(apply(edges(T), e -> (m_(e, states#(e_0), states#(e_1)))))
				)
			);

	return phi
);


T = digraph {{4,1}, {4,2}, {4,3}}
phiList = gmmTreeParam(T);
f = map(ring(phiList_0), QQ[q_(0,0,0)..q_(3,3,3)], phiList);

componentsOfIdeal(f,5)




A = transpose mingens image transpose (maxGrading phi)_(toList(0..63))






B = matrix({toList(64:1)}) || A^(toList(0..(numrows(A)-2)))
C = coneFromVData(B)
I = -1*halfspaces(C)
v = transpose matrix {toList(numrows(I):0)}
E = matrix({{1} | toList(numcols(I)-1:0)}) || hyperplanes(C)
d = 5
w = transpose matrix {{d} | toList(numrows(E)-1:0)}
P = polyhedronFromHData(I, v, E, w)
lats = latticePoints(P);














