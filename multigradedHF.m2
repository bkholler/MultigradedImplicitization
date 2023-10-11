needsPackage "Polyhedra"

sortLatPts = (lats, D, A) -> (
    zer := matrix toList(numRows(lats#0) : {0});
    N := numColumns(D);
    steps := apply(N, i -> D_{i});
    unsorted := delete(zer, lats);
    sorted := {zer};
    while unsorted != {} do (
	for v in unsorted do (
	    toSortNext = all(N, i -> member(v - steps#i, sorted) or any(flatten entries(A*(v - steps#i)), j-> j<0));
	    if toSortNext then (
		sorted = sorted | {v};
		unsorted = delete(v,unsorted);
		)
	    )
	);
    return sorted
    );


--returns all the degrees needed to compute D_pt
--pt = the 'max' degree trying to compute
--D = columns are degrees of variables
--A = rows are inward point normals of cone spanned by D
findLatPts = (pt,D,A) -> (
    I := -A || A;
    m := numRows(A);
    b := map(ZZ^m, ZZ^1, 0) || (A*pt);
    P := polyhedronFromHData(I,b);
    lats := latticePoints(P);
    return sortLatPts(lats,D,A)
    );



--antiDiff map with respect to X in degree d
--Phi : D_d --> D_{d - deg(X)}
antiDiff = (d, R, X) -> (
    --ind := last baseName(X);
    shift := -degree(X);
    deg := flatten entries d;
    colIndex := flatten entries basis(deg,R);
    rowIndex := flatten entries basis(deg + shift,R); 
    return matrix (for a in rowIndex list (
	for b in colIndex list (
	    if a == b // X then (
		1
		) else (
		0
		)
	    )
	)) ** QQ;
    );

coeff = e -> product apply(e, i -> 1/i!);

eval = R -> (
    n := numColumns(vars R);
    return map(QQ,R,toList(n:0))
    );

codom = (v, D, A) -> (
    n := numColumns(D);
    posVals := apply(n, i -> A*(v - D_{i}));
    L := {};
    for i from 0 to n-1 do (
	if all(flatten(entries(posVals#i)), b -> b>=0) then (
	    L = L|{i};
	    );
	); 
    return L
    );

--d = weight of closedness subspace trying to compute
--B = D_(d - D_{i})(I) for all suitable i
--R = QQ[x_0..x_n] is the ambient ring
--D = ith column is degree(x_i)
--A = rows are inward point normals of cone spanned by D
closedSubspace = (d, B, R, D, A) -> (
    deg := flatten entries d;
    N := numColumns(basis(deg, R));
    shifts := {};
    for i from 0 to -1 + numColumns(D) do (
	if all(flatten(entries(A*(d - D_{i}))), j -> j >= 0) then (
	    shifts = shifts|{i};
	    );
	);
    Phi := apply(shifts, i -> antiDiff(d, R, x_i));
    Cd := gens intersect(apply(#B, i -> ker(map(coker(B#i), QQ^N, Phi#i))));
    return Cd;
    );

dualSpaceInDegreeD = (d, G, R, C) -> (
    ed := flatten((flatten entries basis(flatten entries d, R)) / exponents);
    dSpace := image C;
    ev := eval(R);
    for g in G do (
	if degree(g) == flatten(entries(d)) then (
	    dSpace = intersect(dSpace, ker matrix{apply(ed, e -> (ev(diff(coeff(e)*R_e,g))))});
	    );
	);
    return gens dSpace;	
    );



dualSpace = (v, G, R, D, A) -> (
    LL := findLatPts(v,D,A);
    MD := new MutableHashTable;
    scan(LL, w -> MD#w = matrix{{1}}**QQ);
    time for i from 1 to -1+#LL do (
	ind := codom(LL#i,D,A);
	B := apply(ind, j -> MD#(LL#i - D_{j}));
	C := closedSubspace(LL#i,B,R,D,A);
	MD#(LL#i) = dualSpaceInDegreeD(LL#i,G,R,C);
	--print((i+1)/#LL);
	);
    return new HashTable from MD;
    );


dSpace = (v, G, D) -> (
    A := facets(coneFromVData(D));
    S := ring(ideal(G));
    n := numColumns vars S;
    R := QQ[x_0..x_(n-1), Degrees => entries transpose D]; -- change back to x's if it doesn't work
    f := map(R,S, flatten entries vars R);
    H := G / f;
    return dualSpace(v,H,R,D,A);
    );

HF = (i, b, G, D) -> (
    d := dSpace(i*b,G,D);
    mHF :=  applyPairs(d, (m,n) -> (flatten entries m, numColumns n));
    return apply(i+1, j -> mHF#(j*(flatten entries b)));
    );


    








