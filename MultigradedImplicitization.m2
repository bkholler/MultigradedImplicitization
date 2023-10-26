needsPackage "gfanInterface";
allowableThreads = maxAllowableThreads;

-- Input : f : dom = kk[x_1..x_n] --> codom = kk[y_1..y_m]
-- Output : lineality space of Groebner fan
maxGrading = phi -> (
    -- set up elimination ideal
    dom := source phi;
    codom := target phi;
    elimRing := dom ** codom;
    X := vars dom;
    n := numgens dom;
    elimIdeal := ideal(sub(X, elimRing) - sub(phi(X), elimRing));
    return (transpose linealitySpace(gfanHomogeneitySpace(elimIdeal)))_(toList(0..n-1))
);

findBasisInDegree = (deg, dom, previousPolys, basisHash) -> (

    if #previousPolys == 0 then (
        return basis(deg, dom, Strategy => Default);
    );

    -- otherwise, we shift G in all possible ways to land in R_deg

    L := for g in previousPolys list(

            checkDegree := deg - degree(g);

            if basisHash#?checkDegree then flatten entries (g*basisHash#checkDegree) else g*(flatten entries basis(deg - degree(g), dom, Strategy => Default))
        );
    
    -- stick em all in a matrix
    mat := sub(matrix({flatten L}), dom);

    -- and collect coefficients.
    (mons, coeffs) := coefficients(mat);

    -- find the independent linear relations
    coeffs = try mingens(image sub(coeffs, QQ)) else print(deg);

    -- remove monomials corresponding to pivots
    badMonomials := apply(pivots coeffs, i -> mons_(0,i#0));
    monomialBasis := flatten entries basis(deg, dom);
    scan(badMonomials, m -> monomialBasis = delete(m, monomialBasis));

    return matrix{monomialBasis}
);

-- G = known generators of degree less than deg
-- this guy fucks with all Z^k-gradings
componentOfIdeal = (deg, dom, phi, previousPolys, basisHash) -> (

    -- The span of monomialBasis is all you need to search for, since G is known
    monomialBasis := basisHash#deg;

    -- collect coefficients into a matrix
    (mons, coeffs) := try coefficients(phi(sub(monomialBasis, source phi))) else print("THIS FAILED");

    -- find the linear relations among coefficients
    K := gens ker sub(coeffs, QQ);

    newGens := flatten entries (monomialBasis * K);

    return newGens
);

-- this guy only works if it has usual Z-grading
-- things to add to this function:
--      make it work when it doesn't have usual Z-grading, and check for this
--      make an option to override this maybe if you already have D that you want to use
--      make an option if the ring already has the appropriate grading
--      write it in Oscar lol
componentsOfIdeal = (phi, d) -> (

    n := numgens(source phi);
    omega := maxGrading(phi);
    dom := newRing(source phi, Degrees => omega);

    previousPolys := {};
    basisHash := new MutableHashTable;
    
    -- assumes homogeneous with normal Z-grading
    for i in 1..d do (
        B := sub(basis(i, source phi), dom);
        lats := unique apply(flatten entries B, m -> degree m);

        scan(lats, deg -> basisHash#deg = findBasisInDegree(deg, dom, previousPolys, basisHash));

        
        -- This is the part that should be in parallel
        oldPolys := previousPolys;
        for deg in lats do (
            previousPolys = previousPolys | componentOfIdeal(deg, dom, phi, oldPolys, basisHash);
        );
    );
    return previousPolys
);


parComponentOfIdeal = (deg, dom, phi, basisHash, polyHash) -> () -> (

    previousPolys := flatten for z in gens(dom) list if polyHash#?(deg - degree(z)) then polyHash#(deg - degree(z));

    previousPolys = delete(null, previousPolys);

    monomialBasis := findBasisInDegree(deg, dom, previousPolys, basisHash);

    (mons, coeffs) := coefficients(phi(sub(monomialBasis, source phi)));

    -- find the linear relations among coefficients
    K := gens ker sub(coeffs, QQ);

    newGens := flatten entries (monomialBasis * K);

    return (monomialBasis, newGens)
    )

parComponentsOfIdeal = (phi, d) -> (

    n := numgens(source phi);
    omega := maxGrading(phi);
    dom := newRing(source phi, Degrees => omega);

    basisHash := new MutableHashTable;
    polyHash := new MutableHashTable;

    threadVariable monomialBasis;
    threadVariable mons;
    threadVariable coeffs;
    threadVariable newGens;
    threadVariable previousPolys;
    threadVariable badMonomials;
    threadVariable mat;
    threadVariable K;

    for i in 1..d do(

        allMonomials := sub(basis(i, source phi), dom);
        uniqueDegrees := unique apply(flatten entries allMonomials, m -> degree m);

        tasksPerComponent = new MutableHashTable from apply(uniqueDegrees, deg -> deg => createTask parComponentOfIdeal(deg, dom, phi, basisHash, polyHash));
        tasks = values(tasksPerComponent) / schedule;
        setIOSynchronized();
        print("all tasks scheduled");
        while (not all(tasks / isReady, i -> i)) do(

            sleep 1;
            print(tasks / isReady);
            );

        scan(uniqueDegrees, deg -> (

            (monomialBasis, newGens) := taskResult(tasksPerComponent#deg);

            basisHash#deg = flatten entries monomialBasis;
            polyHash#deg = newGens;
            )
            );
        );

    return polyHash;
    )






