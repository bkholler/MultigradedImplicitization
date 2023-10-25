needsPackage "gfanInterface";

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

findBasisInDegree = (G, R, deg, basisHash) -> (

    if #G == 0 then (
        return basis(deg, R);
    );

    -- otherwise, we shift G in all possible ways to land in R_deg

    L := apply(G, g -> (
            checkDegree := deg - degree(g);
            if basisHash#?checkDegree then (
                g*basisHash#checkDegree
            ) else (
                g*basis(checkDegree, R)
            )
        )
    );
    
    -- stick em all in a matrix
    mat := L#0;
    scan(1..#L-1, i -> mat = mat | L#i);

    -- and collect coefficients.
    (mons, coeffs) := coefficients(mat);

    -- find the independent linear relations
    coeffs = mingens(image sub(coeffs, QQ));

    -- remove monomials corresponding to pivots
    badMonomials := apply(pivots coeffs, i -> mons_(0,i#0));
    monomialBasis := flatten entries basis(deg, R);
    scan(badMonomials, m -> monomialBasis = delete(m, monomialBasis));

    return matrix{monomialBasis}
);

-- G = known generators of degree less than deg
componentOfIdeal = (deg, G, phi, dom, basisHash) -> (
    
    --G in dom
    domG := apply(G, g -> sub(g,dom));

    -- The span of monomialBasis is all you need to search for, since G is known
    monomialBasis := findBasisInDegree(domG, dom, deg, basisHash);
    --monomialBasis = basisHash#deg;

    -- collect coefficients into a matrix
    (mons, coeffs) := coefficients(phi(sub(monomialBasis, source f)));

    -- find the linear relations among coefficients
    K := gens ker sub(coeffs,QQ);

    newGens := flatten entries (monomialBasis * K);

    return newGens
);

-- this guy only works if it has usual Z-grading
-- things to add to this function:
--      make it work when it doesn't have usual Z-grading, and check for this
--      make an option to override this maybe if you already have D that you want to use
--      make an option if the ring already has the appropriate grading
--      write it in Oscar
componentsOfIdeal = (phi, d) -> (

    n := numgens(source phi);
    omega := maxGrading(phi);
    dom := newRing(source phi, Degrees => omega);

    G := {};
    basisHash := new MutableHashTable;
    
    -- assumes homogeneous with normal Z-grading
    for i in 0..d do (
        B := sub(basis(i, source phi), dom);
        lats := unique apply(flatten entries B, m -> degree m);

        scan(lats, i -> basisHash#i = findBasisInDegree(G,dom,i,basisHash));

        
        -- This is the part that should be in parallel
        oldG = G;
        for deg in lats do (
            G = G | componentOfIdeal(deg, oldG, phi, dom, basisHash);
        );
    );
    return G
    -- return G / (g -> sub(g, source phi))
);