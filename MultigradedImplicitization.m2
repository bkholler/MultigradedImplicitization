needsPackage "gfanInterface";

-- Input : f : dom = kk[x_1..x_n] --> codom = kk[y_1..y_m]
-- Output : lineality space of Groebner fan
maxGrading = f -> (
    -- set up elimination ideal
    dom := source f;
    codom := target f;
    elimRing := dom ** codom;
    X := vars dom;
    elimIdeal := ideal(sub(X, elimRing) - sub(f(X), elimRing));
    return transpose linealitySpace(gfanHomogeneitySpace(elimIdeal))
);

findBasisInDegree = (G, R, deg) -> (  
    -- if G = {}, then you just have to use all of R_deg
    if #G == 0 then (
        return basis(deg, R);
    );

    -- otherwise, we shift G in all possible ways to land in R_deg
    L := apply(G, g -> g*basis(deg - degree(g), R));
    
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

-- Input: 
--      deg = homogeneous degree of kernel you want to comput
--      G = known generators of degree less than deg
--      f : R --> S
--      D = Degrees
-- Output: 
--      newGens = generators in R_deg which are not linear combinations of things in G
componentOfIdeal = (deg, G, f, D) -> (
    -- number variables of R
    n := numColumns(vars source f);

    -- apply degrees to R if not already done
    dom := newRing(source f, Degrees => D_(toList(0..n-1)));

    -- substitute G into graded ring
    domG := apply(G, g -> sub(g,dom));

    -- The span of monomialBasis is all you need to search for, since G is known
    monomialBasis := findBasisInDegree(domG, dom, deg);

    -- collect coefficients into a matrix
    (mons, coeffs) := coefficients(f(sub(monomialBasis, source f)));

    -- find the linear relations among coefficients
    K := gens ker sub(coeffs,QQ);

    -- substitute back into polynomials and return
    newGens := flatten entries (monomialBasis * K);

    return newGens
);


