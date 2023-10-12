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

    if #G == 0 then (
        return basis(deg, R);
    );

    L := apply(G, g -> g*basis(deg - degree(g), R));
    mat := L#0;
    scan(1..#L-1, i -> mat = mat | L#i);
    (mons, coeffs) := coefficients(mat);
    coeffs = mingens(image sub(coeffs, QQ));
    badMonomials := apply(pivots coeffs, i -> mons_(i#0));
    monomialBasis := flatten entries basis(deg, R);
    scan(badMonomials, m -> monomialBasis = delete((flatten entries m)#0, monomialBasis));

    return matrix{monomialBasis}
);


-- G = non-empty known generators of degree less than deg
componentOfIdeal = (deg, G, f) -> (
    
    --G in dom
    domG := apply(G, g -> sub(g,dom));

    -- basis for V where  dom_deg = span(G) \oplus V
    monomialBasis := findBasisInDegree(domG, dom, deg);

    (mons, coeffs) := coefficients(f(sub(monomialBasis, source f)));

    K := gens ker sub(coeffs,QQ);

    print(K);

    newGens := flatten entries (monomialBasis * K);

    return newGens
)


componentsOfIdeal = (degs, f, D) -> (

    n := numgens(source phi);
    dom := newRing(source f, Degrees => B_(toList(0..n-1)));
    codom := target f;

    return for deg in degs list componentOfIdeal(flatten entries deg, {}, phi)
)

