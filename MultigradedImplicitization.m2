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

-- does linear algebra, no optimization
-- deg is the degree as a list
-- f : dom --> codom
componentOfIdeal = {MaxGrading => null} >> opts -> (deg, f) -> (
    -- maximal grading of ideal
    D := if opts.MaxGrading === null then maxGrading(f, dom, codom) else opts.MaxGrading; 

    dom := source f;
    codom := target f;

    -- introduce gradings to domain
    n := numColumns(vars dom);
    degreesDom := D_(toList(0..n-1));
    grDom := newRing(dom, Degrees => degreesDom);

    -- Basis of R_deg
    B := basis(deg, grDom);

    -- image of R_deg under f
    fB := f(sub(B,dom));
    
    --gather coefficients
    (mons, coeffs) := coefficients(fB);

    -- solve coeffs = 0
    K := gens ker sub(coeffs,QQ);

    -- return basis of R_deg \cap ker(f)
    return sub(B,dom) * K
);

findBasisInDegree = (G, R, deg) -> (  
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


-- G = known generators of degree less than deg
restrictedComponentOfIdeal = (deg, G, f, D) -> (
    --domain and codomain of f
    n := numColumns(vars source f);
    dom := newRing(source f, Degrees => D_(toList(0..n-1)));
    codom := target f;

    --G in dom
    domG := apply(G, g -> sub(g,dom));

    -- basis for V where  dom_deg = span(G) \oplus V
    monomialBasis := findBasisInDegree(domG, dom, deg);

    (mons, coeffs) := coefficients(f(sub(monomialBasis, source f)));

    K := gens ker sub(coeffs,QQ);

    newGens := flatten entries (monomialBasis * K);

    return newGens
)