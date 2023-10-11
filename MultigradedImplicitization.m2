needsPackage "gfanInterface";

-- Input : f : dom = kk[x_1..x_n] --> codom = kk[y_1..y_m]
-- Output : lineality space of Groebner fan
maxGrading = (f, dom, codom) -> (
    -- set up elimination ideal
    elimRing := dom ** codom;
    X := vars dom;
    elimIdeal := ideal(sub(X, elimRing) - sub(f(X), elimRing));
    return transpose linealitySpace(gfanHomogeneitySpace(elimIdeal))
);

-- does linear algebra, no optimization
-- deg is the degree as a list
-- f : dom --> codom
componentOfIdeal = (deg, f, dom, codom) -> (
    -- maximal grading of ideal
    D := maxGrading(f, dom, codom);

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