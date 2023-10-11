needsPackage "gfanInterface";

-- f : dom = kk[x_1..x_n] --> codom = kk[y_1..y_m]
maxGrading = (f, dom, codom) -> (
    -- set up elimination ideal
    elimRing := dom ** codom;
    elimIdeal := ideal(apply(flatten entries vars dom, x -> sub(x,elimRing) - sub(f(x),elimRing)));
    -- return lineality space of the Groebner fan
    return transpose linealitySpace(gfanHomogeneitySpace(elimIdeal))
)

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
)












-- ezpz example
R = QQ[x,y,z]
S = QQ[s,t]
f = map(S,R,{t^8*(1-t)^(10)*s^(15), t^5*(1-t)^7*s^(10), t^3*(1-t)^4*s^6});


timing componentOfIdeal({30}, f, R, S)

timing ker f



