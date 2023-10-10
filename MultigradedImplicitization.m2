needsPackage "gfanInterface";

-- f : dom = kk[x_1..x_n] --> codom = kk[y_1..y_m]
maxGrading = (f,dom,codom) -> (
    -- set up elimination ideal
    elimRing := dom ** codom;
    elimIdeal := ideal(apply(flatten entries vars dom, x -> sub(x,elimRing) - sub(f(x),elimRing)));
    -- return lineality space of the Groebner fan
    return linealitySpace(gfanHomogeneitySpace(elimIdeal))
)



--example
R = QQ[x,y,z]
S = QQ[s,t]
f = map(S,R,{t^8*(1-t)^(10)*s^(15), t^5*(1-t)^7*s^(10), t^3*(1-t)^4*s^6})
maxGrading(f,R,S)
