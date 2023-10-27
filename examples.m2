
-- Example 1
-- this example is broken at the moment. i think it's because of the non-standard Z-grading
restart
needsPackage "MultigradedImplicitazation"
R = QQ[x,y,z];
S = QQ[s,t];

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im);









--Example 2
--still fine
restart
needsPackage "MultigradedImplicitization"
load "~/Documents/my_scripts/sunlets/sunletQuadGens.m2"

n = 5;
R = qRing n;
f = sunletParam n;
S = ring(f#0);
f = map(S,R,f);

G = componentsOfKernel(3,f)

G = new HashTable from G;

G = apply(flatten values(G), g -> sub(g, R))








--Example 3
restart 
needsPackage "MultigradedImplicitization"
needsPackage "Polyhedra"

R = QQ[x..z]
S = QQ[u,v]

f1 = map(S,R,{u^2,u*v,v^2})
D1 = maxGrading f1
I1 = ker f1

f2 = map(S,R,{(u+v)^2,(u+v)*(u-v),(u-v)^2})
D2 = maxGrading f2
I2 = ker f2

I1 == I2
D1 != D2

isSubset(I1, toricMarkov(D2,R))



--Example 4
--Grassmannian
restart
needsPackage "MultigradedImplicitization"


grass = (n,k) -> (
    ind := subsets(0..n-1,k) / (i -> toSequence i);
    R := QQ[apply(ind, i -> p_i)];
    S := QQ[x_(0,0)..x_(k-1,n-1)];
    
    A := new MutableMatrix from map(S^k, S^n ,0);
    scan(0..k-1, i -> scan(0..n-1, j -> A_(i,j) = x_(i,j)));
    A = matrix A;
    phiList := toList apply(ind, i -> det (A_(toList i)));
    return map(S,R,phiList)
)

grassWithMI = (n,k) -> (
    phi := grass(n,k);
    return componentsOfKernel(2,phi)
)

grassWithMICoarseGrading = (n,k) -> (
    phi := grass(n,k);
    D := matrix{toList(binomial(n,k):1)};
    return componentsOfKernel(2, phi, Grading => D)
)

grassWithGB = (n,k) -> (
    return ker grass(n,k)
)

timing scan(4..7, n -> scan(1..n//2, k -> grassWithMI(n,k))) -- uses componentsOfKernel with maximal grading
timing scan(4..7, n -> scan(1..n//2, k -> grassWithMICoarseGrading(n,k))) -- uses componentsOfKernel with coarse grading
timing scan(4..7, n -> scan(1..n//2, k -> grassWithGB(n,k))) -- just gb



G = grassWithMI(6,3)
G = new HashTable from G;
I36 = ideal flatten values(G)
dim I36 == 3*(6-3) + 1
isPrime I36

