
-- Example 1
restart
needsPackage "MultigradedImplicitization"
R = QQ[x,y,z,w];
S = QQ[u,v,t];

phiList = {u^(8)*(1-u)^(10)*v^(15)*t^(16), u^(5)*(1-u)^(7)*v^(10)*t^(11), u^(3)*(1-u)^(4)*v^(6)*t^(7), t};
phi = map(S,R,phiList);

D = maxGrading phi

G = componentsOfKernel(11,phi)

G = new HashTable from G;

G = apply(flatten values(G), g -> sub(g, R))

G / factor




--Example 2
restart
needsPackage "MultigradedImplicitization"
load "~/Documents/my_scripts/sunlets/sunletQuadGens.m2"

n = 5;
R = qRing n;
f = sunletParam n;
S = ring(f#0);
f = map(S,R,f);

G = componentsOfKernel(2,f)

G = new HashTable from G;

G = apply(flatten values(G), g -> sub(g, R))








--Example 3
restart 
needsPackage "MultigradedImplicitization"
needsPackage "Polyhedra"

R = QQ[x..z]
S = QQ[u,v]

phi = map(S,R,{u^2,u*v,v^2})
A = maxGrading phi
I = ker phi

G = componentsOfKernel(10,phi)
G = new HashTable from G;
G = apply(flatten values G, g -> sub(g,R))




psi = map(S,R,{(u+v)^2,(u+v)*(u-v),(u-v)^2})
B = maxGrading psi
J = ker psi 

H = componentsOfKernel(10,psi)
H = new HashTable from H;
H = apply(flatten values H, g -> sub(g,R))


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



phi = grass(6,3);
R = source phi

G = grassWithMI(6,3)
G = new HashTable from G;
G = apply(flatten values G, g -> sub(g,R))
I36 = ideal G
dim I36 == 3*(6-3) + 1
isPrime I36

G / phi

