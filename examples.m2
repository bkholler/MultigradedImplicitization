
-- Example 1
restart
load "MultigradedImplicitization.m2"
R = QQ[x,y,z];
S = QQ[s,t];

im = {t^8*(1-t)^(10)*s^(15), t^5*(1-t)^(7)*s^(10), t^3*(1-t)^4*s^6};
f = map(S,R,im);

D = maxGrading(f);

G = {};
L = for i in 0..60 list (
    G = G | componentOfIdeal({i}, G, f, D;
    {i, G}
);
netList L


--Example 2
restart
allowableThreads = maxAllowableThreads;
load "MultigradedImplicitization.m2"
load "/Users/bkholler/My Drive/GroupBasedPhylogeneticNetworks/CFN_Network_Ideals/code/suppMaterials/sunletQuadGens.m2"
--setIOSynchronized()

n = 6;
R = qRing n;
phi = sunletParam n;
S = ring(phi#0);
phi = map(S,R,phi);
parComponentsOfIdeal(phi, 2)

D = matrix{toList(2^(n-1):1)};

--Example 3
restart 
load "MultigradedImplicitization.m2"
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

-- Example 4
restart
load "MultigradedImplicitization.m2"

B = transpose matrix {{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1}};
A = (B | 0*B) || (0*B | B)
R = QQ[x_(1,1)..x_(2,2), y_(1,1)..y_(2,2)]
S = QQ[t_1..t_(numrows(A))]
phi = map(S, R, apply(numcols(A), i -> S_(flatten entries A_i)))
dom = newRing(R, Degrees => A)
I = ker phi;
G = I_* / (i -> sub(i, dom))


n = numgens(source phi);
polyHash = new MutableHashTable from apply(gens(dom), i -> degree(i) => {})
basisHash = new MutableHashTable from apply(gens(dom), i -> degree(i) => matrix{{i}})
scan(unique apply(flatten entries(sub(basis(2, source phi), dom)), m -> degree m), i -> basisHash#i = basis(i, dom))
scan(unique apply(flatten entries(sub(basis(3, source phi), dom)), m -> degree m), i -> basisHash#i = basis(i, dom))
allowableThreads = 4;

for d from 2 to 3 do(

    allMonomials = sub(basis(d, source phi), dom);
    degs = unique apply(flatten entries allMonomials, m -> degree m);
    k = floor(#degs/4);
    degsPerThread = apply(4, i -> if i!= 4 -1 then degs_(toList(k*i..k*(i+1)-1)) else degs_(toList(k*i..#degs-1)));
    tasks = new MutableHashTable;

    for i from 0 to (#degsPerThread - 1) do(

        for j from 0 to (#degsPerThread_i - 1) do(

            curDeg = (degsPerThread_i)_j;
            if j > 0 then prevDeg = (degsPerThread_i)_(j-1);

            tasks#curDeg = createTask(componentOfIdeal(curDeg, dom, phi, {}, basisHash));

            if j == 0 then continue else addStartTask(tasks#prevDeg, tasks#curDeg);
            );
        );

    degsPerThread / (i -> schedule tasks#(i_0));

    print("all tasks scheduled");

    while (not all(values(tasks) / isReady, i -> i)) do(

        sleep 1;
        print(values(tasks) / isReady);
        );

    print("all tasks done");


    for i in degs do(

    out = taskResult(tasks#i);
    print({i, out});
    polyHash#i = out;
    );

    print(concatenate("degree done:", toString(d)));
    print("\n");
    )

peek basisHash


outputTask = createTask retrieveOut(values(tasksPerComponent));

addStartTask(values(tasksPerComponent), outputTask);

