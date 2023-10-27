restart
needsPackage "MultigradedImplicitization"

--Ring definitions: 
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,4)];
S=sort elements (set {1,2,3,4})^**4/splice/splice;
var=toList apply(S,i->(symbol x)_i);
Rx=K[var];

--Retrieve tensor on our desired basis
use R;
pbar=value get "4leaves_tensor_FinalBasis.txt";

--Different types of entries in the tensor in the new basis
nonMonomial=select(S,i->(length terms pbar_(position(S,j->j==i),0)>1));

nonZeroEntries=S_(positions(flatten entries pbar,i->i!=0));

monomialNonZeroEntries=select(nonZeroEntries,i->(length terms pbar_(position(S,j->j==i),0)==1));


-- substitute in a random root distribution
pbar = sub(pbar, {p_1 => 1/9, p_2 => 2/63, p_3 => 11/21, p_4 => 1/3});
nonZeroPBar = delete(null, apply(flatten entries pbar, m -> if m != 0 then m));

--Parametrization: almost monomial map

S = QQ[l_(1,1)..l_(5,4)];
T = QQ[apply(nonZeroEntries, ind -> x_ind)];
f = map(S,T, apply(nonZeroPBar, i -> sub(i,S)));

--grading
D = maxGrading f

--dimension of grading
rank(D)

--compute the quadratics
G = componentsOfKernel(2,f);

H = (flatten delete({},values(G))) / (g -> sub(g, source f))

--double check in the ideal
(unique (H / f)) == {0}


--stick in a file
fileName = "TN93_quartet_quadrics" << "";

for g in H do (
    fileName << g << endl
);

fileName << close;
