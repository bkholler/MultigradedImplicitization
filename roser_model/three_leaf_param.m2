load "~/Documents/my_scripts/MultigradedImplicitization/MultigradedImplicitization.m2";

R = QQ[q_(0,0,0)..q_(3,3,3)]
S = QQ[b_0..d_2,p_1..p_4]

makeMatrix = (k,S) -> (
    m = new MutableMatrix from map(S^4,S^4,0);
    for i from 0 to 2 do (
       
        for j from i+1 to 3 do (
            
            if i == 0 and j == 1 then (
                m_(i,j) = p_2*c_k;
                m_(j,i) = p_1*c_k;
            );

            if i <= 1 and j >= 2 then (
                m_(i,j) = p_(j+1)*b_k;
                m_(j,i) = p_(i+1)*b_k;
            );

            if i == 2 and j == 3 then (
                m_(i,j) = p_4*d_k;
                m_(j,i) = p_3*d_k;
            );

        );
    );

    ones = sub(transpose(matrix {{1,1,1,1}}),S);
    rowSums = (matrix m) * ones;
    for i from 0 to 3 do (
        m_(i,i) = 1 - rowSums_(i,0);
    );

    return matrix m
)

M = apply(3, i -> makeMatrix(i,S))

states = apply(flatten entries vars R, x -> last baseName x);

im = apply(flatten entries vars R, 
    x -> (
        leaf_state = last baseName x;

        -- i = internal state of root
        sum for i from 0 to 3 list (
            p_(i+1) * product(toList apply(0..2, j -> (M#j)_(i,leaf_state#j)))
        )
    )
)

phi = map(S,R,im)

maxGrading phi




