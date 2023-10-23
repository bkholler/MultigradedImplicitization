restart
needsPackage "FourTiTwo"
load "MultigradedImplicitization.m2"


allowableThreads = maxAllowableThreads

temp = i -> () -> (

	B := transpose matrix {{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1}};
	A := B**B**matrix(toList(i:{1,1}));
	S := QQ[x_1..x_(numcols A)];
	R := QQ[t_1..t_(numrows A)];
	phi := map(R, S, apply(numcols A, i -> R_(flatten entries A_i)));
	
	return ker phi;
	)

test = i -> (

	B := transpose matrix {{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1}};
	A := B**B**matrix(toList(i:{1,1}));
	S := QQ[x_1..x_(numcols A)];
	R := QQ[t_1..t_(numrows A)];
	phi := map(R, S, apply(numcols A, i -> R_(flatten entries A_i)));
	
	return ker phi;
	)

n = 15
tasks = apply(2..n, i -> createTask(temp(i)));
schedTasks = tasks / schedule;
count = 0;
while (not all(schedTasks / isReady, i -> i)) do(

	count = count + 1;
	sleep 1
	);
parallelOut = schedTasks / taskResult;
count


out = time apply(2..n, i -> test(i))


certs = for i from 0 to n-2 list(

	I = parallelOut_i;
	J = out_i;

	I == sub(J, ring I)
	)

all(certs, i -> i)