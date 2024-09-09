await = method(Dispatch => Thing)
await Thing := identity
await BasicList := L -> apply(L, await)
await HashTable := H -> applyValues(H, await)
-- TODO: applyValues should take a mutable hash table
await MutableHashTable := H -> ( scan(keys H, k -> H#k = await H#k); H )
await Task := await @@ taskResult

async = method(Dispatch => Thing)
async Function := f -> (x -> schedule(f, x))

end--

repeat = (n, R) -> (
    if n == 1 then return n;
    r := {repeat(n//2, R), repeat(n//2, R)};
    basis(10, ideal vars R);
    r)
repeat = async repeat

end--
restart
needs "parallel.m2"

R = ZZ/101[x_0..x_8]

f = x -> (sleep 2; 2*x)
f 2

repeat(4, R)

g = async f
r = g 2;

await r

elapsedTime repeat(8, R)



elapsedTime await repeat(allowableThreads*5, R)
elapsedTime await repeat(allowableThreads*5, R)

allowableThreads = maxAllowableThreads

taskResult o2


R = ZZ/101[x_0..x_10]
repeat = n -> (
    if n == 1 then return n;
    r := {repeat(n//2), repeat(n//2)};
    basis(10, ideal vars R);
    r)
repeat 5
