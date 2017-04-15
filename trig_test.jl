


function test_small(r)

    println("small number test")
    @benchmark 10*r ;
    @benchmark sin(r);
    @benchmark sin(mod(r, 2*pi));
end

function test_big(r)
    println("Big number test")
    r = length(r)*r
    @benchmark 10*r ;
    @benchmark sin(r);
    @benchmark sin(mod(r, 2*pi));
end

@fastmath function test_big_fast(r)
    println("Big number test")
    r = length(r)*r
    @benchmark 10*length(R)r ;
    @benchmark sin(r);
    @benchmark sin(mod(r, 2*pi));
end

N = 10^6
const r = 2*pi*rand(N) ;

test_small(r) ;

test_big(r) ;

test_big_fast(r) ;
