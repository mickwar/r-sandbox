f = function(x, tol=1e-12){
    curr = exp(x)
    test = exp(-curr)
    while (abs(curr - test) > tol){
        curr = test
        test = f(-curr)
        }
    return (test)
    }

x = f(0)
