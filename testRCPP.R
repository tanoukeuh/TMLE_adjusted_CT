## inline to compile, load and link the C++ code
require(inline)

## we need a pure C/C++ function as the generated function
## will have a random identifier at the C++ level preventing
## us from direct recursive calls
incltxt <- '
double powerCabral(int base, int exponent){
    return pow(base,exponent);
}

int fibonacci(const int x) {
if (x == 0) return(0);
if (x == 1) return(1);
return fibonacci(x - 1) + fibonacci(x - 2) + powerCabral(x,2);
}'

## now use the snippet above as well as one argument conversion
## in as well as out to provide Fibonacci numbers via C++
fibRcpp <- cxxfunction(signature(xs="int"),
                       plugin="Rcpp",
                       incl=incltxt,
                       body
                       ='
                       int x = Rcpp::as<int>(xs);
                       return Rcpp::wrap( fibonacci(x) );
                       ')
