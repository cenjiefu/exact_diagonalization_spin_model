#include "ed_lanczos/common.h"

// Brent's method from https://rosettacode.org/wiki/Roots_of_a_function#C.2B.2B
double brents_fun(std::function<double (double)> f, double lower, double upper, double tol, unsigned int max_iter)
{
  double a = lower;
  double b = upper;
  double fa = f(a); // calculated now to save function calls
  double fb = f(b); // calculated now to save function calls
  double fs = 0;    // initialize 
 
  if (!(fa * fb < 0))
  {
    std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
    return -11;
  }
 
  if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
 
  double c = a;     // c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;     // precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;    // boolean flag used to evaluate if statement later on
  double s = 0;     // Our Root that will be returned
  double d = 0;     // Only used if mflag is unset (mflag == false)
 
  for (unsigned int iter = 1; iter < max_iter; ++iter)
  {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tol)
    {
      std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
      return s;
    } // end if
 
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
        + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
        + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
 
      // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (  ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
        ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
        ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
        ( mflag && (std::abs(b-c) < tol) ) ||
        ( !mflag && (std::abs(c-d) < tol))  )
    {
      // bisection method
      s = (a+b)*0.5;
 
      mflag = true;
    }
    else
    {
      mflag = false;
    }
 
    fs = f(s);  // calculate fs
    d = c;    // first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;    // set c equal to upper bound
    fc = fb;  // set f(c) = f(b)
 
    if ( fa * fs < 0) // fa and fs have opposite signs
    {
      b = s;
      fb = fs;  // set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;  // set f(a) = f(s)
    }
 
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);   // swap a and b
      std::swap(fa,fb); // make sure f(a) and f(b) are correct after swap
    }
 
  } // end for
 
  std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
  return -1;
 
}

// Binary search of the k_table which is sorted by construction
int64_t findIndexInKTable(uint32_t representative_state, vector<uint32_t> & k_table){
    uint32_t index_min = 0;
    uint32_t index_max = k_table.size()-1;
    while(true){
        uint32_t index = index_min + (index_max - index_min)/2;
        if ( representative_state < k_table[index] ) {
            if ( index == 0 )
                return -1;
            index_max = index - 1;
        } else if ( representative_state > k_table[index] ){
            index_min = index + 1;
        } else {
            return index;
        }
        if ( index_min > index_max ){
            return -1;
        }
    }
};


uint32_t qua2dec(vector<unsigned char>& qua){
    uint32_t dec = 0;
    // vector<unsigned char> temp = qua;
    // std::reverse_copy(qua.begin(), qua.end(),temp.begin());
    std::reverse(qua.begin(), qua.end());
    for (size_t j=0; j<qua.size(); j++){
        dec = dec + uint32_t(qua[j])*pow(4,j);
    }
    std::reverse(qua.begin(), qua.end());
    return dec;
}

vector<unsigned char> dec2qua(uint32_t decimal, uint32_t N_digit){
    vector<unsigned char> qua(N_digit, 0);
    size_t j = N_digit-1;
    while(decimal>0){
        qua[j] = decimal%4;
        decimal = decimal/4;
        j = j-1;
    }
    return qua;
}

void printState(vector<unsigned char>& qua){
  for (size_t k=0; k< qua.size();k++){
        cout << double(qua[k]) << " ";
    }
  cout << endl;
}

void printState(uint32_t decimal, uint32_t N_digit){
  vector<unsigned char> binary(N_digit, 0);
  size_t j = N_digit-1;
  while(decimal>0){
      binary[j] = decimal%2;
      decimal = decimal/2;
      j = j-1;
  }

  for (size_t k=0; k< binary.size();k++){
        cout << double(binary[k]) << " ";
    }
  cout << endl;
}

double maxAbsRowSum(Eigen::Ref<SpMatCD> A){
  double maxsum = 0;
  for (int k=0; k < A.outerSize(); ++k){
    double rowsum = 0;
    for (Eigen::Ref<SpMatCD>::InnerIterator it(A,k); it; ++it){
      rowsum += std::abs(it.value());
    }
    if (rowsum >= maxsum){
      maxsum = rowsum;
    }
  }
  return maxsum;
}