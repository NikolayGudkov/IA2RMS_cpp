# Overview

This is the C++ implementation of the Independent Doubly Adaptive Rejection Metropolis Sampling (IA2RMS) withing Gibbs Sampling
algorithm. The algorithm overcomes the drawbacks of the [ARS and ARMS algorithms](http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html) and allows to perform efficient sampling from a broad spectrum of univariate distributions.
It can be used as a stand-alone algorithm as well as a within-Gibbs sampler.

The original algorithm is described in [L. Martino, J. Read, D. Luengo, "Independent Doubly Adaptive Rejection Metropolis Sampling within Gibbs Sampling",
IEEE Transactions on Signal Processing, Volume 63, Issue 12, Pages: 3123-3138, 2015.](http://www.lucamartino.altervista.org/TSP_IA2RMS.pdf)

We implement three types of proposal construction with:
*  uniform pieces (type 'u'),
*  exponential pieces (type 'e'),
*  linear pieces (type 'l').

The code was tested on two examples from [L. Martino, J. Read, D. Luengo (2015)](http://www.lucamartino.altervista.org/TSP_IA2RMS.pdf):
1. Multimodal Target: Mixture of Gaussians,
1. Heavy-Tailed Distribution: Levy distribution.

# Content

* __unif.cpp__ - defines simulation of uniform random variables
* __Interval.hpp__ - defines class Interval. We assign several private and public members to it. In particular, we are interested in such properties of an interval as its boundaries, values of the target and proposal distribution at the points located in the interval, area/weight associated with the interval. Objects of this class are used as building blocks to construct the proposal distribution which approximates the target distribution.
* __IA2RMS.cpp__ and __IA2RMS.hpp__ - realisation of the IA2RMS algorithm.

# Usage 

## Description
Function **IA2RMS** simulates *N* points from the target continuous distribution provided by *target_function*.
### Input:
* __*target_function*__ - continuous probability distribution function known up to a constant of proportionality. *Note that the function has to be defined in the log-scale.*
* __*S0*__ - initial set of support points. 
*Note that to avoid numerical issues, specify sufficient number (e.g. 2-5) of initial points, S0, in the area where the target distribution has the largest mass of probability.*
* __*N*__ - number of random samples required.
* __*x*__ - vector to store simulated random samples.
* __*s_min*__ and __*s_max*__ - bounds of the domain where the target function is defined (e.g. s_min=-inf, s_max=+inf).
*Note that, if s_min (s_max) is -inf (inf), by default the left (right) tail is defined by the exponential construction (i.e. type 'e' is selected). 
Also, if the target function is not well defined in either s_min or s_max, numerical issues may occur. In this case, use finite values instead (e.g. use 10^(-9) instead of 0).*
* __*t*__ - character which specifies the type of pieces of the proposal distribution (exponential 'e', uniform 'u', or linear 'l').

### Output: 
* area under the proposal distribution curve.

There is also an option to create the target probability density curve and the proposal distribution curve constructed via the IA2RMS algorithm. These values are saved in a file, which can be used by the gnuplot.

## Example of usage

### Simulation from the Gaussian distribution 
```
#include <iostream>
#include "IA2RMS.hpp"

using namespace std;

//##############################/
// Target distribution function /
//##############################/
double target_function(const double& x){
    // This is the non-normalised distribution function in the log-scale we want to sample from.
    // We asume that this function is given in the log scale. That is, if f(x) is the target density function defined in the normal scale, use log(f(x)) as a target function.
    // For example, to sample from the normal distribution with density function proportional to exp(-(x-mu)^2/sigma^2),
    // where mu is mean and sigma^2 is variance, use the target -(x-mu)^2/(sigma*sigma).
    
    return -(x-1.)*(x-1.)/(2.*2.);
}


// Define the boundaries of the target function domain
 const double s_min=-inf, s_max=inf;

int main(){
    srand(1);
// Define the required number of simulations
    const long N1=10, N2=10000;
 
// Define a vector to store simulated values
    vector<double> x(N1);
 
// Define a vector of initial support points
    vector<double> S0={-10,-5,0,5,10};
 
// Call the IA2RMS function to simulate random values
    IA2RMS(target_function,S0,N1,x,s_min,s_max,'l');
 
// Print out the simulated values
    for (auto const& value: x) {cout<<value<<endl;}

// Extend the size of x
    x.resize(N2);
     
// Call the IA2RMS function to simulated more random values
    IA2RMS(target_function,S0,N2,x,s_min,s_max,'l');

// Print out the average of the simulated values
    cout<<"Mean of x is "<<accumulate(x.begin(),x.end(),0.)/N2<<endl;
 
return 0;
}
```
### Output:
```
Simulated values:
-8.18659
1.35211
2.21899
0.0892406
-1.90427
2.05544
0.890752
1.37205
2.17782
1.36979

Mean of x is 0.99525
Program ended with exit code: 0
```

### Examples from [L. Martino, J. Read, D. Luengo (2015)](http://www.lucamartino.altervista.org/TSP_IA2RMS.pdf)
```
#include <iostream>
#include "IA2RMS.hpp"

using namespace std;

int main() {

   // Multimodal Target:
   // Mixture of Gaussians
     double a, b, tmp;
     const long N_sim=2000;
     const int N=5000;
     vector<double> x(N);
     double sum=0., S0_min, S0_max;
     
     S0_min=-10; S0_max=10;
    
     // Specify the mixture of Gaussians as a target distribution
     auto target_1 = [](const double& x) {return log(0.3*exp(-(x+5)*(x+5)/2)+0.3*exp(-(x-1)*(x-1)/2)+0.4*exp(-(x-7)*(x-7)/2));};
    
     for (long i=0; i<N_sim; i++){
         a=(S0_max-S0_min)*randu()+S0_min; b=(S0_max-S0_min)*randu()+S0_min;
         if (a>b) {tmp=b; b=a; a=tmp;}
         
         // Specify the initial support points
         vector<double> S0={S0_min,a,b,S0_max};
         
         // Call the simulation algorithms
         IA2RMS(target_1,S0,N,x,-inf,inf,'l');

         // Compute the sum of simulated values
         sum+=accumulate(x.begin(),x.end(),0.)/N;
       }
    cout<<"Mean of Mixture of Gaussians: "<<endl;
    cout<<"Mean="<<sum/double(N_sim)<<endl;

    cout<<endl;
    
    // Heavy-Tailed Distribution:
    // Le'vy distribution
    vector<double> C(N_sim), MSE(N_sim), mu;
    
    // Specify the Le'vy distribution as a target
    auto target_2 = [](const double& x) {return log(exp(-1/x)/pow(x,1.5));};
    
    // Define boundaries of the target function domain
    S0_min=pow(10,-5); S0_max=10.;
    
    for (long i=0; i<N_sim; i++){
        a=(S0_max-S0_min)*randu()+S0_min; b=(S0_max-S0_min)*randu()+S0_min;
        if (a>b) {tmp=b; b=a; a=tmp;}
        
        // Specify the initial support points
        vector<double> S0={S0_min,a,b};
        
        // Obtained the total area under the proposal distribution curve
        C[i]=IA2RMS(target_2,S0, N, x, pow(10,-7), inf,'l');
        
        // Compute the Mean Squared Error (MSE)
        MSE[i]=pow(C[i]-1./sqrt(M_PI),2.);
    }
    
    cout<<"Estimation of the constant in the Le'vy distribution:"<<endl;
    cout<<"1/c_p="<<accumulate(C.begin(),C.end(),0.)/N_sim<<endl;
    cout<<"MSE="<<sqrt(accumulate(MSE.begin(),MSE.end(),0.)/N_sim)<<endl;
    
    cout<<endl;

    return 0;
}
```

## Output
```
Mean of Mixture of Gaussians: 
Mean=1.60188

Estimation of the constant in the Le'vy distribution:
1/c_p=0.563965
MSE=0.0012882

Program ended with exit code: 0
```
