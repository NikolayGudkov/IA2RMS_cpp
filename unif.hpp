//
//  unif.hpp
//
//  Created by Nikolay Gudkov on 27/4/20.
//  Copyright © 2020 Nikolay Gudkov. All rights reserved.
//
// Simulate a value from the uniform diftribution in [0,1].

#ifndef unif_h
#define unif_h
#include <random>

inline double randu(){
    return rand()/double(RAND_MAX);
};

#endif /* unif_h */
