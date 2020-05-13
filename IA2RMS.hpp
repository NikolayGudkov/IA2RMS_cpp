//
//  IA2RMS.hpp
//
//  Created by Nikolay Gudkov on 19/4/20.
//

#ifndef IA2RMS_hpp
#define IA2RMS_hpp
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
#include <fstream>
#include "Interval.hpp"

using namespace std;

// Select an interval from a vector I using their areas/weights and the total sum of the areas/weights
long simulate_index(const long& m, const vector<Interval<double> >& I, const double& sum_w);

// Build the proposal from the initial set of points S0 using the interpolation defined by the type ('u','e','l')
void build_proposal(double (*target_function)(const double&), vector<Interval<double> >& I,double& sum_w, vector<double>& f_S, vector<double>& S0, long& m, const double& s_min, const double& s_max, const char& t);

// Update the left tail using the slope and intercept of the interval next to it
void update_left_tail (double (*target_function)(const double&), vector<Interval<double> >& I, double& sum_w, long& j_LEFT, long& m, long& j_ntl, const long& n);

// Update the right tail using the slope and intercept of the interval next to it
void update_right_tail (double (*target_function)(const double&), vector<Interval<double> >& I, double& sum_w, long& j_RIGHT, long& m, long& j_ntr, const long& n);

// Update the proposal by inserting point x into interval I[j]
void update_proposal (double (*target_function)(const double&), const double& x, const long& j, const double& f_x, vector<Interval<double> >& I, double& sum_w, long& m, long& j_ntl, long& j_ntr, long& j_LEFT, long& j_RIGHT, const long& n, double& p_n, long& j_n, Interval<double>& I_n, const double& x_n, const double& s_min, const double& s_max, const char& t);

// Simulates from the target distribution using the IA2RMS algorithm
double IA2RMS(double (*target_function)(const double&), vector<double>& S0, const long& N, vector<double>& x, const double& s_min, const double& s_max, const char& t);


#endif /* IA2RMS_hpp */
