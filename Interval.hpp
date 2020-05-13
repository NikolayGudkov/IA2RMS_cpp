//
//  Interval.hpp
//
//  Created by Nikolay Gudkov on 20/4/20.
//
// Here, we define the class Interval. We assign several private and public members to it. In particular, we are interested in such properties of an interval as its boundaries, values of the target and proposal distribution, area/weight associated with the interval.
// Objects of this class are used as building blocks to construct the proposal distribution which approximates the target distribution.

#ifndef Interval_hpp
#define Interval_hpp
#include <iostream>
#include <limits>
#include "unif.hpp"

using namespace std;

const double inf = numeric_limits<double>::infinity();

template <class T>
class Interval{
private:
    T s_l;                  // left boudary of the interval
    T s_r;                  // right boundary of the interval
    T f_l;                  // value of the target at the left boundary (we assume that target is given in the log-scale)
    T f_r;                  // value of the target at the right boundary (we assume that target is given in the log-scale)
    T l;                    // value of the intercept of the proposal line
    T k;                    // value of the slope of the proposal line
    T w;                    // "unnormalised" weight (area) of the interval
    
    char t;                 // type of the interval 'u' uniform (flat), 'l' linear (logarithmic in log-scale), 'e' exponential (linear                         in the log-scale)
    bool left=false;        // this is a flag that indicates that the interval is left in the grid. The default value is false
    bool right=false;       // this is a flag that indicates that the interval is right in the grid. The default value is false
    
public:
    // Default constructor
    Interval() {t='u',s_l=0.; s_r=0.; f_l=0.; f_r=0.; l=0.; k=0.; w=0.;}
    
    // Constructor of the interval using parameters: type (t_c), left boundary (s_l_c), right boundary (s_r_c), values of the target on the boundaries of the interval f_l_c (left boundary) and f_r_c (right boundary)
    Interval(const char& t_c, const T& s_l_c, const T& s_r_c, const T& f_l_c, const T& f_r_c) {
        s_l=s_l_c; s_r=s_r_c; f_l=f_l_c; f_r=f_r_c;
        // This is the special case when the interval is far in the tail
        if (abs(f_l)==inf&&abs(f_r)==inf) {t='u'; k=0.; l=max(f_l,f_r); w=0.; return;}
        // Uniform interval
        if ((t_c!='e'&&t_c!='l')||f_l==f_r) {t='u'; k=0.; l=max(f_l,f_r); w=exp(l)*(s_r-s_l); return;}
        // Exponential interval (linear in log-scale)
        if (t_c=='e') {t=t_c; k=(f_r-f_l)/(s_r-s_l); l=f_l_c-k*s_l; w=(exp(f_r)-exp(f_l))/k; return;}
        // Linear interval (logarithmic in log-scale)
        if (t_c=='l') {t=t_c; k=(exp(f_r)-exp(f_l))/(s_r-s_l); l=exp(f_r)-k*s_r; w=0.5*(s_r-s_l)*(exp(f_r)+exp(f_l)); return;}
    }
    
    // Constructor of the interval using the parameter: left boundary (s_l_c), right boundary (s_r_c), values of the target on the boundaries f_l_c (left boundary), f_r_c (right boundary), slope (k_c) and intercept (l_c). We will use this constructor for leftmost and the rightmost intervals which bound domain where the proposal is defined. These intervals may inherit slope (k) and intercept (l) from their neighbours if either s_l_c=-inf or s_r_c=inf.
    Interval(const T& s_l_c, const T& s_r_c, const T& f_l_c, const T& f_r_c, const T& k_c, const T& l_c) {
        t='e'; s_l=s_l_c; s_r=s_r_c; f_l=f_l_c; f_r=f_r_c; k=k_c; l=l_c; w=(exp(f_r)-exp(f_l))/k;
        if (abs(f_l)==inf&&abs(f_r)==inf) {t='u'; k=0.; l=max(f_l,f_r); w=0.;}
    }
    
    // Methods to get information about parameters of the interval
    char get_t() const {return t;}
    T get_s_l() const {return s_l;}
    T get_s_r() const {return s_r;}
    T get_f_l() const {return f_l;}
    T get_f_r() const {return f_r;}
    T get_k() const {return k;}
    T get_l() const {return l;}
    
    // Derive value of the proposal at point x
    T get_p(const T& x) const {
        if (!(*this).is_member(x))
        {cout<<"Point is not in the interval"<<endl; return -1;}
        else{
            // Exponential interval
            if (t=='e') {if (abs(k)==inf) return -inf; return (k*x+l);}
            // Linear interval
            if (t=='l') {if (abs(k)==inf) return -inf; return log(k*x+l);}
            // Uniform interval
            return l;
        }
    }
    
    T get_w() const {return w;}
    bool is_left() const {return left;}
    bool is_right() const {return right;}
    
    // Methods to set values of the parameters of the interval
    void set_t(const char& t_c) {t=t_c;}
    void set_s_l(const T& s_l_c) {s_l=s_l_c;}
    void set_s_r(const T& s_r_c) {s_r=s_r_c;}
    void set_k(const T& k_c) {k=k_c;}
    void set_l(const T& l_c) {l=l_c;}
    void set_f_l(const T& f_l_c) {f_l=f_l_c;}
    void set_f_r(const T& f_r_c) {f_r=f_r_c;}
    void set_right() {right=true;}
    void set_left() {left=true;}
    void release_right() {right=false;}
    void release_left() {left=false;}
    
    // Check if the interval is valid 1) s_l<=s_r; 2) "left" interval should have a positive slope k>0; 3) "right" interval should have a negative slope k<0; 4) weight/area should be nonnegative w>=0;
    bool is_valid() const {
        if (s_l>s_r) {
            cout<<"Left side is greater than right side."<<endl; return false;}
        if (left&&(k<=0.)) {
            cout<<"Left interval with non-positive slope."<<endl; return false;}
        if (right&&(k>=0.)) {
            cout<<"Right interval with non-negative slope."<<endl; return false;}
        if (w<0.) {
            cout<<"Weight of the interval is negative."<<endl; return false;}
        return true;
    }
    
    // Method to simulate a point from the interval
    T simulate() const {
        if (!(*this).is_valid()) {cout<<"Interval is not valid. Cannot simulate a point."<<endl; return 0;}
        else{
            // Simulate from an exponential interval
            if(t=='e')  {
                return log(exp(k*s_l)+randu()*(exp(k*s_r)-exp(k*s_l)))/k;
            }
            // Simulate from a linear interval
            if (t=='l') {
                T C=k*s_l*s_l/2.+l*s_l+w*randu();
                return (-l+sqrt(l*l+2.*k*C))/k;
            }
            // Simulate from a uniform interval
            return (s_l+randu()*(s_r-s_l));
        }
    }
    
    // Check if a point is in the interval
    bool is_member(const T& x) const {return ((s_l<=x)&&(x<=s_r));}

    // Override operator =
    void operator=(const Interval<T>& that){
        s_l=that.get_s_l();
        s_r=that.get_s_r();
        f_l=that.get_f_l();
        f_r=that.get_f_r();
        w=that.get_w();
        l=that.get_l();
        k=that.get_k();
        left=that.is_left();
        right=that.is_right();
        t=that.get_t();
    }

    // Override operator <
    bool operator<(const Interval<T>& that) const{ if (s_l<that.s_l) {return true;} return false;}
    
    // Override operator ==
    bool operator==(const Interval<T>& that) const{
        if((s_l==that.get_s_l())&&(s_r==that.get_s_r())&&(f_l==that.get_f_l())&&(f_r==that.get_f_r())&&(k==that.get_k())&&(l==that.get_l())&&(w==that.get_w())&&(t==that.get_t())&&(left==that.is_left())&&(right==that.is_right())) return true;
        return false;}
    
    // Override operator !=
    bool operator!=(const Interval<T>& that) const{ return !(*this==that);}
    
    // Insert point x into interval (s_l,s_r). This results in splitting this interval into two parts (s_l,x) and (x,s_r).
    // Find parameters of the new intervals. Old interval acquires parameters of the left intervals (s_l,x), function returns right interval (x,s_r)
    Interval<T> insert(const T& x, const T& f_x){
        Interval<T> I_r(t,x,s_r,f_x,f_r);
        s_r=x; f_r=f_x;
        if (t=='e'){
            k=(f_r-f_l)/(s_r-s_l); l=f_l-k*s_l; w=(exp(f_r)-exp(f_l))/k;
        }
        else {
            if (t=='l'){
                k=(exp(f_r)-exp(f_l))/(s_r-s_l); l=exp(f_l)-k*s_l; w=0.5*(s_r-s_l)*(exp(f_l)+exp(f_r));
            }
            else{
                k=0.; l=max(f_l,f_x); w=exp(l)*(s_r-s_l);
            }
        }
        if (!is_valid())         cout<<"Left sub-interval is invalid after a point has been inserted."<<endl;
        if (!I_r.is_valid())     cout<<"Right sub-interval is invalid after a point has been inserted."<<endl;
        return I_r;
    }
};

// Override << operator
template<class T>
inline ostream& operator<<(ostream& os, const Interval<T>& I){
    os<<"(t="<<I.get_t()<<", s_l="<<I.get_s_l()<<", s_r="<<I.get_s_r()<<", f_l="<<I.get_f_l()<<", f_r="<<I.get_f_r()<<", k="<<I.get_k()<<", l="<<I.get_l()<<" ,w="<<I.get_w()<<", left="<<I.is_left()<<", right="<<I.is_right()<<")";
    return os;
}

#endif /* Interval_hpp */
