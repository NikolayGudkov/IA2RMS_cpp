//
//  IA2RMS.cpp
//
//  Created by Nikolay Gudkov on 19/4/20.
//

#include "IA2RMS.hpp"

// These parameters are used in the update of left and right tails
const double tau=0.01, beta=0.95;

// This constant is used to reserve memory to hold vectors of intervals. Increase it if you expect that the number of intervals to be used to approximate the target is higher and want to avoid memory reallocations.
const int m_max=2000;

//#################/
// Select interval /
//#################/
long simulate_index(const long& m, const vector<Interval<double> >& I, const double& sum_w){
    
// Select an interval I[i] from a vector of intervals {I[0],...,I[m-1]} using their non-normalised weights (areas) {w_0,...,w_{m-1}} and their total sum, {sum_w}, which is equal to the total area under the proposal curve.

    double cs_w=I[0].get_w()/sum_w;
    double r=randu();
    long i=0;
    
    while (r>cs_w){
        i++;
        cs_w+=(I[i].get_w()/sum_w);
    }
    return i;
}

//#############################/
// Build proposal distribution /
//#############################/
void build_proposal(double (*target_function)(const double&), vector<Interval<double> >& I,double& sum_w, vector<double>& f_S, vector<double>& S0, long& m, const double& s_min, const double& s_max, const char& t){
    
// This function aims to construct the initial proposal distribution by filling the vector of intervals (I) and computing the total sum of areas/weights of all intervals sum_w.
// Input: pointer to the target function from which we aim to simulate random values, bounds of the target function domain s_min and s_max, the initial set of support points S0, type of "inner" intervals t, m=size of S0, m+1 is the number of intervals including the left and right intervals.
// Output: the function fills the vector of intervals I, a vector of target function values at the support points f_S, and total area under the proposal curve (sum of weights) sum_w.
// Note that on some occasions we may add points to S0. It happens if the side intervals have incorrect slopes. In these cases, S0, f_S and m will represent the parameters of the "extended" set of support points.
      
    // Compute values of target function at all initial support points SHOULD WE ADD -50 HERE????
    for (long i=0;i<m-1;i++) f_S.push_back(target_function(S0[i]));
    
    // Define the "left" interval (i.e. the one which contains s_min)
    if (s_min==-inf){ // If s_min is -inf. The slope of the left interval is defined by the slope of the interval next to it.
        double k_left=(f_S[1]-f_S[0])/(S0[1]-S0[0])*(1-beta*exp(-tau));
        
        // If the slope is not positive, there will be numerical issues, so we have to search for an alternative point
            if (k_left<=0){
                // cout<<"Left interval is not correct. Searching for a new interval"<<endl;
                   double ds=0.001, s_tmp=S0[0], k_tmp, f_tmp;
                   do{
                       s_tmp-=ds;
                       f_tmp=target_function(s_tmp);
                       k_tmp=(f_S[0]-f_tmp)/(S0[0]-s_tmp);
                   }
                   while(k_tmp<=0);
                
                // Adjust the slope and compute the intercept of the left interval
                k_left=k_tmp*(1-beta*exp(-tau));
                double l_left=f_tmp-s_tmp*k_left;
                
                I.push_back(Interval<double>(-inf,s_tmp,target_function(-inf),f_tmp,k_left,l_left));
                I.push_back(Interval<double>(t,s_tmp,S0[0],f_tmp,f_S[0],f_tmp));
                // Add s_tmp and f_tmp to the beginning of the vectors S0 and f_S. Note that for "good" initial points S0 we don't have to search for alternative point s_tmp. So this operation will not be performed.
                S0.insert(S0.begin(),s_tmp); f_S.insert(f_S.begin(),f_tmp);
                
                // Increase the length m
                m++;
             }
        // If the slope is positive then we accept the candidate as a left interval and add it to the container I
            else{
                double l_left=f_S[0]-S0[0]*k_left;
                I.push_back(Interval<double>(-inf,S0[0],target_function(-inf),f_S[0],k_left,l_left));
            }
    }
    else{ // If s_min is not -inf, we push the interval into the container I
        Interval<double> I0(t,s_min,S0[0],target_function(s_min),f_S[0]);
        if (I0.is_valid()) I.push_back(I0);
        else {cout<<"The left interval is not valid. Choose other support points"<<endl; return;}
    }
    
    // Set the left interval
    I[0].set_left();
    
    // Fill the container I with intervals defined on points {S0[0],...,S0[m-2]}
    Interval<double> Ii;
    for (long i=0;i<=m-3;i++){
        Ii=Interval<double>(t,S0[i],S0[i+1],f_S[i],f_S[i+1]);
        if (Ii.is_valid()) I.push_back(Ii);
        else {cout<<"The interval is not valid. Choose other support points"<<endl; return;}
    }

    // Define the "right" interval (i.e. the one which contains s_max)
    if (s_max==inf){ // In s_max is inf. The slope of the right interval is defined by the slope of the interval next to it.
        double k_right=(f_S[m-2]-f_S[m-3])/(S0[m-2]-S0[m-3])*(1-beta*exp(-tau));
               
               // If the slope is not negative, there will be numerical issues, so we have to search for an alternative point
                   if (k_right>=0){
                          double ds=0.001, s_tmp=S0[m-2], k_tmp, f_tmp;
                          do{
                              s_tmp+=ds;
                              f_tmp=target_function(s_tmp);
                              k_tmp=(f_tmp-f_S[m-2])/(s_tmp-S0[m-2]);
                          }
                          while(k_tmp>=0);
                       
                       k_right=k_tmp*(1-beta*exp(-tau));
                       double l_right=f_tmp-s_tmp*k_right;
                       
                       I.push_back(Interval<double>(t,S0[m-2],s_tmp,f_S[m-2],f_tmp));
                       I.push_back(Interval<double>(s_tmp,inf,f_tmp,target_function(inf),k_right,l_right));
                       
                       // Add s_tmp and f_tmp to the end of the vectors S0 and f_S. Note that for "good" initial points S0 we don't have to search for alternative point s_tmp. So this operation will not be performed.
                       S0.push_back(s_tmp); f_S.push_back(f_tmp);
                       
                       // Increase the length m
                       m++;
                    }
               // If the slope is negative then we accept the candidate as a right interval and add it to the container I
                   else{
                       double l_right=f_S[m-2]-S0[m-2]*k_right;
                       I.push_back(Interval<double>(S0[m-2],inf,f_S[m-2],target_function(inf),k_right,l_right));
                   }
           }
    else{  // If s_max is not inf, we push the interval into the container I
               Interval<double> Im(t,S0[m-2],s_max,f_S[m-2],target_function(s_max));
               if (Im.is_valid()) I.push_back(Im);
               else {cout<<"The left interval is not valid. Choose other support points"<<endl; return;}
    }
    
    // Set the right interval
    I[m-1].set_right();
    
    // Fill the vector of weights and comute the total sum
    sum_w=0.;
    for (long i=0;i<m;i++) sum_w+=I[i].get_w();
}

//##################/
// Update left tail /
//##################/
void update_left_tail (double (*target_function)(const double&), vector<Interval<double> >& I, double& sum_w, long& j_LEFT, long& m, long& j_ntl, const long& n){
    
// This function is called when the left interval is extended to -inf. It helps to update the left tail of the proposal distribution using the slope and intercept of the interval next to it.
    
    double k_left=(I[j_ntl].get_f_r()-I[j_ntl].get_f_l())/(I[j_ntl].get_s_r()-I[j_ntl].get_s_l())*(1-beta*exp(-tau*n));
    
    // If the slope is not positive, there will be numerical issues, so we have to search for an alternative point
         if (k_left<=0){
            double ds=0.001, s_tmp=I[j_ntl].get_s_l(), k_tmp, f_tmp;
            do{
                s_tmp-=ds;
                f_tmp=target_function(s_tmp);
                k_tmp=(I[j_ntl].get_f_l()-f_tmp)/(I[j_ntl].get_s_l()-s_tmp);
            }
            while(k_tmp<=0);
             
            // Adjust the slope and compute the intercept of the left interval
            k_left=k_tmp*(1-beta*exp(-tau*n));
            double l_left=f_tmp-s_tmp*k_left;
             
            // Update left interval
            I[j_LEFT]=Interval<double>(-inf,s_tmp,target_function(-inf),f_tmp,k_left,l_left);
             
            // Add new interval to the vector I
            I.push_back(Interval<double>(I[j_ntl].get_t(),s_tmp,I[j_ntl].get_s_l(),f_tmp,I[j_ntl].get_f_l()));
            j_ntl=m;
            sum_w+=I[m].get_w();
             
            // Increase the length m
            m++;
         }
         else{
            double l_left=I[j_ntl].get_f_l()-k_left*I[j_ntl].get_s_l();
            I[j_LEFT]=Interval<double>(-inf,I[j_ntl].get_s_l(),target_function(-inf),I[j_ntl].get_f_l(),k_left,l_left);
         }
    I[j_LEFT].set_left();
}

//###################/
// Update right tail /
//###################/
void update_right_tail (double (*target_function)(const double&), vector<Interval<double> >& I, double& sum_w, long& j_RIGHT, long& m, long& j_ntr, const long& n){
    
// This function is called when the right interval is extended to +inf. It helps to update the right piece of the proposal distribution using the slope and intercept of the interval next to it.

    double k_right=(I[j_ntr].get_f_r()-I[j_ntr].get_f_l())/(I[j_ntr].get_s_r()-I[j_ntr].get_s_l())*(1-beta*exp(-tau*n));

      // If the slope is not negative, there will be numerical issues, so we have to search for an alternative point
      if (k_right>=0){
          double ds=0.001, s_tmp=I[j_ntr].get_s_r(), k_tmp, f_tmp;
             do{
                 s_tmp+=ds;
                 f_tmp=target_function(s_tmp);
                 k_tmp=(f_tmp-I[j_ntr].get_f_r())/(s_tmp-I[j_ntr].get_s_r());
             }
             while(k_tmp>=0);
          
          // Adjust the slope and compute the intercept of the right interval
          k_right=k_tmp*(1-beta*exp(-tau));
          double l_right=f_tmp-s_tmp*k_right;
          
          // Update right interval
          I[j_RIGHT]=Interval<double>(s_tmp,inf,f_tmp,target_function(inf),k_right,l_right);
          
          // Add new interval to the vector I
          I.push_back(Interval<double>(I[j_ntr].get_t(),I[j_ntr].get_s_r(),s_tmp,I[j_ntr].get_f_r(),f_tmp));
          
          // Increase the length m
          m++;
          
          j_ntr=m;
          sum_w+=I[m].get_w();
          
      }
      else{
          double l_right=I[j_ntr].get_f_r()-k_right*I[j_ntr].get_s_r();
          I[j_RIGHT]=Interval<double>(I[j_ntr].get_s_r(),inf,I[j_ntr].get_f_r(),target_function(inf),k_right,l_right);
      }
      I[j_RIGHT].set_right();
}

//###############################/
// Updated proposal distribution /
//###############################/
void update_proposal (double (*target_function)(const double&), const double& x, const long& j, const double& f_x, vector<Interval<double> >& I, double& sum_w, long& m, long& j_ntl, long& j_ntr, long& j_LEFT, long& j_RIGHT, const long& n, double& p_n, long& j_n, Interval<double>& I_n, const double& x_n, const double& s_min, const double& s_max, const char& t){

// This function updates the proposal function by inserting point x into interval I[j]. In addition to this it updates the left and right intervals when those are extended to -inf and +inf respectively.

    // Since we update the interval I[j], its weight will change, so we substract it from the total sum
    sum_w-=I[j].get_w();
    
    // Update internal interval
    if (!I[j].is_left()&&!I[j].is_right()){
        
        // Update I[j] by changing its parameters and puting the right cut of the interval j-(x_st,s_r[j]) at the back of the vector I
        I.push_back(I[j].insert(x,f_x));
        
        // Add the weights of the new intervals to the total sum
        sum_w+=(I[m].get_w()+I[j].get_w());

        // Check if x_n is located in either of the "updated" I[j] and I[m] intervals and update p_n accordingly.
        // Note that f_n is not affected by changes in the interval I[j]
        if (j_n==j){
            if (I[j].is_member(x_n)) {p_n=I[j].get_p(x_n); j_n=j; I_n=I[j];}
            else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}}
        }
        
        // In case we have updated the interval located next to the right interval, the right cut is moved to the end of the array, so its index is changed to m. We need to keep track of it because it provides its slope and intercept parameters to the right interval, which contains inf.
        if (j==j_ntr) j_ntr=m;
           
        // Update slope and intercept of the left interval if it contains -inf
        if (I[j_LEFT].get_s_l()==-inf) {
            // Since we update the interval I[j_LEFT], its weight will change, so we substract it from the total sum
            sum_w-=I[j_LEFT].get_w();
            update_left_tail (target_function, I, sum_w, j_LEFT, m, j_ntl, n);
            sum_w+=I[j_LEFT].get_w();
        }

        if (j_n==j_LEFT){
                 if (I[j_LEFT].is_member(x_n)) {p_n=I[j_LEFT].get_p(x_n); j_n=j_LEFT; I_n=I[j_LEFT];}
                 else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}}
        }
        
        // Update slope and intercept of the right interval if it contains -inf
        if (I[j_RIGHT].get_s_r()==inf) {
            // Since we update the interval I[j_RIGHT], its weight will change, so we substract it from the total sum
            sum_w-=I[j_RIGHT].get_w();
            update_right_tail (target_function, I, sum_w, j_RIGHT, m, j_ntr, n);
            sum_w+=I[j_RIGHT].get_w();
        }
        
        if (j_n==j_RIGHT){
                 if (I[j_RIGHT].is_member(x_n)) {p_n=I[j_RIGHT].get_p(x_n); j_n=j_RIGHT; I_n=I[j_RIGHT];}
                 else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}}
        }
    }
    // Update the side intervals
    else{
        // Update the left interval j=j_LEFT
        if (I[j].is_left()){
            // Create new next to left interval an push in to the back of the vector I
            I.push_back(Interval<double>(I[j_ntl].get_t(),x,I[j].get_s_r(),f_x,I[j].get_f_r()));
            j_ntl=m;
            sum_w+=I[m].get_w();
            
            if (I[j_LEFT].get_s_l()==-inf){
                update_left_tail(target_function, I, sum_w, j_LEFT, m, j_ntl, n);
            }
            else{
                Interval<double> Ij(I[j_ntl].get_t(),I[j].get_s_l(),x,I[j].get_f_l(),f_x);
                I[j]=Ij; I[j].set_left();
            }
             sum_w+=I[j_LEFT].get_w(); j_LEFT=j;
            
            // Check if x_n is located in either of the "updated" I[j] and I[m] intervals and update p_n accordingly.
              // Note that f_n is not affected by changes in the interval I[j]
              if (j_n==j){
                  if (I[j].is_member(x_n)) {p_n=I[j].get_p(x_n); j_n=j; I_n=I[j];}
                  else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}}
              }
            
            // Update slope and intercept of the right interval if it contains -inf
            if (I[j_RIGHT].get_s_r()==inf) {
                // Since we update the interval I[j_RIGHT], its weight will change, so we substract it from the total sum
                sum_w-=I[j_RIGHT].get_w();
                update_right_tail (target_function, I, sum_w, j_RIGHT, m, j_ntr, n);
                sum_w+=I[j_RIGHT].get_w();
            }
            
            if (j_n==j_RIGHT){
                     if (I[j_RIGHT].is_member(x_n)) {p_n=I[j_RIGHT].get_p(x_n); j_n=j_RIGHT; I_n=I[j_RIGHT];}
                     else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}}
            }

        }
        // Update the right interval j=j_RIGHT
        else{
            Interval<double> Ij(I[j_ntr].get_t(),I[j].get_s_l(),x,I[j].get_f_l(),f_x);

            if (I[j_RIGHT].get_s_r()==inf){
                double k_right=(Ij.get_f_r()-Ij.get_f_l())/(Ij.get_s_r()-Ij.get_s_l())*(1-beta*exp(-tau*n));
                
                // If the slope is not negative, there will be numerical issues, so we have to search for an alternative point
                    if (k_right>=0){
                        double ds=0.001, s_tmp=x, k_tmp, f_tmp;
                           do{
                               s_tmp+=ds;
                               f_tmp=target_function(s_tmp);
                               k_tmp=(f_tmp-I[j_ntr].get_f_r())/(s_tmp-I[j_ntr].get_s_r());
                           }
                           while(k_tmp>=0);
                        
                        // Adjust the slope and compute the intercept of the right interval
                        k_right=k_tmp*(1-beta*exp(-tau));
                        double l_right=f_tmp-s_tmp*k_right;
                        
                        // Add new interval to the vector I
                        I.push_back(Interval<double>(I[j_ntr].get_t(),I[j_ntr].get_s_r(),s_tmp,I[j_ntr].get_f_r(),f_tmp));
                        j_ntr=m;
                        sum_w+=I[m].get_w();
                        
                        // Increase the length m
                        m++;
                        
                        // Update right interval
                        I.push_back(Interval<double>(s_tmp,inf,f_tmp,target_function(inf),k_right,l_right));
                    }
                    else{
                        double l_right=f_x-k_right*x;
                        I.push_back(Interval<double>(x,inf,f_x,target_function(inf),k_right,l_right));
                        j_ntr=j;
                    }
            }
            else{
                I.push_back(Interval<double>(I[j].get_t(),x,I[j].get_s_r(),f_x,I[j].get_f_r()));
                j_ntr=j;
            }
            I[j]=Ij; j_RIGHT=m; I[m].set_right(); sum_w+=(I[m].get_w()+I[j].get_w());

           if (j_n==j){
                     if (I[j].is_member(x_n)) {p_n=I[j].get_p(x_n); j_n=j; I_n=I[j];}
                     if (I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}
                     if (I[m-1].is_member(x_n)) {p_n=I[m-1].get_p(x_n); j_n=m-1; I_n=I[m-1];}
            }

            // Update slope and intercept of the left interval if it contains -inf
            if (I[j_LEFT].get_s_l()==-inf) {
                // Since we update the interval I[j_LEFT], its weight will change, so we substract it from the total sum
                sum_w-=I[j_LEFT].get_w();
                update_left_tail (target_function,I, sum_w, j_LEFT, m, j_ntl, n);
                sum_w+=I[j_LEFT].get_w();
            }
            
            if (j_n==j_LEFT){
                     if (I[j_LEFT].is_member(x_n)) {p_n=I[j_LEFT].get_p(x_n); j_n=j_LEFT; I_n=I[j_LEFT];}
                     else{if(I[m].is_member(x_n)) {p_n=I[m].get_p(x_n); j_n=m; I_n=I[m];}
                     }
            }
        }
    }
}

//###################/
// IA2RMS simulation /
//###################/
double IA2RMS(double (*target_function)(const double&), vector<double>& S0, const long& N, vector<double>& x, const double& s_min, const double& s_max, const char& t){
    
    // Samples from the target distribution function via the IA2RMS algorithm
    //    target_function - the log-scale probability distribution function which we know up to the constant of proportionality
    //    S0 - vector initial set of support points
    //    N  - number of sample points
    //    x  - vector where we store sample points
    //    s_min and s_max - bounds of the target function domain
    //    t - character which defined the type of inner intervals ('u' uniform, 'e' exponential, 'l' linear)
    
    // Count the number of intervals in the support set s
    long m=S0.size()+1;
    
    vector<double> f_S; f_S.reserve(m);
    
    // Total sum of non-normalised weights (areas) of all intervals
    double sum_w=0.;
    
    // Define a container for intervals under the proposal distribution and reserve memory to hold it
    vector<Interval<double> > I; I.reserve(m_max);
    
    build_proposal(target_function, I, sum_w, f_S, S0, m, s_min, s_max, t);

    long j_LEFT=0, j_ntl=1, j_ntr=m-2, j_RIGHT=m-1; //define indices of the intervals which contain s_min (LEFT) and s_max (RIGHT) and also intervals next to them ntl (next to the left) and ntr (next to the right) which we are going to use to determine slopes of the LEFT and RIGHT intervals
    
    // Simulate initial point
    long j=simulate_index(m,I,sum_w); //select the interval
    x[0]=I[j].simulate(); //simulated initial point
    unsigned long n=0; // counts how many points we have simulated minus 1
    
    // Evaluate the target and proposal functions at x[0]
      double f_n=target_function(x[0]), p_n=I[j].get_p(x[0]), x_st, p_st, y, alp, f_st, f_y;
      long j_y, j_n=j;
    
    // Define objects of type Interval which will store:
    //  the interval I_st is used to sample the current candidate (x_st)
    //  the interval I_y is the one, where the auxiliary variable (y) is located
    //  the interval I_n is the one, where the current simulated value (x_n) is located
    Interval<double> I_y, I_n=I[j];
  
   // simulate values from the proposal which approximate the target distribution using the IA2RMS algorithm
while(n<N-1){
        
        // Simulate candidate x_st
        j=simulate_index(m,I,sum_w);        // Select the interval from which we sample a candidate
        x_st=I[j].simulate();               // Simulate value from the selected interval
        f_st=target_function(x_st);         // Value of the tagret at x_st
        p_st=I[j].get_p(x_st);              // Value of the proposal at simulated point

        //%%%%%%%%%%%%%%%%%%%%
        //% Rejection test %%%
        //%%%%%%%%%%%%%%%%%%%%
            if (randu()>exp(f_st-p_st)){    // If the candidate is rejected we insert it in the support grid s.
                                            // For this, we have to calculate parameters of the support function in the two new intervals I[j]=[x_j, x_st] and I[m]=[x_st, x_j+1]
                        
            // Insert x_st into interval j and update the proposal accordingly
                update_proposal(target_function,x_st,j,f_st,I,sum_w,m,j_ntl,j_ntr,j_LEFT,j_RIGHT,n,p_n,j_n,I_n,x[n], s_min, s_max, t);
        
                m++;   //count new point in the grid
                
            }
            else {
            // If the candidate x_st is not rejected we perform the MH test and decide if we accept it as a new state or we keep the old value
            //%%%%%%%%%%%%%%%%%
            //%   MH test     %
            //%%%%%%%%%%%%%%%%%
                alp=exp((f_st-p_st)*(f_st>=p_st)+(p_n-f_n)*(p_n<=f_n));
                if (randu()<=min(1.,alp)){
                    y=x[n]; I_y=I_n; f_y=f_n; j_y=j_n;
                    x[n+1]=x_st; f_n=f_st; I_n=I[j]; j_n=j; p_n=p_st;
                }
                else{
                     y=x_st; f_y=f_st; I_y=I[j]; j_y=j; x[n+1]=x[n];
                }
                n++;                      // count new simulated point

            //%%%%%%%%%%%%%%%%%
            //% Control test  %
            //%%%%%%%%%%%%%%%%%
            if (randu()>exp(I_y.get_p(y)-f_y)){  // Test if the auxilirary point y is added to the grid
     
                // Insert x_st into interval j and update the proposal accordingly
                update_proposal(target_function,y,j_y,f_y,I,sum_w,m,j_ntl,j_ntr,j_LEFT,j_RIGHT,n, p_n, j_n, I_n, x[n], s_min, s_max, t);
                //count new point in the grid
                m++;
            }
            }
}

// Indicator if we want to plot the target probability function as well as the proposal distribution constructed.
    if (0){
// Since the constructed intervals are stored in the vector I not in the sorted order, we sort them before plotting
    sort(I.begin(),I.end());
        
// Define the grid of points which are used to plot the density curves
    long N_sim=100000;
    double S1=-10.0+I[0].get_s_r(), S2=10.0+I[m-1].get_s_l(), ds=(S2-S1)/(N_sim-1);
        
    j=0;
    ofstream dens;
    // We store the generated points in the file. We can use for example gnuplot to read this file and plot curves using the values stored there.
    dens.open("simulated_density.txt");
    for (long i=0;i<N_sim;i++){
        // Store the values of the grid (horizontal axis)
        dens<<S1+i*ds<<' ';
        
        // Store the values of the proposal constructed for a given interval
        if(I[j].is_member(S1+i*ds)) dens<<exp(I[j].get_p(S1+i*ds))<<' ';
        // Move to the next interval if the point is not in the current interval
        else{do {++j;} while(!I[j].is_member(S1+i*ds)); dens<<exp(I[j].get_p(S1+i*ds))<<' ';}
        
        // Store the values of the target
        dens<<exp(target_function(S1+i*ds))<<endl;
    }
    dens.close();
   }

    return 1/sum_w;
}
