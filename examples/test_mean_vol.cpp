// Copyright (c) 2016 Université libre de Bruxelles, Belgium.
//
// This file is part of MeanVol.
//
// MeanVol is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// MeanVol is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with MeanVol,
// see <http://www.gnu.org/licenses/>.
//
// Author(s)     : Vissarion Fisikopoulos

#include <algorithm>
#include <parallel/algorithm>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <iostream>

#include "../include/mean_volume_functions.hpp"
#include "../include/volume_functions.hpp"



////////////////////////////////////////////////////////////////////////

//Test the algorithms
int main(const int argc, const char** argv){

    //////////////////////////////////////////////////////////////////  
    // SET the parameters
    
    int D,N,k;
    bool test_brute_force=false;
    
    if (argc == 5){
      D = atoi(argv[1]);
      N = atoi(argv[2]);
      k = atoi(argv[3]);
      test_brute_force = atoi(argv[4]); 
    } else {
      std::cout << "Usage " << argv[0] << " D N k (dimension, number of points, parameter k) " <<
                   std::endl;
      exit(1);
    }    
    
    //Set print options
    bool experiments=true; //print just the timing for experiments
    bool check=true; //print also the values to check 
    bool print=false; //debug print

    //Experiments    
    //for(D=2; D<10; D++){
    //  for(N=10; N<200; N+=10){
        //for(k=D+1; k<N/2; k++){    
    //    k=D+2;
    if(1){if(1){    
      
    //////////////////////////////////////////////////////////////////  
    // GENERATE the input points
    
    std::vector<Point> points;
    
    // Toy example 2D
/*
    NT a0[] = {0,0};
    NT a1[] = {5,0};
    NT a2[] = {0,5};
    NT a3[] = {5,5};
    NT a4[] = {7,6};
    
    std::vector<NT> b0 (a0, a0 + sizeof(a0) / sizeof(NT) );
    std::vector<NT> b1 (a1, a1 + sizeof(a1) / sizeof(NT) );
    std::vector<NT> b2 (a2, a2 + sizeof(a2) / sizeof(NT) );
    std::vector<NT> b3 (a3, a3 + sizeof(a3) / sizeof(NT) );
    std::vector<NT> b4 (a4, a4 + sizeof(a4) / sizeof(NT) );
    
    Point p0(D,b0.begin(),b0.end());
    points.push_back(p0);
    Point p1(D,b1.begin(),b1.end());
    points.push_back(p1);
    Point p2(D,b2.begin(),b2.end());
    points.push_back(p2);
    Point p3(D,b3.begin(),b3.end());
    points.push_back(p3);
    Point p4(D,b4.begin(),b4.end());
    points.push_back(p4);
    
*/    
    // Random generation
    //CGAL::Random_points_in_cube_d<Point> rand_it(D, 1.0);
    CGAL::Random_points_on_sphere_d<Point> rand_it(D, 1.0);
    CGAL::cpp11::copy_n(rand_it, N, std::back_inserter(points));
    
    // Assign indices to points
    size_t index=0;
    for (std::vector<Point>::iterator pit=points.begin();
         pit != points.end(); ++pit){
      //std::cout << *pit << "[" << index << "]"<< std::endl;
      pit->set_index(index++);
    }
  
    double t1 = (double)clock()/(double)CLOCKS_PER_SEC;
    
    //////////////////////////////////////////////////////////////////  
    // TEST k-mean volume algorithm
    
    NT mean_vol = MeanVol::mean_volume(points.begin(),points.end(),k);
     
    double t2 = (double)clock()/(double)CLOCKS_PER_SEC;
  

    //////////////////////////////////////////////////////////////////  
    // TEST brute force
    
    NT  mean_vol_test=0;
    if(test_brute_force){
      mean_vol_test = MeanVol::mean_volume_bf(points.begin(),points.end(),k);
    }
    
    double t3 = (double)clock()/(double)CLOCKS_PER_SEC;

    //////////////////////////////////////////////////////////////////  
    // DISPLAY results
    
    if(check){
      std::cout<<"\nTotal Mean "<<k<<"-Volume="<<mean_vol<<std::endl;
      std::cout<<"\nTotal Mean(double) "<<k<<"-Volume="<<CGAL::to_double(mean_vol)<<std::endl;
      std::cout<<"\n================================\n";
    }
    
    if(check){
      std::cout<<"\nBRUTE FORCE: Total Mean  "<<k<<"-Volume="<<mean_vol_test<<std::endl;
      std::cout<<"\nBRUTE FORCE: Total Mean(double)  "<<k<<"-Volume="<<CGAL::to_double(mean_vol_test)<<std::endl;
      std::cout<<"\nTime (mean) = "<<t2-t1<<std::endl;
      std::cout<<"\nTime (brute-force) = "<<t3-t2<<std::endl;
    }
    
    if(experiments) std::cout<<D<<" "
                             <<N<<" "
                             <<k<<" "
                             <<t2-t1<<" "
                             <<t3-t2<<" "
                             <<MeanVol::global_compare
                             <<std::endl;
    //}
    }
    }   
    return 0;

}
