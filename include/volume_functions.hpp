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

#ifndef VOLUME_FUNCTIONS_H
#define VOLUME_FUNCTIONS_H

#include <CGAL/Cartesian_d.h>  
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>

#include "external/combination.h"

typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag >    K_filtered;
typedef CGAL::Triangulation<K_filtered>                 Triangulation;

namespace MeanVol{

template <class Point>
int convert_copy_point(Point p1, Triangulation::Point& p2){
  
  std::vector<double> coords;
  for (CGAL::Cartesian_d<NT>::Cartesian_const_iterator_d pit=p1.cartesian_begin(); pit!=p1.cartesian_end();++pit){
    coords.push_back(CGAL::to_double(*pit));
  }
  p2 = Triangulation::Point(coords.begin(),coords.end());
  return 1;
}

bool is_generic(){
  return true;
}

//Compute the volume of a simple by computing the determinant
template<class FT, class PointInputIterator>
FT simplex_vol(PointInputIterator pfirst,
                   PointInputIterator plast){
  
  typedef Eigen::Matrix<FT,Eigen::Dynamic,Eigen::Dynamic> MT;
  int D=pfirst->dimension();
  MT m(D,D);
  
  PointInputIterator s = pfirst;
  ++s;
  for( int j = 0; j < D; ++s, ++j )
    for( int i = 0; i < D; ++i )
      m(i,j) = s->cartesian(i) - pfirst->cartesian(i);
  
  FT det = m.determinant();
  
  // Return the absolute value
  det = det<0 ? -1*det : det;
  
  return det/factorial(D);
}


//Compute the volume of a polytope by computing a triangulation and 
//adding the volumes of simplices
template<class FT, class PointInputIterator>
FT volume(PointInputIterator pfirst,
          PointInputIterator plast){
  
    int D = pfirst->dimension();
    Triangulation t(D);                      // create triangulation
    CGAL_assertion(t.empty());
    
    t.insert(pfirst, plast);  // compute triangulation
    CGAL_assertion( t.is_valid() );
          
    int i=0;
    typedef Triangulation::Finite_full_cell_iterator Finite_full_cell_iterator;
    
    FT vol=0;
    
    for(Finite_full_cell_iterator cit = t.finite_full_cells_begin();
         cit != t.finite_full_cells_end(); ++cit )
      {
        std::vector<Triangulation::Point> cell_points;
        
        for(Triangulation::Triangulation_ds::Full_cell::Vertex_handle_iterator vit = cit->vertices_begin(); vit!=cit->vertices_end(); vit++){
              //std::cout<<**vit<<std::endl;
              cell_points.push_back((**vit).point());
        }
        vol += simplex_vol<NT>(cell_points.begin(),cell_points.end());
        
        ++i;
      }
    //std::cout << "There are " << i << " cells on the convex hull."<< std::endl;
    //std::cout << "The volume is " << vol << std::endl;
    
    return vol;
}

// Brute force algorithm for computing the k-mean volume 
template<class PointInputIterator>
NT mean_volume_bf(PointInputIterator pfirst,PointInputIterator plast, int k){
  
      int D = pfirst->dimension();
      int N = std::distance(pfirst,plast);
      
      NT mean_vol=0; //the mean volume of K-sets
      
      // - - - - - - - - - - 
      // Generate the input points, copy original points
      std::vector<Triangulation::Point> points_t;
      for(std::vector<Point>::iterator pit_original=pfirst; pit_original!=plast; pit_original++){      
        Triangulation::Point p_temp;
        convert_copy_point(*pit_original,p_temp);
        points_t.push_back(p_temp);
      }
  
      std::vector<int> cat(N) ; // vector with ints.
      for(int i=0; i<N; i++){cat[i]=i;}
      
      std::vector<int> cbt(k) ; // vector with ints.
      for(int i=0; i<k; i++){cbt[i]=i;}
      
      //Iterate through all possible combinations of D-1 points
      do{
        std::vector<Triangulation::Point> selection_points;
        for(std::vector<int>::iterator it=cbt.begin(); it<cbt.end(); it++){
          selection_points.push_back(points_t[*it]);
        }
        //std::cout<<"\nSelected points= ";
        //display(selection_points.begin(),selection_points.end());
        
        NT vol = volume<NT>(selection_points.begin(),selection_points.end());
        //std::cout<<"volume="<<vol<<std::endl;
        mean_vol+=vol;
      }
      while(stdcomb::next_combination(cat.begin(),cat.end(),cbt.begin(),cbt.end()));
      
      return mean_vol;
}

} // namespace MeanVol

#endif // VOLUME_FUNCTIONS_H
