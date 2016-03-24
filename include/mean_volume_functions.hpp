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

#ifndef MEAN_VOLUME_FUNCTIONS_H
#define MEAN_VOLUME_FUNCTIONS_H

#include <CGAL/algorithm.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
//#include <CGAL/Gmpz.h>

#include <Eigen/Eigen>

#include "external/combination.h"

//typedef CGAL::Gmpq                                      NT;
typedef double                                          NT;
typedef CGAL::Cartesian_d<NT>                           KE;

typedef KE::Point_d                                     Point;
typedef KE::Vector_d                                    Vector;
typedef KE::Hyperplane_d                                Hyperplane;
typedef KE::Orientation_d                               Orientation;


namespace MeanVol{
  
//factorial function
template<class FT>
FT factorial(FT x){
  if(x<2) return 1;
  FT fac=1;
  for(FT i=2; i<=x; i+=1){
    fac*=i;
  }
  return fac;
}

//number of combinations
template<class FT>
FT comb(FT n, FT r){
  if(n==0) return 0;
  if(n<r) return 0;
  return factorial(n)/(factorial(n-r)*factorial(r));
}

// for use with next_combination examples!
template<class PointInputIterator>
void display(PointInputIterator pfirst,PointInputIterator plast)
{
  for (PointInputIterator it=pfirst;it!=plast;++it)
      std::cout<<*it<<" ";
  std::cout<<std::endl;
}


//Test if a point is in the list
template<class PointInputIterator, class Point_gen>
int not_in_selection_points(PointInputIterator pfirst,
                            PointInputIterator plast,
                            Point_gen p){
  for (PointInputIterator it=pfirst;it!=plast;++it)
      if(*it==p)
        return 0;
  return 1;
}          

   
int global_compare=0;

// sort using a custom function object based on orientation along a subspace
// defined by the set of points in selection_points
template<class FT>
struct compare_points{
  
    compare_points(std::vector<Point> points, Point ref) : 
                      selection_points(points),
                      ref(ref){};
  
    bool operator()(Point p, Point q)
    {   
        global_compare++;
        if(p==ref) return true; //ref is the largest element
        if(q==ref) return false; //ref is the largest element
        
        std::vector<Point> points_with_p_ref(selection_points);
        points_with_p_ref.push_back(p);
        points_with_p_ref.push_back(ref);
        FT ori_p_ref = Orientation()(points_with_p_ref.begin(),points_with_p_ref.end());
        
        //std::cout<<p<<","<<ref<<"="<<ori_p_ref<<std::endl;
        
        std::vector<Point> points_with_q_ref(selection_points);
        points_with_q_ref.push_back(q);
        points_with_q_ref.push_back(ref);
        FT ori_q_ref = Orientation()(points_with_q_ref.begin(),points_with_q_ref.end());
        
        //std::cout<<q<<","<<ref<<"="<<ori_q_ref<<std::endl;
        
        std::vector<Point> points_with_pq(selection_points);
        points_with_pq.push_back(p);
        points_with_pq.push_back(q);
        FT ori_pq = Orientation()(points_with_pq.begin(),points_with_pq.end());
        
        //std::cout<<p<<","<<q<<"="<<ori_pq<<std::endl;
        //std::cout<<std::endl;
        
        if(p!=q && ori_pq==0 && !p.get_is_reflected() && !q.get_is_reflected()){
          std::cout<<*selection_points.begin()<<","<<p<<","<<q<<"="<<ori_pq<<std::endl;
        
          std::cout<<"INPUT NOT GENERIC"<<std::endl;
          exit(0);
        }
        
        //if from the same side of ref hyperplane
        if(ori_p_ref*ori_q_ref>=0) return ori_pq>0;
        
        return ori_p_ref<0;
        
    }   
    
    private:
      std::vector<Point> selection_points;
      Point ref;
};

////Barycenter
//std::cout<<"\n\nPoints="<<std::endl;
//display(selection_points.begin(),selection_points.end());
//Point barycenter = barycenter_d(selection_points.begin(),selection_points.end()); 
//std::cout<<"Barycenter="<<barycenter<<std::endl;

//Compute the barycenter of a d-point
template <class PointInputIterator>
Point barycenter_d(PointInputIterator pfirst,
               PointInputIterator plast){
  PointInputIterator it=pfirst;
  Vector barycenter = *it - CGAL::Origin();
  ++it;
  for (;it!=plast;++it){
    barycenter += *it - CGAL::Origin(); 
  }
  barycenter/=std::distance(pfirst,plast);
  return CGAL::Origin()+barycenter;
}

//Compute cross product as Laplace expansion of a matrix determinant
template<class FT, class PointInputIterator>
Vector cross_product(PointInputIterator pfirst,
                     PointInputIterator plast){
  //std::cout<<"cross product: ";
  //display(pfirst,plast);
  
  //Create eigen matrix and compute the determinant
  //TODO: do not construct the matrix every time
  typedef Eigen::Matrix<FT,Eigen::Dynamic,Eigen::Dynamic> MT;
  int D=pfirst->dimension();
  MT m(D-1,D-1);
  //std::cout << "d=" << D << std::endl << std::endl;
  
  //The matrix is constructed by subtitution of the first point from all
  //the other points and then add one row with the coefficients of the 
  //cross-product. The expansion happens along the later row. 
  std::vector<FT> cross_vec;  
  for(int idx=0; idx<D; idx++){
    PointInputIterator s = pfirst;
    ++s;
    for(int j=0; j<D-1; ++s,++j){
      for(int i=0, coord=0; i<D;  i++){
        if(i!=idx){
          m(coord++,j) = s->cartesian(i) - pfirst->cartesian(i);
          //std::cout<<"("<<coord<<","<<j<<")="<<s->cartesian(i) - pfirst->cartesian(i)<<std::endl;
        }
      }
    }
    cross_vec.push_back(m.determinant()*std::pow(-1,idx));
  }
  return Point(D,cross_vec.begin(),cross_vec.end())-CGAL::Origin();
}

// Algorithm for computing the k-mean volume of a set of points   
template<class PointInputIterator>
NT mean_volume(PointInputIterator pfirst,PointInputIterator plast, int k){
  
    int D = pfirst->dimension();
    int N = std::distance(pfirst,plast);
      
    NT mean_vol=0; //the mean volume of K-hulls
    
    bool print=false;
    
    // - - - - - - - - - - 
    // Iterate through combinations N choose D+1 points and 
    // sort by angle
    
    std::vector<bool> select(N);
    std::fill(select.begin(), select.end() - N + D-1, true);

    // Compute combinations using std::permutations
    //do {
    //   for (int i = 0; i < N; ++i) {
    //       if (select[i]) {
    //           std::cout << points[i] << " ";
    //       }
    //   }
    //   std::cout << "\n";
    //} while (std::prev_permutation(select.begin(), select.end()));
    
    // Maybe more efficient way for combinations
    std::vector<int> ca(N) ; // vector with ints.
    for(int i=0; i<N; i++){ca[i]=i;}
    
    std::vector<int> cb(D-1) ; // vector with ints.
    for(int i=0; i<D-1; i++){cb[i]=i;}
     
    //Iterate through all possible combinations of D-1 points
    do{
      std::vector<Point> selection_points;
      for(std::vector<int>::iterator it=cb.begin(); it<cb.end(); it++){
        selection_points.push_back(*(pfirst+*it));
      }
      if(print){
        std::cout<<"Selected points\n";
        display(selection_points.begin(),selection_points.end());
        std::cout<<"----"<<std::endl;
      }
      //Construct rest_of_points list
      //Add reflection points and use annotation to declare 
      //between reflection and original points
      //Do not consider points that belongs to selection points
      std::vector<Point> rest_of_points;
      std::vector<Point>::iterator pit=pfirst;
      //std::vector<int>::iterator it=cb.end();
      //--it;
      //std::advance(pit,*it+1);//start after the last selection point 
                              //to avoid duplicates

      //if(pit==points.end()){break;}
      for(; pit!=plast; pit++){
        if(not_in_selection_points(selection_points.begin(),
                                   selection_points.end(),                           
                                   *pit)){
          Point p_original(*pit);
          p_original.set_is_reflected(false);
          rest_of_points.push_back(p_original);
          
          //Compute reflection point wrt selection points
          Point p_reflected(*selection_points.begin());
          Vector v_reflected = *pit-p_reflected;
          v_reflected*=-1;
          p_reflected+=v_reflected;
          p_reflected.set_is_reflected(true);
          rest_of_points.push_back(p_reflected);
        }
      }
      
      //std::cout<<"Rest of points\n";
      //for (std::vector<Point>::iterator it=rest_of_points.begin();it!=rest_of_points.end();++it)
      //  std::cout<<*it<<" ";
      //std::cout<<"----"<<std::endl;
      
      //Shuffle vector before sorting
      //std::srand ( unsigned ( std::time(0) ) );
      //std::random_shuffle(rest_of_points.begin(), rest_of_points.end());
      
      // SORTING
      //First find a reference point for sorting; it should not be reflected
      std::vector<Point>::iterator refp = rest_of_points.begin(); 
      while(refp->get_is_reflected()){refp++;}
      if(print) std::cout<<"reference point="<<*refp<<std::endl;
      Point refp_reflection(*(++refp));
      if(print) std::cout<<"reflection of reference point="<<refp_reflection<<std::endl;
      refp--;
      
      //Sort points according to angle wrt selected set of points
      std::sort(rest_of_points.begin(), rest_of_points.end(), 
                compare_points<NT>(selection_points,*refp));
      //__gnu_parallel::sort(rest_of_points.begin(), rest_of_points.end(), 
      //          compare_points<NT>(selection_points,*refp));
      
      //Print points
      //for(std::vector<Point>::iterator pit=rest_of_points.begin(); pit!=rest_of_points.end(); pit++){
      //  std::cout<<*pit<<std::endl;
      //}
      //std::cout<<"-----------------------------------\n";

      
      //COUNTING 
      //STEP 1: How many point are separated wrt to reference point/hyperplane
      int total_count=0;
      {
        std::vector<Point>::iterator pit=rest_of_points.begin();
        while(*pit != refp_reflection){
          ++pit;
          if(!pit->get_is_reflected()) ++total_count;
        }
      }
      //std::cout<<"Polytopes count start="<<total_count<<"\n";
     
      //STEP 2: Update total_count
      //Geometrically we rotate (ccw) a hyperplane around selection points
      //Events occur at original points and we count the reflection 
      //points betweem two events 
      //Every point in rest_of_points together with selection points
      //define a facet for which we want to compute the normal vector and 
      //it contribution to the polytope volume
      //We use the divergence thoerem to compute the volume of the polytope
      //as a sum (over facet normal vector) of dot products 
      int refl_count = 1;
      //Point reference_point = *rest_of_points.begin();
      Point last_selection_point = selection_points.back();
      std::vector<Point>::iterator pait = rest_of_points.begin();
      for( ;pait!=rest_of_points.end(); pait++){
        NT vol;
        //Compute facet contribution to volume
        if(!pait->get_is_reflected()){//for original points
          //Add the current point
          selection_points.push_back(*pait);
          if(print) std::cout<<"\nCurrent="<<*pait<<std::endl;
          Vector norm = cross_product<NT>(selection_points.begin(),selection_points.end());
          vol = ((*(selection_points.begin())-CGAL::Origin())*norm)/factorial(D);
          if(print) std::cout<<vol<<std::endl;
          //Remove current point
          selection_points.pop_back();
        }  
        
        //Count in how many polytopes the current facet belongs
        //if(pait!=rest_of_points.begin()){//skip the reference point
          if(pait->get_is_reflected()){//count the reflected points
            refl_count++;
          }else{                      //stop at an original point
            //total_count: counts the number of points separeted by
            //reference point and have positive ccw orientation
            total_count = total_count + refl_count - 1; 
            refl_count=0;
            if(print){
              std::cout<<"Polytopes="<<total_count<<",("
                     <<comb(NT(total_count),NT(k-D))<<")\n";
              std::cout<<"Indices="<<pait->index()<<","
                     <<last_selection_point.index()<<std::endl;
            }
            if(pait->index() > last_selection_point.index()){//to avoid duplicates 
              if(print){
                std::cout<<NT(comb(NT(total_count),NT(k-D)))<<","
                         <<NT(comb(NT(N-total_count-D),NT(k-D)))<<std::endl;
              }
              //ccw combinations of points
              mean_vol += NT(comb(NT(total_count),NT(k-D)))*vol; 
              //cw combinations of points
              mean_vol += NT(comb(NT(N-total_count-D),NT(k-D)))*vol*NT(-1); 
            }
          }
      }
      if(print){
        std::cout<<"\nMean "<<k<<"-Volume="<<mean_vol<<std::endl;
        std::cout<<"\n================================\n";
      }
    }
    while(stdcomb::next_combination(ca.begin(),ca.end(),cb.begin(),cb.end()));
    
    //Correct the sign of the volume 
    mean_vol *= (D%2==0) ? 1 : -1; 
    
    return mean_vol;
}    

} // namespace MeanVol

#endif // MEAN_VOLUME_FUNCTIONS_H
