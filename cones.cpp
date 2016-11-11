/*
 * multiplane.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: bmetcalf & mpellejero
 *
 *      This program is for testing the adaptive griding on a random field by adapting to the
 *      high convergence regions.
 */

#include <slsimlib.h>
#include <standard.h>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
//#include <geometry.h>
#include "elliptic.h"
#include "gridmap.h"
#include "lightcone_construction.h"

using namespace std;

static std::mutex barrier;

int main(int arg,char **argv){
  
  COSMOLOGY cosmo(Planck1yr);
  Point_3d xo,v;
  
  std::vector<std::string> filenames;
  std::vector<double> redshifts(1,0.0);
  
  std::string dir = "/data1/users/gustavo/BigMD/2.5_3840_Planck1/ROCKSTAR/";
  
  //filenames.push_back(dir +"out_80p.list");  // snapshot
  redshifts.push_back(0.04603);                       // maximum redshift
  filenames.push_back(dir +"out_77p.list");
  //redshifts.push_back(0.1058);
  //filenames.push_back(dir +"out_74p.list");
  redshifts.push_back(0.1131);
  filenames.push_back(dir +"out_73p.list");
  redshifts.push_back(0.1636);
  filenames.push_back(dir +"out_64p.list");
  redshifts.push_back(0.2279);
  filenames.push_back(dir +"out_54p.list");
  //redshifts.push_back(0.2464);
  //filenames.push_back(dir +"out_52p.list");
  redshifts.push_back(0.2849);
  filenames.push_back(dir +"out_49p.list");
  redshifts.push_back(0.3256);
  //filenames.push_back(dir +"out_46p.list");
  //redshifts.push_back(0.3581);
  filenames.push_back(dir +"out_44p.list");
  redshifts.push_back(0.4156);
  filenames.push_back(dir +"out_40p.list");
  redshifts.push_back(0.4916);
  filenames.push_back(dir +"out_34p.list");
  //redshifts.push_back(0.5053);
  //filenames.push_back(dir +"out_33p.list");
  //redshifts.push_back(0.547);
  //filenames.push_back(dir +"out_30p.list");
  redshifts.push_back(0.5618);
  filenames.push_back(dir +"out_29p.list");
  redshifts.push_back(0.6383);
  filenames.push_back(dir +"out_24p.list");
  //redshifts.push_back(0.6714);
  //filenames.push_back(dir +"out_22p.list");
  redshifts.push_back(0.7053);
  filenames.push_back(dir +"out_20p.list");
  //redshifts.push_back(0.7232);
  //filenames.push_back(dir +"out_19p.list");
  redshifts.push_back(0.7976);
  filenames.push_back(dir +"out_15p.list");
  redshifts.push_back(0.8868);
  filenames.push_back(dir +"out_12p.list");
  redshifts.push_back(1.0);
  filenames.push_back(dir +"out_11p.list");
  redshifts.push_back(1.445);
  filenames.push_back(dir +"out_10p.list");
  redshifts.push_back(2.145);
  filenames.push_back(dir +"out_9p.list");
  redshifts.push_back(2.484);
  filenames.push_back(dir +"out_8p.list");
  redshifts.push_back(2.891);
  filenames.push_back(dir +"out_7p.list");
  redshifts.push_back(3.0);
  
  
  // shift the redshifts to between the snapshots
  for(int i=1;i<redshifts.size()-1;++i){
    redshifts[i] = 0.5*(redshifts[i] + redshifts[i+1]);
  }
  
  time_t to,t1;
  Utilities::RandomNumbers_NR ran(12709432);
  
  // random observer in the box
  double boxwidth = 3.688948e3;
  
  // random direction

  const int Ncones = 5;
  std::vector<Point_3d> observers(Ncones);
  std::vector<Point_3d> directions(Ncones);
  for(auto &xo : observers){
    xo[0] = ran()*boxwidth;
    xo[1] = ran()*boxwidth;
    xo[2] = ran()*boxwidth;
  }
  for(auto &v : directions){
    v[0] = ran.gauss();
    v[1] = ran.gauss();
    v[2] = ran.gauss();
    v /= v.length();
  }
  
  
  /**********************************************
  // *** for tests ***
  filenames.clear();
  filenames.push_back("file_example.dat");
  redshifts.resize(2);
  redshifts[1] = 2.5;
  // aim the first 4 cones at the first four halos
  Point_3d tmp;
  tmp[0]=3.54674; tmp[1] = 1.12347; tmp[2]=13.69941;
  tmp /= 0.677700;
  directions[0] = tmp - observers[0];
  std::cout << " r = " << directions[0].length() << std::endl;
  tmp[0]=0.92366; tmp[1] =  9.52555; tmp[2]= 16.51622;
  tmp /= 0.677700;
  directions[1] = tmp - observers[1];
  std::cout << " r = " << directions[1].length() << std::endl;
  tmp[0]=9.45172; tmp[1] =  3.61210; tmp[2]= 8.19601;
  tmp /= 0.677700;
  directions[2] = tmp - observers[2];
  std::cout << " r = " << directions[2].length() << std::endl;
  tmp[0]=3.42866; tmp[1] =  2.09676; tmp[2]= 13.78150;
  tmp /= 0.677700;
  directions[3] = tmp - observers[3];
  std::cout << " r = " << directions[3].length() << std::endl;
  // **********************************************/
  
  MultiLightCone mcone(degreesTOradians,observers,directions);
  std::vector<std::vector<LightCone::DataRockStar> > conehalos(Ncones);
  for(auto &c : conehalos) c.reserve(100000);
  
  time(&to);
  for(int i=0 ; i < filenames.size() ; ++i){
    
    std::cout << "Reading from catalog: " << filenames[i] << " z = " << redshifts[i] << " to "
    << redshifts[i+1] << std::endl;
    mcone.ReadBoxRockStar(filenames[i]
                         ,cosmo.coorDist(redshifts[i])
                         ,cosmo.coorDist(redshifts[i+1])
                         ,conehalos);
    
  }
  
  for(int i=0;i<Ncones;++i){
    std::cout << "Number of halos in the cones: " << conehalos[i].size() << std::endl;
    LightCone::WriteLightCone("cone_output_p" + std::to_string(i) + ".csv", conehalos[i]);
  }
  time(&t1);
  
  std::cout << "time for constructing cone : " << difftime(t1,to)/60 << " min"
  << std::endl;
  
}


