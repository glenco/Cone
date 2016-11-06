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
  
  std::vector<LightCone::DataRockStar> conehalos;
  conehalos.reserve(100000);
  LightCone cone(degreesTOradians);
  
  std::vector<std::string> filenames;
  std::vector<double> redshifts(1,0.0);
  
  std::string dir = "../Catalogs/";
  
  filenames.push_back(dir +"out_80p_unweight.list");  // snapshot
  redshifts.push_back(0.04603);                       // maximum redshift
  filenames.push_back(dir +"out_77p_unweight.list");
  redshifts.push_back(0.1058);
  filenames.push_back(dir +"out_74p_unweight.list");
  redshifts.push_back(0.1131);
  filenames.push_back(dir +"out_73p_unweight.list");
  redshifts.push_back(0.1636);
  filenames.push_back(dir +"out_64p_unweight.list");
  redshifts.push_back(0.2279);
  filenames.push_back(dir +"out_54p_unweight.list");
  redshifts.push_back(0.2464);
  filenames.push_back(dir +"out_52p_unweight.list");
  redshifts.push_back(0.2849);
  filenames.push_back(dir +"out_49p_unweight.list");
  redshifts.push_back(0.3256);
  filenames.push_back(dir +"out_46p_unweight.list");
  redshifts.push_back(0.3581);
  filenames.push_back(dir +"out_44p_unweight.list");
  redshifts.push_back(0.4156);
  filenames.push_back(dir +"out_40p_unweight.list");
  redshifts.push_back(0.4916);
  filenames.push_back(dir +"out_34p_unweight.list");
  redshifts.push_back(0.5053);
  filenames.push_back(dir +"out_33p_unweight.list");
  redshifts.push_back(0.547);
  filenames.push_back(dir +"out_30p_unweight.list");
  redshifts.push_back(0.5618);
  filenames.push_back(dir +"out_29p_unweight.list");
  redshifts.push_back(0.6383);
  filenames.push_back(dir +"out_24p_unweight.list");
  redshifts.push_back(0.6714);
  filenames.push_back(dir +"out_22p_unweight.list");
  redshifts.push_back(0.7053);
  filenames.push_back(dir +"out_20p_unweight.list");
  redshifts.push_back(0.7232);
  filenames.push_back(dir +"out_19p_unweight.list");
  redshifts.push_back(0.7976);
  filenames.push_back(dir +"out_15p_unweight.list");
  redshifts.push_back(1.0);
  filenames.push_back(dir +"out_11p_unweight.list");
  redshifts.push_back(1.445);
  filenames.push_back(dir +"out_10p_unweight.list");
  redshifts.push_back(2.145);
  filenames.push_back(dir +"out_9p_unweight.list");
  redshifts.push_back(2.484);
  filenames.push_back(dir +"out_8p_unweight.list");
  redshifts.push_back(2.891);
  filenames.push_back(dir +"out_7p_unweight.list");
  redshifts.push_back(3.0);
  
  
  //filenames.clear();
  //filenames.push_back("file_example.dat");
  //redshifts.resize(2);
  //redshifts[1] = 0.1;
  
  // shift the redshifts to between the snapshots
  for(int i=1;i<redshifts.size()-1;++i){
    redshifts[i] = 0.5*(redshifts[i] + redshifts[i+1]);
  }
  
  time_t to,t1;
  Utilities::RandomNumbers_NR ran(-1928376);
  
  // random observer in the box
  double boxwidth = 1.0e3;
  xo[0] = ran()*boxwidth;
  xo[1] = ran()*boxwidth;
  xo[2] = ran()*boxwidth;
  
  // random direction
  v[0] = ran.gauss();
  v[1] = ran.gauss();
  v[2] = ran.gauss();
  
  
  time(&to);
  for(int i=0 ; i < filenames.size() ; ++i){
    
    cone.ReadBoxRockStar(filenames[i],xo,v
                         ,cosmo.coorDist(redshifts[i])
                         ,cosmo.coorDist(redshifts[i+1])
                         ,conehalos);
    
    std::cout << "Number of halos in the cone: " << conehalos.size() << std::endl;
  }
  
  LightCone::WriteLightCone("cone_output", conehalos);
  
  time(&t1);
  
  std::cout << "time for constructing cone : " << difftime(t1,to)/60 << " min"
  << std::endl;
  
}


