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
  
  { // testing LightCone::ReadLightConeParticles
    COSMOLOGY cosmo(Planck1yr);

    std::vector<LensHaloParticles *> halovec;
    
    LightCone::ReadLightConeParticles("Cones/cone_particles0.csv", cosmo, halovec,10,1.0e10,0.0,false,true);
    
    long seed =-1927;
    Lens lens(&seed,3,cosmo);
    
    for(auto ptr_halo : halovec){
      lens.insertMainHalo(ptr_halo,true);
    }
    
    std::cout << "Number of planes : " << lens.getNplanes() << std::endl;
    exit(0);
  }
  
  COSMOLOGY cosmo(Planck1yr);
  Point_3d xo,v;
  
  std::vector<std::string> filenames;
  std::vector<double> redshifts(1,0.0);
  // random observer in the box
  double boxwidth = 3.688948e3,hubble = 0.677700;

  //std::string dir = "/data1/users/gustavo/BigMD/2.5_3840_Planck1/ROCKSTAR/out_";
  //std::string suffix = "p.list";
  std::string dir = "/data1/users/gustavo/BigMD/2.5_3840_Planck1/DM_SAMPLES/dm_particles_snap_0";
  std::string suffix = ".dat";
  
  std::vector<std::string> num = {"77","73","64","54","49","44","40","34","29","24","20","15","12","11","10","09","08","07"};

  for(int i = 0; i < num.size() ; ++i ) filenames.push_back(dir + num[i] + suffix);
  
  redshifts.push_back(0.04603);                       // maximum redshift
  redshifts.push_back(0.1131);
  redshifts.push_back(0.1636);
  redshifts.push_back(0.2279);
  redshifts.push_back(0.2849);
  redshifts.push_back(0.3256);
  redshifts.push_back(0.4156);
  redshifts.push_back(0.4916);
  redshifts.push_back(0.5618);
  redshifts.push_back(0.6383);
  redshifts.push_back(0.7053);
  redshifts.push_back(0.7976);
  redshifts.push_back(0.8868);
  redshifts.push_back(1.0);
  redshifts.push_back(1.445);
  redshifts.push_back(2.145);
  redshifts.push_back(2.484);
  redshifts.push_back(2.891);
  redshifts.push_back(3.0);
  
  
  // shift the redshifts to between the snapshots
  for(int i=1;i<redshifts.size()-1;++i){
    redshifts[i] = 0.5*(redshifts[i] + redshifts[i+1]);
  }
  
  time_t to,t1;
  Utilities::RandomNumbers_NR ran(12709432);
  
  
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
  
  MultiLightCone hcone(degreesTOradians,observers,directions);
  //std::vector<std::vector<LightCone::DataRockStar> > conehalos(Ncones);
  std::vector<std::vector<Point_3d> > cone_particles(Ncones);
  for(auto &c : cone_particles) c.reserve(100000);
  
  time(&to);
  for(int i=0 ; i < filenames.size() ; ++i){
    
    std::cout << "Reading from catalog: " << filenames[i] << " z = " << redshifts[i] << " to "
    << redshifts[i+1] << std::endl;
    //hcone.ReadBoxRockStar(filenames[i]
    //                     ,cosmo.coorDist(redshifts[i])
    //                     ,cosmo.coorDist(redshifts[i+1])
    //                     ,conehalos);
    
    hcone.ReadBoxXYZ(filenames[i]
                     ,cosmo.coorDist(redshifts[i])
                     ,cosmo.coorDist(redshifts[i+1])
                     ,cone_particles
                     ,hubble,boxwidth);
  }
  
  //for(int i=0;i<Ncones;++i){
  //  std::cout << "Number of halos in the cones: " << conehalos[i].size() << std::endl;
  //  LightCone::WriteLightCone("cone_output_p" + std::to_string(i) + ".csv", conehalos[i]);
  //}
  for(int i=0;i<Ncones;++i){
    std::cout << "Number of halos in the cones: " << cone_particles[i].size() << std::endl;
    LightCone::WriteLightCone("cone_output_p" + std::to_string(i) + ".csv", cone_particles[i]);
  }
  time(&t1);
  
  std::cout << "time for constructing cone : " << difftime(t1,to)/60 << " min"
  << std::endl;
  
}


