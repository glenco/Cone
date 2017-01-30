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
  
  // output directory
  const std::string outdir = "Output/";
  
  const COSMOLOGY cosmo(BigMultiDark);

  //// stuff about observers and fields
  const int Ncones = 10;
  const double range = 6*degreesTOradians;
  const double angular_resolution = range/512;
  std::vector<double> zsources = {2.297,1.075,0.4892};
  std::vector<std::vector<PixelMap> > maps;


  /// stuff about the input snapshots
  const double BoxLength = 2.5e3/cosmo.gethubble();
  const float particle_mass = 2.359e10/cosmo.gethubble()/0.005;   // 0.5% of particles are  used
  std::vector<std::string> snap_filenames;
  std::vector<float> snap_redshifts;
  
  
  // random observer in the box

  /*{
    std::string dir = "/data1/users/gustavo/BigMD/2.5_3840_Planck1/DM_SAMPLES/dm_particles_snap_0";
    std::string suffix = ".dat";
    
    std::vector<std::string> num = {"77","73","64","54","49","44","40","34","29","24","20","15","12","11","10","09","08","07"};
    
    for(int i = 0; i < num.size() ; ++i ) snap_filenames.push_back(dir + num[i] + suffix);
    
    snap_redshifts.push_back(0.04603);
    snap_redshifts.push_back(0.1131);
    snap_redshifts.push_back(0.1636);
    snap_redshifts.push_back(0.2279);
    snap_redshifts.push_back(0.2849);
    snap_redshifts.push_back(0.3256);
    snap_redshifts.push_back(0.4156);
    snap_redshifts.push_back(0.4916);
    snap_redshifts.push_back(0.5618);
    snap_redshifts.push_back(0.6383);
    snap_redshifts.push_back(0.7053);
    snap_redshifts.push_back(0.7976);
    snap_redshifts.push_back(0.8868);
    snap_redshifts.push_back(1.0);
    snap_redshifts.push_back(1.445);
    snap_redshifts.push_back(2.145);
    snap_redshifts.push_back(2.484);
    snap_redshifts.push_back(2.891);
    snap_redshifts.push_back(3.0);
  }*/
  /*{
    std::string dir = "Data/dm_particles_snap_0";
    std::string suffix = ".dat";
    
    std::vector<std::string> num = {"77"};
    
    for(int i = 0; i < num.size() ; ++i ) snap_filenames.push_back(dir + num[i] + suffix);
    
    snap_redshifts.push_back(0.04603);
  }*/
  {
    snap_filenames.push_back("Data/head.dat");
    snap_redshifts.push_back(0.04603);
  }
  
  time_t to,t1;
  time(&to);

  
  // random number generator
  Utilities::RandomNumbers_NR ran(12709432);
  
  // random direction and observers
  std::vector<Point_3d> observers(Ncones);
  std::vector<Point_3d> directions(Ncones);
  for(auto &x : observers){
    x[0] = ran()*BoxLength;
    x[1] = ran()*BoxLength;
    x[2] = ran()*BoxLength;
  }
  for(auto &v : directions){
    v[0] = ran.gauss();
    v[1] = ran.gauss();
    v[2] = ran.gauss();
    v /= v.length();
  }
  
  LightCones::FastLightCones(cosmo,zsources,maps,range
                                  ,angular_resolution
                                  ,observers
                                  ,directions
                                  ,snap_filenames
                                  ,snap_redshifts
                                  ,BoxLength
                                  ,particle_mass
                                  ,true);
/*  void FastLightCones(
                      const COSMOLOGY &cosmo
                      ,const std::vector<double> &zsources
                      ,std::vector<std::vector<PixelMap> > &maps
                      ,double range
                      ,double angular_resolution
                      ,std::vector<Point_3d> &observers
                      ,std::vector<Point_3d> &directions
                      ,const std::vector<std::string> &snap_filenames
                      ,const std::vector<float> &snap_redshifts
                      ,double BoxLength
                      ,double particle_mass
                      ,bool verbose
                      ,bool addtocone
                      ){
*/
  // output the maps
  for(int icone = 0 ; icone < Ncones/2 ; ++icone){
    for(int i=0 ; i < zsources.size() ; ++i){
      maps[icone][i].printFITS("!" + outdir + "kappa_c" + std::to_string(icone) + "z"
                               + std::to_string(zsources[i]) + ".fits" );
    }
  }
  // caclulate and output power spectra
  int N = 50;
  std::vector<PosType> pspectrum(N),multipole(N);
  for(int i=0 ; i < zsources.size() ; ++i){
    std::vector<PosType> totalpower(N,0);
    for(int icone = 0 ; icone < Ncones ; ++icone){
      maps[icone][i].PowerSpectrum(pspectrum,multipole);
      
      {
        std::ofstream ps_file(outdir + "kappa_c" + std::to_string(icone) + "z"
                              + std::to_string(zsources[i]) + "PS" + ".csv");
        ps_file << "l,PS" << endl;
        for(int ii=0;ii<pspectrum.size();++ii){
          ps_file << multipole[ii] << "," << pspectrum[ii] << endl;
        }
        ps_file.close();
      }
      
      for(int ii=0 ; ii<N ; ++ii ) totalpower[ii] += pspectrum[ii]/Ncones;
    }
    
    {
      std::ofstream ps_file(outdir + "kappaAve_z"
                            + std::to_string(zsources[i]) + "PS" + ".csv");
      ps_file << "l,PS" << endl;
      for(int ii=0;ii<N;++ii){
        ps_file << multipole[ii] << "," << totalpower[ii] << endl;
      }
      ps_file.close();
    }
  }
  
  
  time(&t1);
  
  std::cout << "time for constructing cone : " << difftime(t1,to)/60 << " min"
  << std::endl;
  
}


