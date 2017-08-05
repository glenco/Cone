
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
#include <time.h>

//#include <geometry.h>
#include "elliptic.h"
#include "gridmap.h"
#include "lightcone_construction.h"
#include "profiles.h"

using namespace std;

static std::mutex barrier;

//**************** test ************************

struct dTNFW{
  dTNFW(double c){
    tnfw = new Profiles::TNFW2D(c);
  }
  ~dTNFW(){delete tnfw;}
  
  Profiles::TNFW2D* tnfw;
  double operator()(double x){return 2*pi*x*(*tnfw)(x);}
};

struct dNFW{
  dNFW(double my_c):cons(my_c){}
  
  double cons;
  
  Profiles::NFW2D nfw;
  double operator()(double x){return 2*pi*x*nfw(x,cons);}
};

struct dspline{
  double operator()(double q){return 2*pi*q*Profiles::Bspline<2>(q);}
};

//***********************************************/


int main(int arg,char **argv){
  
  /**************** test ************************

  dspline integrand;
  std::cout << Utilities::nintegrate<dspline,double>(integrand,0,5,1.0e-3)
  << std::endl;
  
  double c =10.0;
  
  //Profiles::NFW2D nfwprofile;
  //for(double x=0.00100591;x <c*1.01;x *= 1.01){
  //  cout << x << " " << nfwprofile(x,c) << endl;
  //}

  for(c=1.1;c<21;c *= 1.1){
    dNFW dnfw(c);
    std::cout << c << " " << Utilities::nintegrate<dNFW,double>(dnfw,0.0001,c,1.0e-5)
  << std::endl;
  }
  //cout.precision(17);
  //cout << 10./(7*pi)/4 << endl;
  
  exit(0);
  //***********************************************/

  
  Utilities::print_date();
  
  // output directory
  const std::string outdir = "Output/";
  // simulation directory
  //const std::string datadir = "/home/marcos/LensingCatalogs/LSSandHaloCatalogs/";
  const std::string datadir = "Data/";

  COSMOLOGY cosmo(BigMultiDark);
  
  //// stuff about observers and fields

  const int Ncones = 10;

  const double range = 6*degreesTOradians;
  const double angular_resolution = range/512/2;
  std::vector<double> zsources = {2.297,1.075,0.4892};
  std::vector<std::vector<PixelMap> > mapsLSS;
  std::vector<std::vector<PixelMap> > mapsHALO;


  /// stuff about the input snapshots
//  const double BoxLength = 2.5e3/cosmo.gethubble();
  //const float particle_mass = 2.359e10/cosmo.gethubble()/0.005;   // 0.5% of particles are  used
  const double BoxLength = 2500;
  const float particle_mass = 1.0;   // given in file
  std::vector<std::string> snap_filenamesHALO;
  std::vector<std::string> snap_filenamesLSS;
  std::vector<float> snap_redshifts;
  
  /*{  // from Gusavo's down sampled particle files
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
  
  /*{  // From halo and LSS partical catalogs
   
   //std::string dir = "LSS_catalog/dm_particles_snap_0";
   std::string dir = "/home/marcos/LensingCatalogs/LSSandHaloCatalogs";
   // Use LightCones::FastLightCones<LightCones::ASCII_XMRRT> with halo catalogs
   //std::string suffix = "KFth0.20.dat";

   // Use LightCones::FastLightCones<LightCones::ASCII_XMR> with LSS snapshots
   //std::string suffix = "KFth0.20_lss.dat";
   
   std::vector<std::string> num = {"77","73","64","54","49","44","40","34","29","24","20","15","12","11","10","09","08","07"};
   
    for(int i = 0; i < num.size() ; ++i ) snap_filenamesLSS.push_back(dir + num[i] + "KFth0.20_lss_noS.dat");
    for(int i = 0; i < num.size() ; ++i ) snap_filenamesHALO.push_back(dir + num[i] + "KFth0.20.dat");
   
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
   }/**/
   
   //********************* below are certain test cases without all the input files ******************
  /*{
    std::string dir = "Data/dm_particles_snap_0";
    std::string suffix = ".dat";
    
    std::vector<std::string> num = {"77","73"};
    //std::vector<std::string> num = {"77"};
    
    for(int i = 0; i < num.size() ; ++i ) snap_filenames.push_back(dir + num[i] + suffix);
    
    snap_redshifts.push_back(0.04603);
    snap_redshifts.push_back(0.1131);

  }
  }/**/
  /*{
    snap_filenames.push_back("Data/head.dat");
    snap_redshifts.push_back(0.04603);
  }*/
  

  /*{
    snap_filenamesHALO.push_back(datadir + "dm_particles_snap_007KFth0.20.dat");
    snap_filenamesLSS.push_back(datadir + "dm_particles_snap_007KFth0.20_lss.dat");
    snap_redshifts.push_back(0.04603);
    }*/
  {
    snap_filenamesHALO.push_back("Data/dm_particles_snap_007KFth0.20.dat");
    snap_filenamesLSS.push_back("Data/dm_particles_snap_007KFth0.20_lss.dat");
    snap_redshifts.push_back(0.04603);
  }
  /*{
    snap_filenames.push_back("Data/dm_particles_snap_007KFth0.20.head.dat");
    snap_redshifts.push_back(0.04603);
    }/**/
  //**********************************************************************************
  
  time_t to,t1,t2;
  time(&to);

  
  // random number generator
  Utilities::RandomNumbers_NR ran(12709432);
  
  // random direction and observers
  std::vector<Point_3d> observers(Ncones);
  std::vector<Point_3d> directions(Ncones);
  
  // generate random observers within the box
  LightCones::random_observers(observers,directions,Ncones,BoxLength
                               ,range/sqrt(2),ran);
    
  /*
   Template options are:
   LightCones::ASCII_XV   -- for 6 column ASCII file with position and velocity
   LightCones::ASCII_XM   -- for 4 column ASCII file with position and mass
   LightCones::ASCII_XMR   -- for 5 column ASCII file with position, mass and size
   LightCones::ASCII_XMRRT   -- for 7 column ASCII file with position, mass,Rmax,Rscale and an integer denoting type
   lightCone::ASCII_XMRRT12 -- same as LightCones::ASCII_XMRRT but only takes entries with the 7th column equal to 1 or 2

   and more to come.
   
   */
  int Npower = 50;
  std::vector<PosType> pspectrum(Npower,0),multipole(Npower,0);

  std::cout << "Starting LSS calculation..."  << std::endl;
  time(&t1);
  // This is for LSS particles
  LightCones::FastLightCones<LightCones::ASCII_XMR>(
                                                    cosmo,zsources,mapsLSS,range
                                                    ,angular_resolution
                                                    ,observers
                                                    ,directions
                                                    ,snap_filenamesLSS
                                                    ,snap_redshifts
                                                    ,BoxLength
                                                    ,particle_mass);/**/

  time(&t2);
  std::cout << "time for LSS cones: " << difftime(t2,t1)/60 << " min"
  << std::endl;
  std::cout << "Writing LSS files..."  << std::endl;

  // output the maps of first 5 cones
  for(int icone = 0 ; icone < min(5,Ncones) ; ++icone){
    for(int i=0 ; i < zsources.size() ; ++i){
      mapsLSS[icone][i].printFITS("!" + outdir + "kappa_c" + std::to_string(icone) + "z"
                                  + std::to_string(zsources[i]) + "LSS.fits" );
    }
  }
  
  
  // caclulate and output power spectra
  for(int i=0 ; i < zsources.size() ; ++i){
    std::vector<PosType> powerLSS(Npower,0);
    
    for(int icone = 0 ; icone < Ncones ; ++icone){
      mapsLSS[icone][i].PowerSpectrum(pspectrum,multipole);
      for(int ii=0 ; ii<Npower ; ++ii ) powerLSS[ii] += pspectrum[ii]/Ncones;
    }
    
    std::ofstream ps_file(outdir + "kappaAve_z"
                          + std::to_string(zsources[i]) + "PowLSS.csv");
    ps_file << "l,PS" << endl;
    for(int ii=0;ii<Npower;++ii){
      ps_file << multipole[ii] << "," << powerLSS[ii] << endl;
    }
    ps_file.close();
    
  }

  std::cout << "Starting halo calculation..."  << std::endl;
  time(&t1);
  // This is for halos.  You can also use ightCones::ASCII_XMRRT12
  LightCones::FastLightCones<LightCones::ASCII_XMRRT>(
                                                      cosmo,zsources,mapsHALO,range
                                                      ,angular_resolution
                                                      ,observers
                                                      ,directions
                                                      ,snap_filenamesHALO
                                                      ,snap_redshifts
                                                      ,BoxLength
                                                      ,particle_mass);
  
  /*LightCones::FastLightCones<LightCones::HDF5_XMRRT12>(
                                                      cosmo,zsources,mapsHALO,range
                                                      ,angular_resolution
                                                      ,observers
                                                      ,directions
                                                      ,snap_filenamesHALO
                                                      ,snap_redshifts
                                                      ,BoxLength
                                                      ,particle_mass);
   */

  time(&t2);
  std::cout << "time for Halo cones: " << difftime(t2,t1)/60 << " min"
  << std::endl;
  
  
  // output the maps of first 5 cones
  for(int icone = 0 ; icone < min(5,Ncones) ; ++icone){
    for(int i=0 ; i < zsources.size() ; ++i){
      mapsHALO[icone][i].printFITS("!" + outdir + "kappa_c" + std::to_string(icone) + "z"
                                   + std::to_string(zsources[i]) + "HALO.fits" );
    }
  }
  
  std::cout << "Writing HALO files..."  << std::endl;
 // caclulate and output power spectra
  for(int i=0 ; i < zsources.size() ; ++i){
    std::vector<PosType> powerHALO(Npower,0);
    
    for(int icone = 0 ; icone < Ncones ; ++icone){

      mapsHALO[icone][i].PowerSpectrum(pspectrum,multipole,1,true,true);
      for(int ii=0 ; ii<Npower ; ++ii ) powerHALO[ii] += pspectrum[ii]/Ncones;
    }
    
    {
      std::ofstream ps_file(outdir + "kappaAve_z"
                   + std::to_string(zsources[i]) + "PowHALO.csv");
      ps_file << "l,PS" << endl;
      for(int ii=0;ii<Npower;++ii){
        ps_file << multipole[ii] << "," << powerHALO[ii] << endl;
      }
      ps_file.close();
    }
    
    
  }
  
  // add maps for total, very inofficient
  std::vector<std::vector<PixelMap> > mapsTOTAL(Ncones);
  
  for(int icone = 0 ; icone < Ncones ; ++icone){
    for(int i=0 ; i < zsources.size() ; ++i){
      mapsTOTAL[icone].push_back(mapsLSS[icone][i]);
      mapsTOTAL[icone].back() += mapsHALO[icone][i];
    }
  }
  
  
  
  // output the maps of first 5 cones
  for(int icone = 0 ; icone < min(5,Ncones) ; ++icone){
    for(int i=0 ; i < zsources.size() ; ++i){
     mapsTOTAL[icone][i].printFITS("!" + outdir + "kappa_c" + std::to_string(icone) + "z"
                                    + std::to_string(zsources[i]) + "TOTAL.fits" );
    }
  }
  
  std::cout << "Writing TOTAL files..."  << std::endl;
  // caclulate and output power spectra
  for(int i=0 ; i < zsources.size() ; ++i){
    std::vector<PosType> powerTOTAL(Npower,0);
    
    for(int icone = 0 ; icone < Ncones ; ++icone){
      
       mapsTOTAL[icone][i].PowerSpectrum(pspectrum,multipole);
      for(int ii=0 ; ii<Npower ; ++ii ) powerTOTAL[ii] += pspectrum[ii]/Ncones;
      
    }
    
    {
      std::ofstream ps_file(outdir + "kappaAve_z"
                   + std::to_string(zsources[i]) + "PowTOTAL.csv");
      ps_file << "l,PS" << endl;
      for(int ii=0;ii<Npower;++ii){
        ps_file << multipole[ii] << "," << powerTOTAL[ii] << endl;
      }
      ps_file.close();
    }
  }
  time(&t1);
  
  std::cout << "time for constructing cone : " << difftime(t1,to)/60 << " min"
  << std::endl;
  
}


