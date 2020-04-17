/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "core/FlexibleBin.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include <iostream>
#include <limits>
#include <iterator>
#include <fstream>

// introducing openMP (OMP) parallelisation
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

using namespace std;

namespace PLMD
{
namespace bias
{

//+PLUMEDOC BIAS SITH
//+ENDPLUMEDOC

class SITH : public Bias
{
private:
  int sithstride;   // Stride to perform the cv clustering and generate a new biasing potential
  int cvstride;     // stride to print the value of the cvs in the cv file (for postprocessing)
  string sithfile;  // The name of the file that will contain the hyperbins found at each sithstride
  string cvfile;    // The name of the file that will contain the cvs involved in the taboo search
  double dist;      //Distance in CV space to generate a new hyperbin
  bool walkers_mpi; // multiple walkers same way as metaD

public:
  explicit SITH(const ActionOptions &);
  void calculate();
  static void registerKeywords(Keywords &keys);
};

PLUMED_REGISTER_ACTION(SITH, "SITH")

void SITH::registerKeywords(Keywords &keys)
{
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory", "SITHSTRIDE", "Frequency with which the snapshots are clustered and the bias is updated");
  keys.add("compulsory", "CVSTRIDE", "stride to print the value of the cvs in the cv file (for postprocessing)");
  keys.add("compulsory", "SITHFILE", "The name of the file that will contain the clusters found at each sithstride");
  keys.add("compulsory", "CVFILE", "The name of the file that will contain the cvs involved in the taboo search (read+write)");
  keys.add("compulsory", "DIST", "Distance in CV space to generate a new hyperbin");
  keys.addFlag("WALKERS_MPI", false, "Switch on MPI version (only version available) of multiple walkers");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias", "default", "the instantaneous value of the bias potential");
  keys.addOutputComponent("force2", "default", "the instantaneous value of the squared force due to this bias potential");
}

SITH::SITH(const ActionOptions &ao) : PLUMED_BIAS_INIT(ao),
                                      walkers_mpi(false)
{
  // Note sizes of these vectors are automatically checked by parseVector :-)
  parse("SITHSTRIDE", sithstride);
  parse("CVSTRIDE", cvstride);
  parse("CVFILE", cvfile);
  parse("SITHFILE", sithfile);
  parse("DIST", dist);
  parseFlag("WALKERS_MPI", walkers_mpi);
  checkRead();

  int rank;
  int nranks;

  if (walkers_mpi)
  {
    rank = multi_sim_comm.Get_rank();
    nranks = multi_sim_comm.Get_size();
    if (rank == 0)
    {
      printf("Multiple walkers active using MPI communication\n");
      cout << "Using " << nranks << " MPI ranks\n";
    }
    log.printf("Simulation corresponding to MPI rank %i.\n", rank);
  }
  else
  {
    rank = 0;
    nranks = 1;
    printf("No WALKERS_MPI\n");
  }

  //exit(0);

  addComponent("bias");
  componentIsNotPeriodic("bias");
  addComponent("force2");
  componentIsNotPeriodic("force2");

  if (rank == 0)
  {
    printf("Initialising SITH sampling protocol\n");
    printf("The CVs are going to be printed every %i steps.\n", cvstride);
    printf("Bias is going to be updated every %i steps.\n", sithstride);

    ofstream wfile;
    wfile.open(cvfile.c_str(), std::ios_base::app);
    wfile << "Rank Time ";
    for (int i = 0; i < getNumberOfArguments(); i++)
    {
      wfile << "CV_" << i << " "; //How do you get the CV labels?
    }
    wfile << endl;
    wfile.close();

    ofstream clustfile;
    clustfile.open(sithfile.c_str(), std::ios_base::app);
    clustfile << "Time_print Rank_clust Time_bin Population "; //Time_print is the time at which it has been printed, Time_bin is the time of the first snapshot in this bin
    for (int i = 0; i < getNumberOfArguments(); i++)
    {
      clustfile << "CV_" << i << " "; //How do you get the CV labels?
    }
    clustfile << endl;
    clustfile.close();
  }
}

vector<unsigned> bins; // Contains the population of each bin
vector<unsigned> bins_step; // Contains the step number of each bin center
vector<unsigned> bins_rank; // Contains the mpi rank of each bin center
unsigned n_bins;
vector<vector<double> > bin_values;


void SITH::calculate()
{
  int rank;
  int nranks;
  if (walkers_mpi)
  {
    rank = multi_sim_comm.Get_rank();
    nranks = multi_sim_comm.Get_size();
  }
  else
  {
    rank = 0;
    nranks = 1;
  }
  double dist2 = pow(dist, 2);
  int step = getStep();
  double time = getTime();

  const double pi = 3.1415926535897;
  double ene = 0.0;
  double totf2 = 0.0;
  double f = 0.;

  vector<double> cv;
  vector<double> cv_m;
  for (unsigned i = 0; i < getNumberOfArguments(); ++i)
  {
    cv.push_back(getArgument(i));
    cv_m.push_back(getArgument(i)); //this is the one that passed through MPI
  }
  vector<double> potentials(getNumberOfArguments(), 0);
  vector<double> forces(getNumberOfArguments(), 0);

  

  ////////////////////////////////////////////////
  //                                            //
  // UPDATE THE HYPERBINS AND THEIR POPULATIONS //
  //                                            //
  ////////////////////////////////////////////////
  int mod = step % cvstride;
  if (mod == 0)
  {
    // Print cvfile: Rank, Time, cv_1, cv_2...
    ofstream wfile;
    wfile.open(cvfile.c_str(), std::ios_base::app);
    wfile << rank << " " << time << " ";
    for (int i = 0; i < getNumberOfArguments(); i++)
    {
      double cv = getArgument(i);
      wfile << cv << " ";
    }
    wfile << endl;
    wfile.close();
    
    //send the value of the cvs to rank 0
    for (unsigned i=0; i<getNumberOfArguments();i++)
    {
    multi_sim_comm.Isend(cv_m[i],rank,0);
    }
    
    if (multi_sim_comm.Get_rank() == 0)
    {
      //At the first step, generate a bin with the first snapshot
      if (step == 0)
      {
        bins.push_back(nranks);
        bin_values.push_back(cv);
      }
      //Afterwards, add the snapshot to the relevant bin or generate a new bin
      else
      {
        //Assume the snapshot will make a new bin
        bool new_bin = true;
        for (unsigned j = 0; j < bins.size(); j++)
        {
          //Measure the distance to the bin
          double r2 = 0;
          for (unsigned i = 0; i < getNumberOfArguments(); i++)
          {
            r2 += pow((cv[i] - bin_values[j][i]), 2);
          }
          //if it is within bin boundaries, add 1 element to the current bin
          if (r2 <= dist2)
          {
            bins[j]++;
            new_bin = false;
          }
        }
        if (new_bin)
        {
          bins.push_back(1);
          bin_values.push_back(cv);
          n_bins = bins.size();
        }

        //Do the same for ranks 1 to m
        for (unsigned m = 1; m < nranks; m++)
        {
          //Assume the snapshot will make a new bin
          bool new_bin = true;
          for (unsigned j = 0; j < bins.size(); j++)
          {
            //Measure the distance to the bin
            double r2 = 0;
            for (unsigned i = 0; i < getNumberOfArguments(); i++)
            {
              multi_sim_comm.Recv(cv_m[i],m,0);
              r2 += pow((cv_m[i] - bin_values[j][i]), 2);
            }
            //if it is within bin boundaries, add 1 element to the current bin
            if (r2 <= dist2)
            {
              bins[j]++;
              new_bin = false;
            }
          }
          if (new_bin)
          {
            bins.push_back(1);
            bin_values.push_back(cv_m);
          }
        }
      }
    }
    //Bring the bins from rank 0
    multi_sim_comm.Bcast(n_bins, 0);
    for (unsigned j=0; j<bins.size();j++)
    {
      multi_sim_comm.Bcast(bins[j], 0);
      for (unsigned i=0; i<getNumberOfArguments();i++)
      {
       multi_sim_comm.Bcast(bin_values[j][i], 0);
      }
    }
    multi_sim_comm.Barrier(); //wait for everyone else
  }

  ////////////////////////////////////////////////
  //                                            //
  // GENERATE THE BIASING POTENTIAL AND FORCE   //
  //                                            //
  ////////////////////////////////////////////////

  //Get the distance between the current frame and all bins
  vector<double> r2_bins(bins.size(), 0);
  vector<double> r_bins(bins.size(), 0);
  vector<double> drbins(bins.size(), 0);
  for (unsigned j=0; j<n_bins;j++)
  {
    r2_bins[j]=0;
    r_bins[j]=0;
    for (unsigned i = 0; i < getNumberOfArguments(); i++)
    {
      r2_bins[j] += pow((cv[i] - bin_values[j][i]), 2);
    }
    r_bins[j]=sqrt(r2_bins[j]);
    drbins[j]=1/(2*r_bins[j]);
  }

  //Calculate the biasing potential and its derivative for each bin
  vector<double> V_r(bins.size(),0.);
  vector<double> dV_dr(bins.size(),0.);
  for (unsigned i = 0; i < getNumberOfArguments(); i++)
  {
    for (unsigned j=0; j<n_bins;j++)
    {
      if (r_bins[j]>dist)
          continue;
      double k_j=bins[j];
      V_r[j]=k_j*pow((r_bins[j]-dist),2);
      dV_dr[j]=2*k_j*(r_bins[j]-dist)*drbins[j];
    }
  }

  // Calculate the potential and the force with respect to each CV
  vector<vector<double> > potFor(2);
  vector<double> V_cv(getNumberOfArguments(),0.);
  vector<double> dV_cv(getNumberOfArguments(),0.);
  for (unsigned i=0; i<getNumberOfArguments();i++)
  {
   for (unsigned j=0; j<bins.size();j++)
   {
    
   }
  }


  
  for (unsigned i = 0; i < getNumberOfArguments(); ++i) // It only iterates one time, but I still don't know how to
  {
    ene += potentials[i];
    f = forces[i];
    totf2 += f * f;
    setOutputForce(i, f);
  }
  getPntrToComponent("bias")->set(ene);
  getPntrToComponent("force2")->set(totf2);
}
} // namespace bias
} // namespace PLMD
