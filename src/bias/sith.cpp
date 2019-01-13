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

#include "newtypes.h"

// introducing openMP (OMP) parallelisation
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS SITH
/*
 Applies a potential with the shape of an inverted gaussian:
 
 V=(1/(SIGMA*sqrt(2*pi)))*(1-exp(-((CV-AT)**2)/(2*(SIGMA**2))))
 
 This kind of potential is useful when the CV fluctuates a lot and harmonic
 potentials give forces too large for the integrator to handle
  
 USAGE:
 restr: SITH ARG=je AT=10.0 SIGMA=2.0 
 
*/
//+ENDPLUMEDOC
    
class SITH : public Bias{
private:
  
  double height; // Factor that will rescale he cluster populations (for bias/force generation))
  int sithstride; // Stride to perform the cv clustering and generate a new biasing potential
  double sithstepsup; // number of steps over which the SITH bias is switched on
  double dc; // dc value for clustering (See Laio2014) )
  int dc_opt; //Frequency with which dc will be optimised
  double delta0; // delta0 value for clustering (see laio2014)(if DC_OPT is specified, this is an initial guess)
  int cvstride; // stride to print the value of the cvs in the file that will be read for clustering
  string sithfile; // The name of the file that will contain the clusters found at each sithstride (only write)
  string cvfile; // The name of the file that will contain the cvs involved in the taboo search (read+write)
  string typot; // The type of potential we will use to bias
  bool walkers_mpi; // multiple walkers same way as metaD
  bool limit_r; // limit the bias to a distance equal to delta0
  int iteration; // SITH iteration you are running
  
public:
  explicit SITH(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(SITH,"SITH")

void SITH::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SITHSTRIDE","Frequency with which the snapshots are clustered and the bias is updated");
  keys.add("optional","SITHSTEPSUP","number of steps over which the SITH bias is switched on");
  keys.add("compulsory","CVSTRIDE","Frequency with which the snapshots are clustered and the bias is updated");
  keys.add("optional","HEIGHT","Factor that scales the height of the gaussians");
  keys.add("compulsory","DC","dc value for the clustering (if DC_OPT is specified, this is an initial guess)");
  keys.add("optional","DCOPT","Frequency with which dc will be optimised");
  keys.add("compulsory","DELTA0","delta0 value for the clustering");
  keys.add("compulsory","SITHFILE","The name of the file that will contain the clusters found at each sithstride (only write)");
  keys.add("compulsory","CVFILE","The name of the file that will contain the cvs involved in the taboo search (read+write)");
  keys.add("compulsory","TYPOT","The type of potential we will use to bias");
  keys.add("optional","ITERATION","Taboo search iteration you are running");
  keys.addFlag("LIMIT_R",false,"Limit the bias to a distance equal to delta0");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version (only version available) of multiple walkers");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}


SITH::SITH(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
//height(getNumberOfArguments(),0),
//sithstride(getNumberOfArguments(),0)
walkers_mpi(false),
limit_r(false)
{
  // Setting defaults for optional keywords
  iteration=0;
  sithstepsup=0;
  height=1;
  dc_opt=0;  
  // Note sizes of these vectors are automatically checked by parseVector :-)
  parse("HEIGHT",height);
  parse("SITHSTRIDE",sithstride);
  parse("SITHSTEPSUP",sithstepsup);
  parse("CVSTRIDE",cvstride);
  parse("DC",dc);
  parse("DCOPT",dc_opt);
  parse("DELTA0",delta0);
  parse("CVFILE",cvfile);
  parse("SITHFILE",sithfile);
  parse("TYPOT",typot);
  parseFlag("WALKERS_MPI",walkers_mpi);
  parseFlag("LIMIT_R",limit_r);
  parse("ITERATION",iteration);
  checkRead();
 
  int rank;
  int nranks;
  
  if(walkers_mpi) 
  {
      rank=multi_sim_comm.Get_rank();
      nranks=multi_sim_comm.Get_size();
      if (rank==0)
      {
        printf("Multiple walkers active using MPI communication\n");
        cout << "Using " << nranks << " MPI ranks\n";
      }
      log.printf("Simulation corresponding to MPI rank %i.\n",rank);
  }
  else
  {
      rank=0;
      nranks=1;
      printf("No WALKERS_MPI\n");
  }
  
  //exit(0);
  
  addComponent("bias"); componentIsNotPeriodic("bias");
  addComponent("force2"); componentIsNotPeriodic("force2");
  
  if (rank==0)
  {
   printf("Initialising SITH sampling protocol at iteration %i\n", iteration);
   printf("The CVs are going to be printed every %i steps.\n", cvstride);
   printf("Clustering is going to be performed every %i steps.\n", sithstride);
   printf("with dc equal to %f and delta0 equal to %f.\n",dc,delta0);
   if (dc_opt!=0)
   {
       printf("The function to optimise dc needs to be debugged so it can't be used yet.\n");
       printf("Please delete keyword DCOPT or set it to 0 and restart the calculation\n");
       printf("EXITING\n");
       exit(0);
       //printf("dc will be optimised every %i steps.\n",dc_opt);
   }
   if (height!=1)
       printf("The population of the clusters will be rescaled by a factor of %f.\n", height);
   if (sithstepsup!=0) 
       printf("After clustering, the bias will be linearly switched on along %f steps.\n", sithstepsup);
   printf("Please read and cite: Rodriguez, A.; Laio, A.; Science (2014) 344(6191) p.1496\n");   
      
   if (iteration==0)
   {
     ofstream wfile;
     wfile.open(cvfile.c_str(),std::ios_base::app);
     wfile << "Iteration Rank Time ";
     for (int i=0; i<getNumberOfArguments();i++)
     {
      wfile << "CV_" << i << " "; //How do you get the CV labels?
     }
     wfile << endl;
     wfile.close();
  
     ofstream clustfile;
     clustfile.open(sithfile.c_str(),std::ios_base::app);
     clustfile << "Iteration Time_print Rank_clust Time_clust Population "; //Time_print is the time at which it has been printed, Time_clust is the time of the cluster center
     for (int i=0; i<getNumberOfArguments();i++)
     {
      clustfile << "CV_" << i << " "; //How do you get the CV labels?
     }
     clustfile << endl;
     clustfile.close();
   }
   
  }
}

/*
  //Data struct we will use to perform the clustering
 struct Laio
  {
    int snapshot;
    double rho;
    int cluster;
    int nnhd;
    double delta;
  };
  
    bool mycmp_rho(const Laio& a, const Laio& b)
  {
      return a.rho > b.rho;
  };

  
  //Data struct we will use to store the time, cv values, population and width of each cluster
  struct values {
  int rank;
  double time;
  vector<double> cvs;
  int  population;
  vector<double> sigma;
  values(int rank, double time,vector<double> cvs, int population, vector<double> sigma) 
    {
      this -> rank = rank;
      this -> time = time;
      this -> cvs = cvs;
      this -> population = population;
      this -> sigma = sigma;
    }
  };
*/

vector<values> getCVs(string CV_file)
{
 vector<values> values_raw;
 FILE* fp=fopen(CV_file.c_str(),"r");
 string line;
 while (Tools::getline(fp, line))
 {
  if (line[0] == 'R') continue;
  istringstream iss(line);
  istream_iterator<string> beg(iss), end;
  vector<string> tokens(beg, end);
  int iteration=atoi(tokens[0].c_str());
  int rank=atoi(tokens[1].c_str());
  double time=atof(tokens[2].c_str());
  vector<double> cvs;
  for (int i=3;i<tokens.size();i++)
  {
      cvs.push_back(atof(tokens[i].c_str()));
  }
  int population=0;
  vector<double> sigma;
  values_raw.push_back(values(iteration,rank,time, cvs, population, sigma));
 }
 //exit(0);
 return values_raw;
}

double optimise_dc(vector<values> & values_raw, double dc)
{
 int rho_avg=0;
 int rho_inst=0;
 double k=1.;
 while (rho_avg<0.01*values_raw.size() or rho_avg>0.02*values_raw.size())
 {
    double dc2=dc*dc;
    double r_avg=0;
    #pragma omp parallel for reduction(+:rho_avg, r_avg)
    for (int i=0; i<values_raw.size();i++)
    {
      for (int j=0;j<values_raw.size();j++)
          {
            if (j==i) continue;
            double r2=0;
            for (int k=0; k<values_raw[i].cvs.size();k++) 
                r2=+pow((values_raw[j].cvs[k]-values_raw[i].cvs[k]),2);
            r_avg += sqrt(r2);
            if (r2<dc2) rho_avg+=1;
          }
    }
    rho_avg /= values_raw.size();
    r_avg /= ((pow(values_raw.size(),2)-values_raw.size())/2);
    //cout << "Average distance is " << r_avg << endl;
    
    if (((rho_inst<0.01*values_raw.size()) and (rho_avg>0.02*values_raw.size()))
            or
        ((rho_inst>0.02*values_raw.size()) and (rho_avg<0.01*values_raw.size())))
    {
             cout << "decreasing k\n";
             k=k/10;
    }       
    
    cout << "Average density with dc = " << dc << ": " << rho_avg << endl;
    
    if (rho_avg<0.01*values_raw.size()) // Since we are comparing r2 to dc2, we need to check that sc didn't go below 0
       {
        rho_inst=rho_avg;
        dc += k * r_avg;
       }
    else if (rho_avg>0.02*values_raw.size())
       {
        rho_inst=rho_avg;
        dc -= k * r_avg;
        if (dc<0)
        {
         dc=0;
         k=k/10;     
        }
       }
 }
 cout << "Average density: " << rho_avg << endl;
 cout << "dc will be " << dc << endl;
 return dc;
}

vector<values> cluster_snapshots(vector<values> & values_raw, double dc, double delta0)
{
    double dc2=dc*dc;
    vector<Laio> vec(values_raw.size());
    
    //initialising vec and calculating rho
    
    #pragma omp parallel for
    for (unsigned i=0;i<values_raw.size();i++)
    {
        vec[i].snapshot=i;
        vec[i].cluster=-1; // Initialised at -1
        vec[i].rho=0;
        vec[i].nnhd=-1;
        vec[i].delta=0;
        for (unsigned j=0; j<values_raw.size();j++)
        {
         if(j==i) continue;
         double r2=0;
         for (unsigned k=0; k<values_raw[i].cvs.size();k++)
             r2 += pow((values_raw[i].cvs[k]-values_raw[j].cvs[k]),2);
         if (r2<dc2) vec[i].rho++;
        }
    }
    
    // Sorting vec by density
    sort(vec.begin(),vec.end(),mycmp_rho);

     
    
    // Findig nearest neigbour of higher density and calculating delta

    #pragma omp parallel for
    for (unsigned i=1;i<vec.size();i++)
    {
        int snap_i=vec[i].snapshot;
        double r2min=999999999999.;
        for (unsigned j=0; j<i; j++) // Points j have a higher density than points i
        {
            int snap_j=vec[j].snapshot;
            double r2=0;
            for (unsigned k=0; k<values_raw[i].cvs.size();k++)
                r2 += pow((values_raw[snap_i].cvs[k]-values_raw[snap_j].cvs[k]),2);
            if (r2<r2min)
            {
               vec[i].nnhd=snap_j;
               vec[i].delta=sqrt(r2);
               r2min=r2;
            }
        }
    }
    
    // For the point with the highest density, delta is the furthest point in the dataset
    int snap0=vec[0].snapshot;
        double r2max=0.;
        for (unsigned j=1; j<vec.size(); j++) 
        {
            int snap_j=vec[j].snapshot;
            double r2=0;
            for (unsigned k=0; k<values_raw[0].cvs.size();k++)
                r2 += pow((values_raw[snap0].cvs[k]-values_raw[snap_j].cvs[k]),2);
            if (r2>r2max)
            {
               vec[0].delta=sqrt(r2);
               r2max=r2;
            }
        }
    
    //Choosing cluster centres
    vector<values> clusters_raw;
    int nclust=0;
    for (unsigned i=0;i<vec.size();i++) //DON'T OMP THIS LOOP!! WILL MESS WITH THE CLUSTER INDICES
    {
        if (vec[i].delta>delta0)
        {
         int snap_i=vec[i].snapshot;
         clusters_raw.push_back(values_raw[snap_i]);
         clusters_raw[nclust].population=1;
         vec[i].cluster=nclust;
         nclust += 1;
        }
    }
    
    //Extreme case in which delta0 is too big
    if (nclust==0)
    {
     int snap_i=vec[0].snapshot;
     clusters_raw.push_back(values_raw[snap_i]);
     clusters_raw[nclust].population=1;
     vec[0].cluster=nclust;
     nclust += 1;
    }
    
    // Assigning points to clusters
    for (unsigned i=0;i<vec.size();i++) // DON'T OMP THIS LOOP. YOU WILL END UP LOOKING AT POINTS WHOSE NNHD HASN'T BEEN ASSIGNED YET
    {
        if (vec[i].cluster!=-1) // if cluster index is not -1 it means it's a centre assigned in the previous bit
            continue;
        int snap_i=vec[i].snapshot;
        //YOU COULD IN THEORY OMP THIS LOOP, BUT YOU WOULD BE CREATING AND DELEATING THREADS VERY OFTEN AND THE OVERHEAD WOULD MAKE IT SLOWER
        for (unsigned j=0; j<i; j++) //points j have higher density than points i;
        {
            int snap_j=vec[j].snapshot;
            if (vec[i].nnhd!=snap_j)
                continue;
            vec[i].cluster=vec[j].cluster;
            clusters_raw[vec[i].cluster].population += 1;
        }
    }
    
    // Finding the furthest point form the centre within each cluster
    for (unsigned i=0; i<clusters_raw.size();i++) // THE NUMBER OF CLUSTERS SHOULD IN PRINCIPLE BE SMALL SO NO POINT IN OMPing THIS LOOP
    {
       double r2max=0.;
       for (unsigned j=0; j<vec.size();j++) //DON'T OMP THIS LOOP!!! WILL MESS UP r2max
       {
           if(vec[j].cluster!=i)
               continue;
           int snap_j=vec[j].snapshot;
           double r2=0;
           for (unsigned k=0; k<values_raw[snap_j].cvs.size();k++)
           {
               r2 += pow((clusters_raw[i].cvs[k]-values_raw[snap_j].cvs[k]),2);
           }
           if (r2>=r2max)
           {
             clusters_raw[i].sigma= values_raw[snap_j].cvs;
             r2max=r2;
           }
       }
       //cout << "Sigma for cluster " << i << " equal to " << sqrt(r2max) << endl;
    }
    
   ofstream rhodelta;
   rhodelta.open("rhodelta.txt");
   rhodelta << "Index Snapshot Cluster Rho NNHD Delta\n";
   for (unsigned i=0; i<vec.size(); i++)    
    {
     rhodelta << i << " " << vec[i].snapshot << " " << vec[i].cluster << " " << vec[i].rho << " " << vec[i].nnhd << " " << vec[i].delta << endl;
    }
    cout << "Found " << nclust << " clusters" << endl;
        
    
    
  //exit(0);  
  return clusters_raw;
}

// JCN Jan2019: I think this was done for MPI related reasons.
vector<values> resize_clusters(int nClusters, int nArgs, int iteration)
{
    vector<values> clusters;
    int rank=-1;
    double time=-1;
    vector<double> cvs(nArgs,0);
    int  population=-1;
    vector<double> sigma(nArgs,0);
    for (unsigned i=0;i<nClusters;i++)
      clusters.push_back(values(iteration,rank,time,cvs,population,sigma));  
    
    return clusters;
}

vector<vector<double> > GenPot(string typot, vector<values> & clusters, double resc, double height,vector<double> & cv, double delta0, bool limit_r)
{
  // Calculate the current and equilibrium distances between the current data point and the cluster centres
  // 1) This gives us r=sqrt[(cv1-cv1_0)**2+...+(cvN-cvN_0)**2] and dr=1/2r
  vector<double> r2(clusters.size(),0);
  vector<double> at2(clusters.size(),0);
  //#pragma omp parallel for shared(r2)
  for (unsigned j=0; j<cv.size();j++)
  {
   for (unsigned i=0; i<clusters.size();i++)
   {
    r2[i] += pow((cv[j]-clusters[i].cvs[j]),2);
    at2[i] += pow((clusters[i].cvs[j]-clusters[i].sigma[j]),2);
   }
  }
  
  vector<double> r(clusters.size(),0);
  vector<double> at(clusters.size(),0); 
  vector<double> dr(clusters.size(),0); // Derivative of r it's always 1/(2*r)
  for (unsigned i=0; i<clusters.size();i++)
  {
    r[i]=sqrt(r2[i]);
    at[i]=sqrt(at2[i])*resc;
    if (limit_r)
       if (at[i]>delta0) at[i]=delta0*resc; 
    dr[i]=1/(2*r[i]);
    //cout << "r[ " << i << "] = " << r[i] << " at" << i << "] = " << at[i] << endl;
  }
  
  // Get V(r) and dV(r)
  // 2) This gives us V(r) and dV(r)=dV*dr
  vector<double> V_r(clusters.size(),0.);
  vector<double> dV_dr(clusters.size(),0.);
  if (typot=="LOWER_WALLS")
  {
    for (unsigned i=0; i<clusters.size();i++)
    {
     if (r[i] >= at[i]) continue;
     double k_i=clusters[i].population*height;
     V_r[i]=k_i*pow((r[i]-at[i]),2);
     dV_dr[i]=2*k_i*(r[i]-at[i])*dr[i];
     //cout << "V_r["<< i << "] is equal to " << V_r[i] << endl;
    }
  }
  else
  {
      cout << "This kind of potential is not implemented. Aborting." << endl;
      exit(0);
  }
  
  // Add the derivatives with respect to every CV
  // dV(r)=(dV(r)/dr)*(dr/dCV_i)*dCV_i
  vector<vector<double> > potFor(2);
  vector<double> V_cv(cv.size(),0.);
  vector<double> dV_cv(cv.size(),0.);
  //#pragma omp parallel for shared(V_cv,dV_cv)
  for (unsigned j=0; j<cv.size();j++)
  {
   for (unsigned i=0; i<clusters.size();i++)
   {
    V_cv[j] += V_r[i]/cv.size(); // This might be violating a few laws of physics
    double cv_j=clusters[i].cvs[j];
    double cv0_j=clusters[i].sigma[j];
    dV_cv[j] += dV_dr[i]*2*(cv_j-cv0_j);
   }
   
   potFor[0].push_back(V_cv[j]);
   potFor[1].push_back(-dV_cv[j]);
  }
  
  /*
  for (unsigned j=0; j<cv.size();j++)
  {
      cout << "biasing potential for CV " << j << " is " << potFor[j][0] << endl;
      cout << "negative derivative of the biasing potential for CV " << j << " is " << potFor[j][1] << endl;
  }
   */
   
  
  return potFor;
  
  //exit(0);
}
 
vector<values> values_raw;
vector<values> clusters;
 
void SITH::calculate(){
  // All this has to go here to initialise the bias at 0 for the first few steps
  int rank;
  int nranks;
  if (walkers_mpi)
  {
   rank=multi_sim_comm.Get_rank(); 
   nranks=multi_sim_comm.Get_size();
  }
  else
  {
   rank=0;
   nranks=1;
  }
  
  
  int step=getStep();
  double time=getTime();
  
  const double pi=3.1415926535897;
  double ene = 0.0;
  double totf2 = 0.0;
  double f=0.;
  
  vector<double> cv;
  for(unsigned i=0;i<getNumberOfArguments();++i) cv.push_back(getArgument(i));
  vector<double> potentials(getNumberOfArguments(),0);
  vector<double> forces(getNumberOfArguments(),0);
  
  
  // This is to transfer the cluster across all MPI ranks
  int nClusters=0;
  
  //Check if values have to be printed in cvfile
  int mod = step % cvstride;
  if (mod==0)
  {
   // cvfile has: Rank, Time, cv_1, cv_2...
   ofstream wfile;
   wfile.open(cvfile.c_str(),std::ios_base::app);
   wfile << iteration << " " << rank << " " << time << " ";
   for (int i=0; i<getNumberOfArguments();i++)
   {
       double cv=getArgument(i);
       wfile << cv << " ";
   }
   wfile << endl;
   wfile.close();
  }
  
  
  // Clustering snapshots and generating new potentials if necessary
  //FIXME:
  /*
   if (iteration>0 and step==0)
   * {
   *  read_the_bias_from_somewhere() // You shouldn't need to recluster
   * }
   * else if (step>=sithstride)
   * {
   *  run_clustering_as_below()
   * }
   */
  if ((step>=sithstride) or (iteration>0 and step==0))
  {   
     int mod = step % sithstride; // ... check if new clusters need to be calculated ...
     if (mod==0)                  // ... and do it if you have to ...
     {
       multi_sim_comm.Barrier();
       
       if (multi_sim_comm.Get_rank()==0) // ... go to rank 0 ...
       {
           
         cout << "generating new SITH potentials at time = " << time << " ps (assuming you are doing stuff in ps)." << endl;
         //cout << "Reading the values in CV file: ";
         values_raw=getCVs(cvfile);
         // Optimising dc
         if ((dc_opt!=0) and ((step==sithstride) or (step%dc_opt)==0)) 
             dc=optimise_dc(values_raw, dc);
         clusters=cluster_snapshots(values_raw,dc,delta0);
         nClusters=clusters.size();
         
         ofstream clustfile;
         clustfile.open(sithfile.c_str(),std::ios_base::app);
         clustfile << "---------------------------------------------------------" << endl;
         for (int i=0; i<clusters.size();i++)
          {
           clustfile << iteration << " " << time << " " << clusters[i].rank << " " << clusters[i].time << " " << clusters[i].population << " ";
           for (unsigned j=0; j<getNumberOfArguments(); j++)
            {
              clustfile << clusters[i].cvs[j] << " ";
            }
           clustfile << endl;
          }
         clustfile.close();
       }
       
       multi_sim_comm.Barrier();
         
       // Get the number of clusters
       multi_sim_comm.Bcast(nClusters,0);
       
       // Resize the clusters vector in ranks other than 0
       if (rank!=0)
           clusters=resize_clusters(nClusters,getNumberOfArguments(),iteration);  
       multi_sim_comm.Barrier();
       
       // Bring the clusters from rank 0 to other ranks
       for (unsigned i=0;i<nClusters;i++)
       {
        //cout << " Bringing time, rank and population of cluster " << i << " to rank " << rank << endl;
        multi_sim_comm.Bcast(clusters[i].iteration,0);   
        multi_sim_comm.Bcast(clusters[i].rank,0);
        multi_sim_comm.Bcast(clusters[i].time,0);
        multi_sim_comm.Bcast(clusters[i].population,0);
        multi_sim_comm.Bcast(clusters[i].cvs,0);
        multi_sim_comm.Bcast(clusters[i].sigma,0);
       }
     }
     multi_sim_comm.Barrier(); // Not sure this barrier is necessary, but it doesn't hurt to have it here (I think)
     
     // you'll rescale the distance to smoothly drag your system to at[i] (see function GenPot) over sithstepsup steps
     double resc=1.0;
     if (sithstepsup!=0) 
        resc=mod/sithstepsup;
     if (resc>1) resc=1;
     vector<vector<double> > potfor=GenPot(typot,clusters,resc,height,cv,delta0,limit_r); // ... and calculate the new bias from the present rank (done at every step)...
     //FIXME: Store information about the bias somewhere so you can read it at the beginning of the next iteration!
     potentials=potfor[0];
     forces=potfor[1];
  }
  
  
  for(unsigned i=0;i<getNumberOfArguments();++i) // It only iterates one time, but I still don't know how to 
  {
    ene += potentials[i];
    f=forces[i];
    totf2 += f*f;
    setOutputForce(i,f);
  }
  getPntrToComponent("bias")->set(ene); 
  getPntrToComponent("force2")->set(totf2);  
}
}
}
