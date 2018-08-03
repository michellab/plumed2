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
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iterator>

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
      
  struct Gaussian {
   vector<double> center;
   vector<double> sigma;
   double height;
   bool   multivariate; // this is required to discriminate the one dimensional case 
   vector<double> invsigma;
   Gaussian(const vector<double> & center,const vector<double> & sigma,double height, bool multivariate ):
   center(center),sigma(sigma),height(height),multivariate(multivariate),invsigma(sigma)
      {
       // to avoid troubles from zero element in flexible hills
       for(unsigned i=0;i<invsigma.size();++i)abs(invsigma[i])>1.e-20?invsigma[i]=1.0/invsigma[i]:0.; 
      }
  };
  
  double height; // Factor that will multiply the height of the gaussians
  int sithstride; // Stride to perform the cv clustering and generate a new biasing potential
  double dc; // dc value for clustering (See Laio2014))
  double delta0; // delta0 value for clustering (see laio2014)
  int cvstride; // stride to print the value of the cvs in the file that will be read for clustering
  string sithfile; // The name of the file that will contain the clusters found at each sithstride (only write)
  string cvfile; // The name of the file that will contain the cvs involved in the taboo search (read+write)
  
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
  keys.add("compulsory","CVSTRIDE","Frequency with which the snapshots are clustered and the bias is updated");
  keys.add("compulsory","HEIGHT","Factor that scales the height of the gaussians");
  keys.add("compulsory","DC","dc value for the clustering");
  keys.add("compulsory","DELTA0","delta0 value for the clustering");
  keys.add("compulsory","SITHFILE","The name of the file that will contain the clusters found at each sithstride (only write)");
  keys.add("compulsory","CVFILE","The name of the file that will contain the cvs involved in the taboo search (read+write)");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}


SITH::SITH(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
//height(getNumberOfArguments(),0),
//sithstride(getNumberOfArguments(),0)
{
  // Note sizes of these vectors are automatically checked by parseVector :-)
  parse("HEIGHT",height);
  parse("SITHSTRIDE",sithstride);
  parse("CVSTRIDE",cvstride);
  parse("DC",dc);
  parse("DELTA0",delta0);
  parse("CVFILE",cvfile);
  parse("SITHFILE",sithfile);
  checkRead();
   
  printf("The CVs are going to be printed every %i steps.\n", cvstride);
  printf("Clustering is going to be performed every %i steps.\n", sithstride);
  printf("The height of the gaussians will be rescaled by a factor of %f.\n", height);
  printf("Please read and cite: Rodriguez, A.; Laio, A.; Science (2014) 344(6191) p.1496\n");
  
  //exit(0);
  addComponent("bias"); componentIsNotPeriodic("bias");
  addComponent("force2"); componentIsNotPeriodic("force2");
  
  ofstream wfile;
  wfile.open(cvfile.c_str());
  wfile << "Time ";
  for (int i=0; i<getNumberOfArguments();i++)
  {
   wfile << "CV_" << i << " "; //How do you get the CV labels?
  }
  wfile << endl;
  wfile.close();
  //exit(0);
}


//Data struct we will use to perform the clustering
 struct Laio
  {
    int snapshot;
    double rho;
    int cluster;
    unsigned nnhd;
    double delta;
  };
  
    bool mycmp_rho(const Laio& a, const Laio& b)
  {
      return a.rho > b.rho;
  };

  
  //Data struct we will use to store the time, cv values, population and width of each cluster
  struct values {
  double time;
  vector<double> cvs;
  int  population;
  vector<double> sigma;
  values(double time,vector<double> cvs, int population, vector<double> sigma) 
    {
      this -> time = time;
      this -> cvs = cvs;
      this -> population = population;
      this -> sigma = sigma;
    }
  };
  
vector<values> values_raw;
vector<values> clusters;

/////////////////////////////////////////////////////////
// this reads the cvs from CVFILE                      //
/////////////////////////////////////////////////////////
vector<values> getCVs(string CV_file)
{
 vector<values> values_raw;
 FILE* fp=fopen(CV_file.c_str(),"r");
 string line;
 while (Tools::getline(fp, line))
 {
  if (line[0] == 'T') continue;
  istringstream iss(line);
  istream_iterator<string> beg(iss), end;
  vector<string> tokens(beg, end);
  double time=atof(tokens[0].c_str()); 
  vector<double> cvs;
  for (int i=1;i<tokens.size();i++)
  {
      cvs.push_back(atof(tokens[i].c_str()));
  }
  int population=0;
  vector<double> sigma;
  values_raw.push_back(values(time, cvs, population, sigma));
 }
 //exit(0);
 return values_raw;
}

/////////////////////////////////////////////////////////
// BEGINNING OF THE LAIO FUNCTION                      //
/////////////////////////////////////////////////////////
vector<values> cluster_snapshots(vector<values> & values_raw, double dc, double delta0)
{
  vector<Laio> vec(values_raw.size());
  vector<values> clusters_raw;
  vector<values> clusters;
  double dc2=dc*dc;
  #pragma omp parallel shared(vec,clusters_raw,clusters,dc2)
  {
    //cout << "calculating density" << endl;
      #pragma omp for
      for (int i=0; i<values_raw.size();i++)
      {
          vec[i].snapshot=i;
          for (int j=0;j<values_raw.size();j++)
          {
              if (j==i) continue;
              double r2=0;
              for (int k=0; k<values_raw[i].cvs.size();k++) 
                  r2=+pow((values_raw[j].cvs[k]-values_raw[i].cvs[k]),2);
              if (r2<dc2) vec[i].rho+=1;
          }
      }
  
  
  // Sorting vector by decreasing density 
  #pragma omp single
   {
    sort(vec.begin(),vec.end(),mycmp_rho);
   }
  
  //cout << "Calculating delta" << endl;
  #pragma omp for
  for (int i=0; i<values_raw.size(); i++)
  {
      vec[i].cluster=-1; // This is initialised at 0 but don't want random lonely points in cluster 0
      vec[i].nnhd=-1;
      int snap_i=vec[i].snapshot;
      if (i==0)
      {
          double maxdist2=0.;
          for (int j=i+1; j<values_raw.size();j++)
          {
           double r2=0;
           for (int k=0; k<values_raw[i].cvs.size();k++) 
              r2=+pow((values_raw[j].cvs[k]-values_raw[i].cvs[k]),2);
           if (r2>maxdist2)
             {
              vec[i].nnhd=snap_i;
              vec[i].delta=sqrt(r2);
              maxdist2=r2;
             }
          }
           continue;
      }
      else
      {
          double mindist2=99999999999999.;
          double nnhd_rho=0.;
          for (int j=0; j<i; j++)
          {
            int snap_j=vec[j].snapshot;
            double r2=0;
            for (int k=0; k<values_raw[i].cvs.size();k++) 
                r2=+pow((values_raw[snap_j].cvs[k]-values_raw[snap_i].cvs[k]),2);
            if (r2<mindist2)
            {
                vec[i].nnhd=snap_j;
                nnhd_rho=vec[j].rho;
                vec[i].delta=sqrt(r2);
                mindist2=r2;
            }
          }
      }
  }
  
  #pragma omp single
  {
 // cout << "Generating an empty values struct" << endl;
  double t_center=0.;
  vector<double> cvs_center(values_raw[0].cvs.size(),0);
  int pop=0.;
  vector<double> sig_center(values_raw[0].cvs.size(),0);
  
  
  //cout << "Choosing cluster centres" << endl;
  int nclust=0;
  for (int i=0; i<values_raw.size();i++)
   {
      if (vec[i].delta>delta0)
       {
          vec[i].cluster=nclust;
          clusters_raw.push_back(values(t_center,cvs_center,pop,sig_center));
          nclust++;
       }
   }
      
  if (clusters_raw.size()==0) // in some cases delta0 can be so high that no clusters are found
   {
     vec[0].cluster=nclust;
     clusters_raw.push_back(values(t_center,cvs_center,pop,sig_center));
   }
  
  
  //cout << "found " << nclust << " cluster centers" << endl;
  //cout << "clusters_raw has " << clusters_raw.size() << " elements" << endl;
   
  //cout << "assigning points to clusters" << endl;
  

    vector<double> maxdist(clusters_raw.size(),0);
    for (int i=0; i<values_raw.size();i++)
    {
        if (vec[i].cluster != -1)
        {
         int snap_center = vec[i].snapshot;
         //cout << "Assigning point " << snap_center << " as center of cluster " << vec[i].cluster << endl;
         //find the snapshot with time vec[i].time
         vector<double> cvs_center = values_raw[snap_center].cvs;
         //cout << "initialising population at 1" << endl;
         int pop = 1;
         //cout << "assigning time" << endl;
         clusters_raw[vec[i].cluster].time = values_raw[snap_center].time;
         //cout << "assigning cvs" << endl;
         clusters_raw[vec[i].cluster].cvs = cvs_center;
         //cout << "assigning pop" << endl;
         clusters_raw[vec[i].cluster].population = pop;
         //cout << "cluster center initialised" << endl;
        }
        else
        {
         int snap_i = vec[i].snapshot;
         //cout << "looking at snapshot " << snap_i << " which corresponds to time: " << values_raw[snap_i].time << endl;
         for (int j=0; j<i;j++) // points j have a higher density than points i
           {
             int snap_j = vec[j].snapshot;
             //cout << "snap_j: " << vec[j].snapshot << " snap_i.nnhd: " << vec[i].nnhd << endl;
             if (vec[i].nnhd != vec[j].snapshot) 
                 continue;
             //cout << "comparing with snapshot " << snap_j << " which corresponds to time: " << values_raw[snap_j].time << endl;
             vec[i].cluster=vec[j].cluster;
             clusters_raw[vec[j].cluster].population += 1;
             double r2=0; // Finding the furthest point in the cluster;
             for (int k=0; k<values_raw[snap_i].cvs.size();k++)
             {   
                 r2 += pow((values_raw[snap_i].cvs[k]-values_raw[snap_j].cvs[k]),2);
             }
             //cout << "r2 for cluster " << vec[j].cluster << "is " << r2 << endl;
             if (r2>maxdist[vec[j].cluster])
             {
              clusters_raw[vec[j].cluster].sigma=values_raw[snap_i].cvs;
             }
           }
         
        }
    }
  }
  //cout << Correcting the sigma of the gaussians
  #pragma omp for
  for (unsigned i=0; i<clusters_raw.size();i++)
   {
    for (int k=0; k<clusters_raw[i].cvs.size();k++)
    {
      if (clusters_raw[i].sigma[k]>0) clusters_raw[i].sigma[k]= abs(clusters_raw[i].sigma[k] - clusters_raw[i].cvs[k]);
    }
   }
  /*
    for (int i=0; i<values_raw.size();i++)
      {
          cout << vec[i].snapshot << " " << vec[i].rho << " " << vec[i].delta << " " << vec[i].nnhd << " " << vec[i].cluster << endl;
      }
  */
  }
  return clusters_raw;
}

void SITH::calculate(){
  
  int step=getStep();
  double time=getTime();

  //Check if values have to be printed in cvfile
  int mod = step % cvstride;
  if (mod==0)
  {
   double time=getTime();
   ofstream wfile;
   wfile.open(cvfile.c_str(),std::ios_base::app);
   wfile << time << " ";
   for (int i=0; i<getNumberOfArguments();i++)
   {
       double cv=getArgument(i);
       wfile << cv << " ";
   }
   wfile << endl;
   wfile.close();
  }
  
  
  // Clustering snapshots and generating new potentials if necessary
  if (step>=sithstride)
  {   
     int mod = step % sithstride;
     if (mod==0)
     {
         //cout << "Reading the values in CV file: ";
         values_raw=getCVs(cvfile);
         //for (int i=0; i<values_raw.size();i++) cout << values_raw[i].time << endl;
         //exit(0);
         //cout << "and performing the clustering to generate new biases" << endl;
         clusters=cluster_snapshots(values_raw,dc,delta0);
         //cout << "------------------------------------" << endl;
         //cout << "Clustering at time = " << time<<endl;
         /*
         for (int i=0; i<clusters.size();i++) 
         {
             
             cout << "Cluster " << i <<":" << endl;
             cout << "Time: " << clusters[i].time << endl;
             cout << "CVs: ";
             for (int j=0;j<clusters[i].cvs.size();j++) cout << clusters[i].cvs[j] << " ";
             cout << endl;
             cout << "SIGMAs: ";
             for (int j=0;j<clusters[i].sigma.size();j++) cout << clusters[i].sigma[j] << " ";
             cout << endl;
             cout << "Population: " << clusters[i].population << endl;
         }
         cout << "------------------------------------" << endl;
          */
     }
   //exit(0);   
  }
  
  const double pi=3.1415926535897;
  double ene = 0.0;
  double totf2 = 0.0;
  double f=0.;
  
  for(unsigned i=0;i<getNumberOfArguments();++i) // It only iterates one time, but I still don't know how to 
  {
    setOutputForce(i,f);
  }
  getPntrToComponent("bias")->set(ene); 
  getPntrToComponent("force2")->set(totf2);  
}
//Close CVfile for writing
}
}
