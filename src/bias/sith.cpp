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
  for (unsigned i=0; i<getNumberOfArguments();i++)
  {
   wfile << "CV_" << i << " "; //How do you get the CV labels?
  }
  wfile << endl;
  wfile.close();
  //exit(0);
}

  struct values {
  double time;
  vector<double> cvs;
  vector<double> center;
  vector<double> sigma;
  values(double time,vector<double> cvs, vector<double> center, vector<double> sigma) 
    {
      this -> time = time;
      this -> cvs = cvs;
      this -> center = center;
      this -> sigma = sigma;
    }
  };
  
vector<values> values_raw;
vector<values> clusters;

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
  for (unsigned i=1;i<tokens.size();i++)
  {
      cvs.push_back(atof(tokens[i].c_str()));
  }
  vector<double> center;
  vector<double> sigma;
  values_raw.push_back(values(time, cvs, center, sigma));
 }
 //exit(0);
 return values_raw;
}


vector<values> cluster_snapshots(vector<values> & values_raw )
{
    
}
 
/*
 double GenGaussian(vector<Gaussian>)
 */

// Open cvfile for writing



void SITH::calculate(){
  
  const double pi=3.1415926535897;
  double ene = 0.0;
  double totf2 = 0.0;
  double f=0;
  
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
   for (unsigned i=0; i<getNumberOfArguments();i++)
   {
       double cv=getArgument(i);
       wfile << cv << " ";
   }
   wfile << endl;
   wfile.close();
  }
  
  
  // Check if the clustering has to be done
  if (step>=sithstride)
  {   
     int mod = step % sithstride;
     if (mod==0)
     {
         cout << "Reading the values in CV file: "<<endl;
         values_raw=getCVs(cvfile);
         //for (unsigned i=0; i<values_raw.size();i++) cout << values_raw[i].time << endl;
         //exit(0);
         cout << "and perform the clustering to generate new biases" << endl;
         clusters=cluster_snapshots(values_raw);
     }
      
  }
  //exit(0);
  
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
