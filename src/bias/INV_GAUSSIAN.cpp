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


using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS INV_GAUSSIAN
/*
 Applies a potential with the shape of an inverted gaussian:
 
 V=KAPPA*(1-exp(-((CV-AT)**2)/(2*(SIGMA**2))))
 
 This kind of potential is useful when the CV fluctuates a lot and harmonic
 potentials give forces too large for the integrator to handle
  
 USAGE:
 restr: INV_GAUSSIAN ARG=je AT=10.0 SIGMA=2.0 
 
*/
//+ENDPLUMEDOC

class INV_GAUSSIAN : public Bias{
  std::vector<double> at;
  std::vector<double> sigma;
  std::vector<double> kappa;
public:
  explicit INV_GAUSSIAN(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(INV_GAUSSIAN,"INV_GAUSSIAN")

void INV_GAUSSIAN::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","AT","The average of the distribution");
  keys.add("compulsory","SIGMA","The standard deviation of the distribution");
  keys.add("compulsory","KAPPA","Value at which the potential becomes flat");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}

INV_GAUSSIAN::INV_GAUSSIAN(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(getNumberOfArguments(),0),
sigma(getNumberOfArguments(),0.0)
{
  // Note sizes of these vectors are automatically checked by parseVector :-)
  parseVector("SIGMA",sigma);
  parseVector("AT",at);
  parseVector("KAPPA",kappa);
  checkRead();

  log.printf("  centered at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  
  log.printf("  with standard deviation");
  for(unsigned i=0;i<sigma.size();i++) log.printf(" %f",sigma[i]);
  log.printf("\n");
  
  log.printf("  becoming flat at");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  addComponent("bias"); componentIsNotPeriodic("bias");
  addComponent("force2"); componentIsNotPeriodic("force2");
}

void INV_GAUSSIAN::calculate(){
  const double pi=3.1415926535897;
  double ene = 0.0;
  double totf2 = 0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
      
    
    const double cv=difference(i,at[i],getArgument(i));
    //const double s=sigma[i];
    //const double k=kappa[i];
    
    ene=kappa[i]*(1-exp(-pow((cv-at[i]),2)/(2*pow(sigma[i],2))));
            
    double f = (1/(pow(sigma[i],3)*sqrt(pi)))*exp(-pow((cv-at[i]),2)/(2*pow(sigma[i],2)))*(cv-at[i]);
    totf2 += f * f;
    
    setOutputForce(i,f);
  }
  getPntrToComponent("bias")->set(ene); 
  getPntrToComponent("force2")->set(totf2);  
}

}
}
