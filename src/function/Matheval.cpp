/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "ActionRegister.h"
#include "Function.h"

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif

using namespace std;

namespace PLMD{
namespace function{


//+PLUMEDOC FUNCTION MATHEVAL
/*
Calculate a combination of variables using a matheval expression.

This action computes an  arbitrary function of one or more precomputed
collective variables. Arguments are chosen with the ARG keyword,
and the function is provided with the FUNC string. Notice that this
string should contain no space. Within FUNC, one can refer to the
arguments as x,y,z, and t (up to four variables provided as ARG).
This names can be customized using the VAR keyword (see examples below).

If you want a function that depends not only on collective variables
but also on time you can use the \subpage TIME action.

\attention
The MATHEVAL object only works if libmatheval is installed on the system and
PLUMED has been linked to it

\par Examples

The following input tells plumed to perform a metadynamics
using as a CV the difference between two distances.
\verbatim
dAB: DISTANCE ARG=10,12
dAC: DISTANCE ARG=10,15
diff: MATHEVAL ARG=dAB,dAC FUNC=y-x PERIODIC=NO
# notice: the previous line could be replaced with the following
# diff: COMBINE ARG=dAB,dAC COEFFICIENTS=-1,1
METAD ARG=diff WIDTH=0.1 HEIGHT=0.5 BIASFACTOR=10 PACE=100
\endverbatim
(see also \ref DISTANCE, \ref COMBINE, and \ref METAD).
Notice that forces applied to diff will be correctly propagated
to atoms 10, 12, and 15.
Also notice that since MATHEVAL is used without the VAR option
the two arguments should be referred to as x and y in the expression FUNC.
For simple functions
such as this one it is possible to use \ref COMBINE, which does
not require libmatheval to be installed on your system.

The following input tells plumed to print the angle between vectors
identified by atoms 1,2 and atoms 2,3
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\verbatim
DISTANCE LABEL=d1 ATOMS=1,2 COMPONENTS
DISTANCE LABEL=d2 ATOMS=2,3 COMPONENTS
MATHEVAL ...
  LABEL=theta
  ARG=d1.x,d1.y,d1.z,d2.x,d2.y,d2.z
  VAR=ax,ay,az,bx,by,bz
  FUNC=acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz))
  PERIODIC=NO
... MATHEVAL
PRINT ARG=theta
\endverbatim
(See also \ref PRINT and \ref DISTANCE).

*/
//+ENDPLUMEDOC


class Matheval :
  public Function
{
  void* evaluator;
  vector<void*> evaluator_deriv;
  vector<string> var;
  string func;
  vector<double> values;
  vector<char*> names;
public:
  explicit Matheval(const ActionOptions&);
  ~Matheval();
  void calculate();
  static void registerKeywords(Keywords& keys);
};

#ifdef __PLUMED_HAS_MATHEVAL
PLUMED_REGISTER_ACTION(Matheval,"MATHEVAL")

void Matheval::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","FUNC","the function you wish to evaluate");
  keys.add("optional","VAR","the names to give each of the arguments in the function.  If you have up to three arguments in your function you can use x, y and z to refer to them.  Otherwise you must use this flag to give your variables names.");
}

Matheval::Matheval(const ActionOptions&ao):
Action(ao),
Function(ao),
evaluator_deriv(getNumberOfArguments()),
values(getNumberOfArguments()),
names(getNumberOfArguments())
{
  parseVector("VAR",var);
  if(var.size()==0){
    var.resize(getNumberOfArguments());
    if(getNumberOfArguments()>3)
      error("Using more than 3 arguments you should explicitly write their names with VAR");
    if(var.size()>0) var[0]="x";
    if(var.size()>1) var[1]="y";
    if(var.size()>2) var[2]="z";
  }
  if(var.size()!=getNumberOfArguments())
    error("Size of VAR array should be the same as number of arguments");
  parse("FUNC",func);
  addValueWithDerivatives(); 
  checkRead();

  evaluator=evaluator_create(const_cast<char*>(func.c_str()));

  if(!evaluator) error("There was some problem in parsing matheval formula "+func);

  char **check_names;
  int    check_count;
  evaluator_get_variables(evaluator,&check_names,&check_count);
  if(check_count!=int(getNumberOfArguments())){
    string sc;
    Tools::convert(check_count,sc);
    error("Your function string contains "+sc+" arguments. This should be equal to the number of ARGs");
  }
  for(unsigned i=0;i<getNumberOfArguments();i++){
    bool found=false;
    for(unsigned j=0;j<getNumberOfArguments();j++){
      if(var[i]==check_names[j])found=true;
    }
    if(!found)
      error("Variable "+var[i]+" cannot be found in your function string");
  }

  for(unsigned i=0;i<getNumberOfArguments();i++)
    evaluator_deriv[i]=evaluator_derivative(evaluator,const_cast<char*>(var[i].c_str()));


  log.printf("  with function : %s\n",func.c_str());
  log.printf("  with variables :");
  for(unsigned i=0;i<var.size();i++) log.printf(" %s",var[i].c_str());
  log.printf("\n");
}

void Matheval::calculate(){
  for(unsigned i=0;i<getNumberOfArguments();i++) values[i]=getArgument(i);
  for(unsigned i=0;i<getNumberOfArguments();i++) names[i]=const_cast<char*>(var[i].c_str());
  setValue(evaluator_evaluate(evaluator,names.size(),&names[0],&values[0]));

  for(unsigned i=0;i<getNumberOfArguments();i++){
    setDerivative(i,evaluator_evaluate(evaluator_deriv[i],names.size(),&names[0],&values[0]));
  }
}

Matheval::~Matheval(){
  evaluator_destroy(evaluator);
  for(unsigned i=0;i<evaluator_deriv.size();i++)evaluator_destroy(evaluator_deriv[i]);
}

#endif

}
}


