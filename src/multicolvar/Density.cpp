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
#include "MultiColvar.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR DENSITY
/*
Calculate functions of the density of atoms as a function of the box.  This allows one to calculate
the number of atoms in half the box.

\par Examples 

The following example calculates the number of atoms in one half of the simulation box. 

\verbatim
DENSITY SPECIES=1-100 LABEL=d
SUBCELL ARG=d XLOWER=0.0 XUPPER=0.5 LABEL=d1
PRINT ARG=d1.* FILE=colvar1 FMT=%8.4f
\endverbatim

*/
//+ENDPLUMEDOC


class Density : public MultiColvar {
public:
  static void registerKeywords( Keywords& keys );
  explicit Density(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
  /// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
  bool isDensity() const { return true; }
  bool hasDifferentiableOrientation() const { return true; }
//  void addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
//                                        const std::vector<double>& weight, MultiColvarFunction* func ){}
  void getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices );
//  unsigned getNumberOfQuantities();
  void getValueForTask( const unsigned& iatom, std::vector<double>& vals );
};

PLUMED_REGISTER_ACTION(Density,"DENSITY")

void Density::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("SPECIES"); 
}

Density::Density(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  int nat=1; readAtoms( nat ); 
  // And check everything has been read in correctly
  checkRead(); 
}

double Density::compute( const unsigned& tindex, AtomValuePack& myvals ) const {
  return 1.0;
}

void Density::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
   indices[jstore]=0; 
}

// unsigned Density::getNumberOfQuantities(){
//    return 2;
// }

void Density::getValueForTask( const unsigned& iatom, std::vector<double>& vals ){
   plumed_dbg_assert( vals.size()==2 ); vals[0]=vals[1]=1.0;
}

}
}

