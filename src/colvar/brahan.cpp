// This document is formatted for Doxygen
/**
\page AddingAColvar How to add a new collective variable
To implement a CV one you need to create a single cpp file called <i>ColvarName</i>.cpp in the directory src/colvar. 
If you use the following template for this file then the manual and the calls to the CV will be looked after automatically.
*/
//\verbatim
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/PDB.h"
#include "kabsch.h" // Kabsch algorithm implementation


#include <string>
#include <cmath>
#include <cassert>
#include <vector>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR NAME

/* \endverbatim
At this point you provide the description of your CV that will appear in the manual along with an description of the input file syntax and an example.  Merging new features of the code into the plumed main branch without proper documentation is punishable by death!  Some instructions as to how to format this information is provided here: \ref usingDoxygen
\verbatim
*/
//+ENDPLUMEDOC

/* We begin by declaring a class for your colvar.  This class inherits everything from the Colvar class.
      This ensures it has a label, a place to store its value, places to the store the values of the derivatives
      and that it can access the various atoms it will employ.
*/
class Brahan : public Colvar {
	double param;
private:
 vector<AtomNumber> reactant_atoms;//list of reactant atoms used for CV
 vector<AtomNumber> product_atoms;//list of product atoms used for CV
 vector<AtomNumber> ts_atoms;//list of ts atoms used for CV
 vector<AtomNumber> binding_atoms;//list of ts atoms used for CV
 vector<Vector> refreactant_pos;// reference coordinates Reactant for alignment
 double refreactant_com[3]; // reference coordinates of the center of mass of the Reactant

public:
 //---- This routine is used to create the descriptions of all the keywords used by your CV
 static void registerKeywords( Keywords& keys ); 
 //---- This is the constructor for your colvar.  It is this routine that will do all the reading.
     //  Hence it takes as input a line from the input file.
 explicit Brahan(const ActionOptions&);
 //---- This is the routine that will be used to calculate the value of the colvar, whenever its calculation is required.
   //    This routine and the constructor above must be present - if either of them are not the code will not compile.
  virtual void calculate();
};

 //------ The following command inserts your new colvar into plumed by inserting calls to your new
      //  routines into the parts of plumed where they are required.  This macro takes two arguments:
      //  The first is the name of your ColvarClass and the second is the keyword for your CV
        //(the first word in the input line for your CV).
PLUMED_REGISTER_ACTION(Brahan,"BRAHAN")

 //----- The following routine creates the documentation for the keyowrds used by your CV
void Brahan::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
 // ActionWithValue::registerKeywords(keys);
 // keys.remove("NUMERICAL_DERIVATIVES");
//\endverbatim
  keys.add("compulsory", "PARAM", "the value of input");
  keys.add("compulsory", "REACTANT", "the pdb of reactant");
  keys.add("compulsory", "PRODUCT", "the pdb of product");
  keys.add("compulsory", "TS", "the pdb of transition state");
  keys.add("compulsory", "BINDING", "the pdb of binding site");
//In here you should add all your descriptions of the keywords used by your colvar as well as descriptions of any components
//that you can use this colvar to calculate. Descriptions as to how to do this can be found here: \ref usingDoxygen

}

 //---- We now write the actual readin (constructor) and calculations routines for the colvar

Brahan::Brahan(const ActionOptions&ao):
 //------ This line sets up various things in the plumed core which colvars rely on.
PLUMED_COLVAR_INIT(ao)
{
  //vector<int> atoms;  /----- You almost always have atoms -----/
  parse("PARAM",param);
  string reactant_file;
  parse("REACTANT",reactant_file);
  string product_file;
  parse("PRODUCT",product_file);
  string ts_file;
  parse("TS",ts_file);
  string binding_file;
  parse("BINDING",binding_file);
  
  PDB reactant_pdb;
  if( !reactant_pdb.read(reactant_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) // using getAtoms from Plumedmain class
    error("missing input file " + reactant_file );

  reactant_atoms = reactant_pdb.getAtomNumbers();
  const std::vector<Vector> reactant_positions = reactant_pdb.getPositions();
  const std::vector<double> reactant_masses = reactant_pdb.getOccupancy();
  cout << " Reactant has ? elements " << reactant_atoms.size() << endl; 
 // cout << " reactant_positions[0][0] " << reactant_positions[0][0] << endl; 
  //cout << " index [0] " << reactant_atoms[0].index() << endl;

  PDB product_pdb;
  if( !product_pdb.read(product_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) // using getAtoms from Plumedmain class
    error("missing input file " + product_file );

  product_atoms = product_pdb.getAtomNumbers();
  const std::vector<Vector> product_positions = product_pdb.getPositions();
  const std::vector<double> product_masses = product_pdb.getOccupancy();
  cout << " Product has ? elements " << product_atoms.size() << endl;


  PDB ts_pdb;
  if( !ts_pdb.read(ts_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) // using getAtoms from Plumedmain class
    error("missing input file " + ts_file );

  ts_atoms = ts_pdb.getAtomNumbers();
  const std::vector<Vector> ts_positions = ts_pdb.getPositions();
  const std::vector<double> ts_masses = ts_pdb.getOccupancy();
  cout << " Transition state has ? elements " << ts_atoms.size() << endl;


  PDB binding_pdb;
  if( !binding_pdb.read(binding_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) // using getAtoms from Plumedmain class
    error("missing input file " + binding_file );

  binding_atoms = binding_pdb.getAtomNumbers();
  const std::vector<Vector> binding_positions = binding_pdb.getPositions();
  const std::vector<double> binding_masses = binding_pdb.getOccupancy();
  cout << " Binding state has ? elements " << binding_atoms.size() << endl;
//Insert code here to read the arguments of the CV here using plumed's parsing functionality.  N.B.  The label is read in already elsewhere.

//\verbatim
  checkRead();     //--- This command checks that everything on the input line has been read properly
//--- The following two lines inform the plumed core that we require space to store the value
    // of the CV and that the CV will act on a particular list of atoms.
  addValueWithDerivatives();
  setNotPeriodic();
// calculate the center of mass of substrate reference file
  double sub_mass_tot = 0.0;
  refreactant_com[0] =0.0;
  refreactant_com[1] =0.0;
  refreactant_com[2] =0.0;
  for (unsigned i=0; i< reactant_atoms.size() ; ++i)
    {
      refreactant_com[0] += reactant_masses[i] * reactant_positions[i][0];
      refreactant_com[1] += reactant_masses[i] * reactant_positions[i][1];
      refreactant_com[2] += reactant_masses[i] * reactant_positions[i][2];
  refreactant_pos.push_back( Vector(reactant_positions[i][0], reactant_positions[i][1], reactant_positions[i][2]) );

      sub_mass_tot += reactant_masses[i];
    }

  refreactant_com[0] /= sub_mass_tot;
  refreactant_com[1] /= sub_mass_tot;
  refreactant_com[2] /= sub_mass_tot;

  cout << " refreactant_com " << refreactant_com[0] << " " << refreactant_com[1] << " " << refreactant_com[2] << endl;




// all_atoms vector = substrate + binding site
  vector<AtomNumber> all_atoms( reactant_atoms.size() + binding_atoms.size() );
   
  for ( unsigned i = 0; i < reactant_atoms.size(); ++i )
    all_atoms[i] = reactant_atoms[i];

  
  for ( unsigned i = 0; i < binding_atoms.size(); ++i )
    all_atoms[reactant_atoms.size()+i] = binding_atoms[i];




  requestAtoms(all_atoms);
   
/*--- For a number of the free energy methods in plumed it is necessary to calculate the
      distance between two points in CV space.  Obviously, for periodic CVs one must take
      periodicities into account when calculating distances and use the minimum image
      convention in distance calculations.  Hence, we set the periodicity of the cv using
      the following two lines.
*/  
 //getValue("")->setPeridodicity(true);  // Set this true if the CV is periodic otherwise set if false.
  // getValue("")->setDomain(min,max);     // This routine is only required if the function is periodic.  It sets the minimum and maximum values of the colvar.
  cout << "PARAM " << param << endl;
  log.printf(" PARAM %f \n" , param);
}

void Brahan::calculate(){
//--- These are the things you must calculate for any cv ---/
  double brahan_val=0.0;
  
  vector<double> R_x; // array with update of x coordinates of each reference substate Reactant according to translation/rotation
  R_x.reserve(reactant_atoms.size());

  vector<double> R_y; // array with update of y coordinates of each reference substate Reactant according to translation/rotation
  R_y.reserve(reactant_atoms.size());

  vector<double> R_z; // array with update of z coordinates of each reference substate Reactant according to translation/rotation
  R_z.reserve(reactant_atoms.size());





 // Tensor boxDerivatives;      /--- The derivative of the cv with respect to the box vectors ----/
 // vector<double> derivatives; /--- The derivative of the cv with respect to the atom positions ---/

//Insert the code to calculate your cv, its derivatives and its contribution to the virial here. Please use, where possible, the library of tools described in \ref TOOLBOX.

//\verbatim
///---- Having calculated the cv, its derivative and the contribution to the virial you now
  //    transfer this information to the plumed core using the following three commands. 
  //for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
  //setBoxDerivatives(boxDerivatives);
  //setValue(param);
  cout << "Beginning  Brahan calculation "<< endl;
  cout << ".............................. "<< endl;
  unsigned n_reactantatoms = reactant_atoms.size();
//  unsigned n_bindingatoms = binding_atoms.size();
  double sub_com[3] = {0.0, 0.0 ,0.0};
  double sub_mass=0.0;

// calculate new center of mass of substate
  for (unsigned j =0; j < n_reactantatoms ; ++j)
    {
      double j_mass = getMass(j);     
      Vector j_pos = getPosition(j);
      sub_com[0] += j_mass * j_pos[0];
      sub_com[1] += j_mass * j_pos[1];
      sub_com[2] += j_mass * j_pos[2];
      sub_mass += j_mass;
    }
  sub_com[0] /= sub_mass;
  sub_com[1] /= sub_mass;
  sub_com[2] /= sub_mass;
  cout << " sub_com " << sub_com[0] << " " << sub_com[1] << " " << sub_com[2] << endl;

  double ref_xlist[n_reactantatoms][3]; // coordinate of referece substrate reactant
  double mov_xlist[n_reactantatoms][3]; // new coordinate of substate
  double mov_com[3] = {0.0,0.0,0.0}; // new center of mass of substrate
  double mov_to_ref[3]; // vector between the com of move and ref substate
  double rotmat[3][3]; //the rotation matrix for least-squares fit
  double rmsd = 0.0;

  double mov_mass_tot=0.0;

  for (unsigned i=0; i < n_reactantatoms ; ++i)
    {
      ref_xlist[i][0] = refreactant_pos[i][0];
      ref_xlist[i][1] = refreactant_pos[i][1];
      ref_xlist[i][2] = refreactant_pos[i][2];
      Vector i_pos = getPosition( i );
      mov_xlist[i][0] = i_pos[0];
      mov_xlist[i][1] = i_pos[1];
      mov_xlist[i][2] = i_pos[2];
      double i_mass = getMass( i );
      mov_mass_tot += i_mass;
      mov_com[0] += i_mass * i_pos[0];
      mov_com[1] += i_mass * i_pos[1];
      mov_com[2] += i_mass * i_pos[2];
    }
  // Set mov_com and mov_to_ref
  mov_com[0] /= mov_mass_tot;
  mov_com[1] /= mov_mass_tot;
  mov_com[2] /= mov_mass_tot;

  mov_to_ref[0] = refreactant_com[0] - mov_com[0];
  mov_to_ref[1] = refreactant_com[1] - mov_com[1];
  mov_to_ref[2] = refreactant_com[2] - mov_com[2];

  rotmat[0][0] = 1.0;
  rotmat[0][1] = 0.0;
  rotmat[0][2] = 0.0;
  rotmat[1][0] = 0.0;
  rotmat[1][1] = 1.0;
  rotmat[1][2] = 0.0;
  rotmat[2][0] = 0.0;
  rotmat[2][1] = 0.0;
  rotmat[2][2] = 1.0;

  calculate_rotation_rmsd( ref_xlist, mov_xlist, n_reactantatoms, mov_com, mov_to_ref, rotmat, &rmsd  );

  cout << "Just called calculated_rotation_rmsd" << endl;
  cout << "rotmat elements  of substate :" << endl;
  cout << rotmat[0][0] << " " << rotmat[0][1] << " " << rotmat[0][2] << endl;
  cout << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << endl;
  cout << rotmat[2][0] << " " << rotmat[2][1] << " " << rotmat[2][2] << endl; 
  cout << "rmsd " << rmsd << endl;

// translate Reactant_ref using new com of mov_substate
  double new_refR_com_x;
  double new_refR_com_y;
  double new_refR_com_z;
  new_refR_com_x = sub_com[0];
  new_refR_com_y = sub_com[1];
  new_refR_com_z = sub_com[2];

 // cout << " New reference R " << new_refR_com_x << " " << new_refR_com_y << " " << new_refR_com_z << endl;

  // Now rotate the reference Reactant position at origin and then translate to new com

  for(unsigned i=0;i< n_reactantatoms;i++)
    {
      R_x[i]  = ( rotmat[0][0]*(refreactant_pos[i][0]) + rotmat[1][0]*(refreactant_pos[i][1]) + rotmat[2][0]*(refreactant_pos[i][2]) ) + new_refR_com_x;
      R_y[i]  = ( rotmat[0][1]*(refreactant_pos[i][0]) + rotmat[1][1]*(refreactant_pos[i][1]) + rotmat[2][1]*(refreactant_pos[i][2]) ) + new_refR_com_y;
      R_z[i]  = ( rotmat[0][2]*(refreactant_pos[i][0]) + rotmat[1][2]*(refreactant_pos[i][1]) + rotmat[2][2]*(refreactant_pos[i][2]) ) + new_refR_com_z;
      cout << "new reference R coordinate of atom " << i << " "  << R_x[i] << " " << R_y[i] << " " << R_z[i] << endl;

     }
  
 



  setValue(brahan_val);
}
//\endverbatim

/*\section multicvs Mult-component CVs
To avoid code duplication, and in some cases computational expense, plumed has functionality so that a single line in input can calculate be used to calculate multiple components for a CV.  For example, PATH computes the distance along the path,\f$s\f$, and the distance from the path, \f$z\f$.  Alternatively, a distance can give one the \f$x\f$, \f$y\f$ and \f$z\f$ components of the vector connecting the two atoms.  You can make use of this functionality in your own CVs as follows:
*
- In the constructor we create an additional value for the CV by adding the call PLMD::addValueWithDerivative("new") as well as PLMD::addValueWithDerivatives(””).  In addtion set any periodicity for our component using getValue("new")->setPeridicity() and getValue("new")->setDomain(min,max).  If this CV is called plum in our input file we can now use both plum and plum.new in any of the functions/methods in plumed.
- Obviously in calculate we now must provide functionality to calculate the values, boxDerivatives and the atom derivatives for both our original plum and its component plum.new. Furthermore, all of this data must be transferred to the plumed core.  This is done by the following code:
Here we transfer the value, box derivatives and atomic derivatives for plum.
\verbatim
for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
setBoxDerivatives(boxDerivatives);
setValue(cv_val);
\endverbatim
Here we transfer the value, box derivatives and atomic derivatives for plum.new.
\verbatim
Value* nvalue=getValue("new");
for(int i=0;i<nderivatives.size();i++){ setAtomsDerivatives(nvalue i,nderivatives[i]); }
setBoxDerivatives(nvalue,nboxDerivatives);
setValue(nvalue,ncv_val);
\endverbatim
Please only use this functionality for CVs that are VERY similar.
*/

}
}

