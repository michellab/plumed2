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
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iterator> // std::istream_iterator

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
class parameter
{
public: 
  void readParamFile(string &param_file, string &RParam_file, string &PParam_file, string &TSParam);
  
  vector<double> R_sigma;// sigma LJ parameter for Reactant
  vector<double> R_epsilon;// epsilon LJ parameter for Reactant
  vector<double> P_sigma;// sigma LJ parameter for Product
  vector<double> P_epsilon;// epsilon LJ parameter for Product
  vector<double> TS_sigma;// sigma LJ parameter for Product
  vector<double> TS_epsilon;// sigma LJ parameter for Product
  vector<double> bindSite_sigma;// sigma LJ parameter for binding site
  vector<double> bindSite_epsilon;// epsilon LJ parameter for binding site
  vector<double> R_charge; //charge of reactant
  vector<double> P_charge;
  vector<double> TS_charge;
  vector<double> bindSite_charge;
  vector<string> R_atomName;// atom name of reactant
  vector<string> P_atomName;
  vector<string> TS_atomName;
  vector<string> bindSite_atomName;

};


void parameter::readParamFile(string &bindingParam_file, string &RParam_file, string &PParam_file, string &TSParam_file)
{

  string files[4] = {bindingParam_file, RParam_file, PParam_file, TSParam_file} ;
  for ( int i = 0; i < 4 ; i++)
  {
  //   cout << "file "<<  files[i] << endl;
  
   
     FILE * pFile;
     pFile = fopen(files[i].c_str(), "r" );
     if (!pFile)
       cout << "Cannot find binding parameter file" << endl;
     else
       cout << "reading parameter file " << files[i] << endl;
       string line;
       while(Tools::getline(pFile, line))
       {
         if (line[0] == '#')
           continue;
           //cout << line[0] << endl;
         istringstream iss(line);
         istream_iterator<string> beg(iss), end;
         vector<string> tokens(beg, end);
         // tokens[0] = index, tokens[1]= name, tokens[2]=type, tokens[3]=sigma, tokens[4]=epsilon, tokens[5]=charge
         //cout <<"index " << tokens[0]<< " sigma " << tokens[3] << " epsilon "<< tokens[4] <<endl;
         if (i==0)
         { 
           bindSite_sigma.push_back( atof(tokens[3].c_str()) );
           bindSite_epsilon.push_back( atof(tokens[4].c_str()) );
           bindSite_charge.push_back( atof(tokens[5].c_str()) );
           bindSite_atomName.push_back( tokens[1] );
           //cout << tokens[0] << endl;
         }
         else if (i==1)
         { 
           R_sigma.push_back( atof(tokens[3].c_str()) );
           R_epsilon.push_back( atof(tokens[4].c_str()) );
           R_charge.push_back( atof(tokens[5].c_str()) );
           R_atomName.push_back( tokens[1] );
           //cout << tokens[0] << endal;
         }
         else if (i==2)
         { 
           P_sigma.push_back( atof(tokens[3].c_str()) );
           P_epsilon.push_back( atof(tokens[4].c_str()) );
           P_charge.push_back( atof(tokens[5].c_str()) );
           P_atomName.push_back( tokens[1] );
           //cout << tokens[0] << endl;
         }
         else if (i==3)
         { 
           TS_sigma.push_back( atof(tokens[3].c_str()) );
           TS_epsilon.push_back( atof(tokens[4].c_str()) );
           TS_charge.push_back( atof(tokens[5].c_str()) );
           TS_atomName.push_back( tokens[1] );
           //cout << tokens[0] << endl;
         }
       }
  }
}

class Brahan : public Colvar {
	double param;
private:
 vector<AtomNumber> reactant_atoms;//list of reactant atoms used for CV
 vector<AtomNumber> product_atoms;//list of product atoms used for CV
 vector<AtomNumber> ts_atoms;//list of ts atoms used for CV
 vector<AtomNumber> binding_atoms;//list of ts atoms used for CV
 vector<Vector> refreactant_pos;// reference coordinates Reactant for alignment
 vector<Vector> refproduct_pos;// reference coordinates Reactant for alignment
 vector<Vector> refts_pos;// reference coordinates Reactant for alignment
 double refreactant_com[3]; // reference coordinates of the center of mass of the Reactant
 double refproduct_com[3]; // reference coordinates of the center of mass of the Reactant
 double refts_com[3]; // reference coordinates of the center of mass of the Reactant
 int LJ_rule; // the combining rule of LJ parameter

 parameter params;// parameter for energy calculation
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
  keys.add("compulsory", "BINDING_PARAM", "the parameter of binding site");
  keys.add("compulsory", "R_PARAM", "the parameter of reactant");
  keys.add("compulsory", "P_PARAM", "the parameter of product");
  keys.add("compulsory", "TS_PARAM", "the parameter of product");
  keys.add("compulsory", "REACTANT", "the pdb of reactant");
  keys.add("compulsory", "PRODUCT", "the pdb of product");
  keys.add("compulsory", "TS", "the pdb of transition state");
  keys.add("compulsory", "BINDING", "the pdb of binding site");
  keys.add("optional", "LJ_RULE", "the combinding rule of LJ parameter (1)defalt, (2) use Lorentz-Berthebt rule");
//In here you should add all your descriptions of the keywords used by your colvar as well as descriptions of any components
//that you can use this colvar to calculate. Descriptions as to how to do this can be found here: \ref usingDoxygen

}

 //---- We now write the actual readin (constructor) and calculations routines for the colvar







Brahan::Brahan(const ActionOptions&ao):
 //------ This line sets up various things in the plumed core which colvars rely on.
PLUMED_COLVAR_INIT(ao)
{
  //vector<int> atoms;  /----- You almost always have atoms -----/
  string bindingParam_file;
  parse("BINDING_PARAM",bindingParam_file);
  string reactant_file;
  parse("REACTANT",reactant_file);
  string product_file;
  parse("PRODUCT",product_file);
  string ts_file;
  parse("TS",ts_file);
  string binding_file;
  parse("BINDING",binding_file);
  string RParam_file;
  parse("R_PARAM",RParam_file);
  string PParam_file;
  parse("P_PARAM",PParam_file);
  string TSParam_file;
  parse("TS_PARAM",TSParam_file);
  string LJ_rule_string;
  parse("LJ_RULE",LJ_rule_string);
  if (LJ_rule_string.length() == 0)
    LJ_rule_string = "null";
  if (LJ_rule_string != "null")
    LJ_rule=atoi(LJ_rule_string.c_str());
  else
    LJ_rule= 0;


  PDB reactant_pdb;
  if( !reactant_pdb.read(reactant_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) // using getAtoms from Plumedmain class
    error("missing input file " + reactant_file );
  
  reactant_atoms = reactant_pdb.getAtomNumbers();
  const std::vector<Vector> reactant_positions = reactant_pdb.getPositions();
  const std::vector<double> reactant_masses = reactant_pdb.getOccupancy();
  cout << " Reactant has ? elements " << reactant_atoms.size() << endl; 
 // cout << " reactant_positions[0][0] " << reactant_positions[0][0] << endl; 
  //cout << " index [0] " << reactant_atoms[0].index() << endl;
 // for (unsigned i=0; i< reactant_atoms.size() ; ++i){
  
//  string Rname = reactant_pdb.getAtomName( reactant_atoms[i] ); 
//  cout << "reactant atom " << i <<  " name  " << Rname << endl;
 
 //  }
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
  double subp_mass_tot = 0.0;
  refproduct_com[0] =0.0;
  refproduct_com[1] =0.0;
  refproduct_com[2] =0.0;
  double subts_mass_tot = 0.0;
  refts_com[0] =0.0;
  refts_com[1] =0.0;
  refts_com[2] =0.0;

  for (unsigned i=0; i< reactant_atoms.size() ; ++i)
    {
      refreactant_com[0] += reactant_masses[i] * reactant_positions[i][0];
      refreactant_com[1] += reactant_masses[i] * reactant_positions[i][1];
      refreactant_com[2] += reactant_masses[i] * reactant_positions[i][2];
      refreactant_pos.push_back( Vector(reactant_positions[i][0], reactant_positions[i][1], reactant_positions[i][2]) );

      sub_mass_tot += reactant_masses[i];

      refproduct_com[0] += product_masses[i] * product_positions[i][0];
      refproduct_com[1] += product_masses[i] * product_positions[i][1];
      refproduct_com[2] += product_masses[i] * product_positions[i][2];
      refproduct_pos.push_back( Vector(product_positions[i][0], product_positions[i][1], product_positions[i][2]) );

      subp_mass_tot += product_masses[i];
     
      refts_com[0] += ts_masses[i] * ts_positions[i][0];
      refts_com[1] += ts_masses[i] * ts_positions[i][1];
      refts_com[2] += ts_masses[i] * ts_positions[i][2];
      refts_pos.push_back( Vector(ts_positions[i][0], ts_positions[i][1], ts_positions[i][2]) );

      subts_mass_tot += ts_masses[i];


    }

  refreactant_com[0] /= sub_mass_tot;
  refreactant_com[1] /= sub_mass_tot;
  refreactant_com[2] /= sub_mass_tot;
   

//  cout << " refreactant_com " << refreactant_com[0] << " " << refreactant_com[1] << " " << refreactant_com[2] << endl;

  refproduct_com[0] /= subp_mass_tot;
  refproduct_com[1] /= subp_mass_tot;
  refproduct_com[2] /= subp_mass_tot;

 // cout << " refproduct_com " << refproduct_com[0] << " " << refproduct_com[1] << " " << refproduct_com[2] << endl;


  refts_com[0] /= subts_mass_tot;
  refts_com[1] /= subts_mass_tot;
  refts_com[2] /= subts_mass_tot;

 // cout << " refts_com " << refts_com[0] << " " << refts_com[1] << " " << refts_com[2] << endl;




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
  // read paramter file
  params.readParamFile(bindingParam_file, RParam_file, PParam_file, TSParam_file);


}

void Brahan::calculate(){
//--- These are the things you must calculate for any cv ---/
  double brahan_val=0.0;
  





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
//  cout << " sub_com " << sub_com[0] << " " << sub_com[1] << " " << sub_com[2] << endl;

  double refR_xlist[n_reactantatoms][3]; // coordinate of referece substrate reactant
  double refP_xlist[n_reactantatoms][3]; // coordinate of referece substrate product
  double refTS_xlist[n_reactantatoms][3]; // coordinate of referece substrate ts
  double mov_xlist[n_reactantatoms][3]; // new coordinate of substate
  double mov_com[3] = {0.0,0.0,0.0}; // new center of mass of substrate
  double movR_to_ref[3]; // vector between the com of move and ref reactant substate
  double movP_to_ref[3]; // vector between the com of move and ref product substate
  double movTS_to_ref[3]; // vector between the com of move and ref ts substate
  double rotmatR[3][3]; //the rotation matrix for least-squares fit
  double rotmatP[3][3]; //the rotation matrix for least-squares fit
  double rotmatTS[3][3]; //the rotation matrix for least-squares fit
  double rmsdR = 0.0;
  double rmsdP = 0.0;
  double rmsdTS = 0.0;

  double mov_mass_tot=0.0;

  for (unsigned i=0; i < n_reactantatoms ; ++i)
    {
      refR_xlist[i][0] = refreactant_pos[i][0];
      refR_xlist[i][1] = refreactant_pos[i][1];
      refR_xlist[i][2] = refreactant_pos[i][2];

      refP_xlist[i][0] = refproduct_pos[i][0];
      refP_xlist[i][1] = refproduct_pos[i][1];
      refP_xlist[i][2] = refproduct_pos[i][2];

      refTS_xlist[i][0] = refts_pos[i][0];
      refTS_xlist[i][1] = refts_pos[i][1];
      refTS_xlist[i][2] = refts_pos[i][2];

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

  movR_to_ref[0] = refreactant_com[0] - mov_com[0];
  movR_to_ref[1] = refreactant_com[1] - mov_com[1];
  movR_to_ref[2] = refreactant_com[2] - mov_com[2];

  movP_to_ref[0] = refproduct_com[0] - mov_com[0];
  movP_to_ref[1] = refproduct_com[1] - mov_com[1];
  movP_to_ref[2] = refproduct_com[2] - mov_com[2];

  movTS_to_ref[0] = refts_com[0] - mov_com[0];
  movTS_to_ref[1] = refts_com[1] - mov_com[1];
  movTS_to_ref[2] = refts_com[2] - mov_com[2];

  // Here for debugging purposes, write oordinates of mov-xlist are indeed those of the alignement atoms
  int step = getStep();

  ofstream wfile000;
  wfile000.open("mov-xlist-bef.xyz");
  int n_align = reactant_atoms.size();
  wfile000 << n_align << endl;
  wfile000 << "comment" << endl;
  for (int i=0; i < n_align; i++)
    {
      wfile000<< params.R_atomName[i] << " " << std::fixed << std::setprecision(5) << mov_xlist[i][0]*10 << " " << mov_xlist[i][1]*10 << " " << mov_xlist[i][2]*10 << endl;
      //wfile000<< "C " << std::fixed << std::setprecision(5) << mov_xlist[i][0]*10 << " " << mov_xlist[i][1]*10 << " " << mov_xlist[i][2]*10 << endl;
    }
  wfile000.close();


  rotmatR[0][0] = 1.0;
  rotmatR[0][1] = 0.0;
  rotmatR[0][2] = 0.0;
  rotmatR[1][0] = 0.0;
  rotmatR[1][1] = 1.0;
  rotmatR[1][2] = 0.0;
  rotmatR[2][0] = 0.0;
  rotmatR[2][1] = 0.0;
  rotmatR[2][2] = 1.0;

  rotmatP[0][0] = 1.0;
  rotmatP[0][1] = 0.0;
  rotmatP[0][2] = 0.0;
  rotmatP[1][0] = 0.0;
  rotmatP[1][1] = 1.0;
  rotmatP[1][2] = 0.0;
  rotmatP[2][0] = 0.0;
  rotmatP[2][1] = 0.0;
  rotmatP[2][2] = 1.0;

  rotmatTS[0][0] = 1.0;
  rotmatTS[0][1] = 0.0;
  rotmatTS[0][2] = 0.0;
  rotmatTS[1][0] = 0.0;
  rotmatTS[1][1] = 1.0;
  rotmatTS[1][2] = 0.0;
  rotmatTS[2][0] = 0.0;
  rotmatTS[2][1] = 0.0;
  rotmatTS[2][2] = 1.0;

  calculate_rotation_rmsd( refR_xlist, mov_xlist, n_reactantatoms, mov_com, movR_to_ref, rotmatR, &rmsdR  );
  calculate_rotation_rmsd( refP_xlist, mov_xlist, n_reactantatoms, mov_com, movP_to_ref, rotmatP, &rmsdP  );
  calculate_rotation_rmsd( refTS_xlist, mov_xlist, n_reactantatoms, mov_com, movTS_to_ref, rotmatTS, &rmsdTS  );
  
  string  xyzfile = "mov-xlist.xyz";
  ostringstream convertstep;
  convertstep << step;
  string xyzstep;
  xyzstep = convertstep.str() + xyzfile;
  // cout << xyzstep << endl;
  ofstream wfile001;
  wfile001.open(xyzstep.c_str());
  wfile001 << n_align << endl;
  wfile001 << "comment" << endl;
  for (int i=0; i < n_align; i++)
    {
      wfile001<< params.R_atomName[i] << " " << std::fixed << std::setprecision(5) << mov_xlist[i][0]*10 << " " << mov_xlist[i][1]*10 << " " << mov_xlist[i][2]*10 << endl;
    }
  wfile001.close();

  ofstream wfile002;
  wfile002.open("refreactant-xlist.xyz");
  wfile002 << n_align << endl;
  wfile002 << "comment" << endl;
  for (int i=0; i < n_align; i++)
    {
      wfile002<< params.R_atomName[i] << " " << std::fixed << std::setprecision(5) << refR_xlist[i][0]*10 << " " << refR_xlist[i][1]*10 << " " << refR_xlist[i][2]*10 << endl;
    }
  wfile002.close();

  ofstream wfile003;
  wfile003.open("refproduct-xlist.xyz");
  wfile003 << n_align << endl;
  wfile003 << "comment" << endl;
  for (int i=0; i < n_align; i++)
    {
      wfile003<< params.P_atomName[i] << " " << std::fixed << std::setprecision(5) << refP_xlist[i][0]*10 << " " << refP_xlist[i][1]*10 << " " << refP_xlist[i][2]*10 << endl;
    }
  wfile003.close();

  ofstream wfile004;
  wfile004.open("refts-xlist.xyz");
  wfile004 << n_align << endl;
  wfile004 << "comment" << endl;
  for (int i=0; i < n_align; i++)
    {
      wfile004<< params.TS_atomName[i] << " " << std::fixed << std::setprecision(5) << refTS_xlist[i][0]*10 << " " << refTS_xlist[i][1]*10 << " " << refTS_xlist[i][2]*10 << endl;
    }
  wfile004.close();


/**
  cout << "Just called calculated_rotation_rmsd" << endl;
  cout << "rotmat elements  of reactant reference :" << endl;
  cout << rotmatR[0][0] << " " << rotmatR[0][1] << " " << rotmatR[0][2] << endl;
  cout << rotmatR[1][0] << " " << rotmatR[1][1] << " " << rotmatR[1][2] << endl;
  cout << rotmatR[2][0] << " " << rotmatR[2][1] << " " << rotmatR[2][2] << endl; 
  cout << "rmsd " << rmsdR << endl;

  cout << "Just called calculated_rotation_rmsd" << endl;
  cout << "rotmat elements  of product reference :" << endl;
  cout << rotmatP[0][0] << " " << rotmatP[0][1] << " " << rotmatP[0][2] << endl;
  cout << rotmatP[1][0] << " " << rotmatP[1][1] << " " << rotmatP[1][2] << endl;
  cout << rotmatP[2][0] << " " << rotmatP[2][1] << " " << rotmatP[2][2] << endl; 
  cout << "rmsd " << rmsdP << endl;

  cout << "Just called calculated_rotation_rmsd" << endl;
  cout << "rotmat elements  of ts reference :" << endl;
  cout << rotmatTS[0][0] << " " << rotmatTS[0][1] << " " << rotmatTS[0][2] << endl;
  cout << rotmatTS[1][0] << " " << rotmatTS[1][1] << " " << rotmatTS[1][2] << endl;
  cout << rotmatTS[2][0] << " " << rotmatTS[2][1] << " " << rotmatTS[2][2] << endl; 
  cout << "rmsd " << rmsdTS << endl;
*/
// calculate energy
  double Coulomb_R = 0.0; // the Coulombic energy of reactant with binding site
  double Coulomb_P = 0.0; // the Coulombic energy of product with binding site
  double Coulomb_TS = 0.0; // the Coulombic energy of transition state with binding site
  double LJ_R = 0.0; // the LJ energy of reactant with binding site
  double LJ_P = 0.0; // the LJ energy of product with binding site
  double LJ_TS = 0.0; // the LJ energy of transition state with binding site
  double Ke = 8.98755e9;
  double qi = 0.0;
  double qjR = 0.0;
  double qjP = 0.0;
  double qjTS = 0.0;
  unsigned n_bindingatoms = binding_atoms.size();
  /// loop over binding site atom
  for (unsigned i= n_reactantatoms; i < n_reactantatoms+n_bindingatoms; ++i)
    {
     ///loop over reactant atom
     Vector binding_pos = getPosition(i); // coordinate of binding site atom i
     qi = params.bindSite_charge[i];
     for (unsigned j=0; j < n_reactantatoms; ++j )
       { 
                 

         // calculate distance atom i and j (magnitude of Rij)
         double Rx = binding_pos[0] - refR_xlist[j][0] ;
         double Ry = binding_pos[1] - refR_xlist[j][1] ;
         double Rz = binding_pos[2] - refR_xlist[j][2] ;
         double distR = sqrt(Rx*Rx + Ry*Ry + Rz*Rz) ;
         // distance atom ref product and binding site
         double Px = binding_pos[0] - refP_xlist[j][0] ;
         double Py = binding_pos[1] - refP_xlist[j][1] ;
         double Pz = binding_pos[2] - refP_xlist[j][2] ;
         double distP = sqrt(Px*Px + Py*Py + Pz*Pz) ;
       //  cout << "Px " << Px << "Py " << Py << "Pz " << Pz << "dist "<< distP << endl;
         // distance atom ref ts and binding site
         double TSx = binding_pos[0] - refTS_xlist[j][0] ;
         double TSy = binding_pos[1] - refTS_xlist[j][1] ;
         double TSz = binding_pos[2] - refTS_xlist[j][2] ;
         double distTS = sqrt(TSx*TSx + TSy*TSy + TSz*TSz) ;
      //   cout << "TSx " << TSx << "TSy " << TSy << "TSz " << TSz << "dist "<< distTS << endl;
         qjR = params.R_charge[j];     //qj = getCharge(j);
         qjP = params.P_charge[j];     //qj = getCharge(j);
         qjTS = params.TS_charge[j];     //qj = getCharge(j);
     // Coulomb_R += ((Ke*qi*qj)/distR) ; // the Coulombic energy of reactant with binding site
     
         Coulomb_P += ((Ke*qi*qjP)/distP) ; // the Coulombic energy of product with binding site
     
         Coulomb_TS += ((Ke*qi*qjP)/distTS) ; // the Coulombic energy of TS with binding site
        
         Coulomb_R += ((Ke*qi*qjTS)/distR) ; // the Coulombic energy of reactant with binding site
        // cout << "dist R " << distR << " distP " << distP <<  " distTS " << distTS << endl; 
         // calculate the Lennard Jones energy
         // V(LJ) = 4(Eij)(sigmaij^12/(Rij)^12 - sigmaij^6/(Rij)^6)
         if (LJ_rule==1)
         {
           double R_sigmaij = 0.5*(params.bindSite_sigma[i-n_reactantatoms] + params.R_sigma[j]);
           double R_epsilonij = sqrt(params.bindSite_epsilon[i-n_reactantatoms]*params.R_epsilon[j]);
           LJ_R += 4*R_epsilonij*(  pow((R_sigmaij/distR),12) - pow( (R_sigmaij/distR),6)  ) ;
       
           //cout << "bind atom " << i-n_reactantatoms <<" sigma " << params.bindSite_sigma[i-n_reactantatoms] << " R atom " << j << " sigma " << params.R_sigma[j] << " sigmaij " << R_sigmaij << endl;
           //cout << "bind atom " << i-n_reactantatoms <<" epsilon " << params.bindSite_epsilon[i-n_reactantatoms] << " R atom " << j << " epsilon " << params.R_epsilon[j] << " epsilonij " << R_epsilonij << endl;
           //cout << "bind atom " << i-n_reactantatoms << " R atom " << j << " sigma ij  " << R_sigmaij << " epsilon ij " <<   R_epsilonij << " dist R " << distR << " LJ energy " << LJ_R << endl;
           double P_sigmaij = 0.5*(params.bindSite_sigma[i-n_reactantatoms] + params.P_sigma[j]);
           double P_epsilonij = sqrt(params.bindSite_epsilon[i-n_reactantatoms]*params.P_epsilon[j]);
           LJ_P += 4*P_epsilonij*(  pow((P_sigmaij/distP),12) - pow( (P_sigmaij/distP),6)  ) ;
      
           double TS_sigmaij = 0.5*(params.bindSite_sigma[i-n_reactantatoms] + params.TS_sigma[j]);
           double TS_epsilonij = sqrt(params.bindSite_epsilon[i-n_reactantatoms]*params.TS_epsilon[j]);
           LJ_TS += 4*TS_epsilonij*(  pow((TS_sigmaij/distTS),12) - pow( (TS_sigmaij/distTS),6)  ) ;
         }
         else
         {
           // VLJ = Aij/Rij^12 - Bij/R^6
           double Ai = 4*params.bindSite_epsilon[i-n_reactantatoms]* ( pow( params.bindSite_sigma[i-n_reactantatoms],12) ) ;
           double Bi = 4*params.bindSite_epsilon[i-n_reactantatoms]* ( pow( params.bindSite_sigma[i-n_reactantatoms],6) ) ;
           double Aj_R = 4*params.R_epsilon[j]* ( pow( params.R_sigma[j],12) ) ;
           double Bj_R = 4*params.R_epsilon[j]* ( pow( params.R_sigma[j],6) ) ;
	   double Aij_R = sqrt( Ai * Aj_R ) ;
	   double Bij_R = sqrt( Bi * Bj_R ) ;

           LJ_R += (Aij_R/pow(distR,12)) - (Bij_R/pow(distR,6))  ;
       
           //cout << "bind atom " << i-n_reactantatoms <<" sigma " << params.bindSite_sigma[i-n_reactantatoms] << " R atom " << j << " sigma " << params.R_sigma[j] << " sigmaij " << R_sigmaij << endl;
           //cout << "bind atom " << i-n_reactantatoms <<" epsilon " << params.bindSite_epsilon[i-n_reactantatoms] << " R atom " << j << " epsilon " << params.R_epsilon[j] << " epsilonij " << R_epsilonij << endl;
           //cout << "bind atom " << i-n_reactantatoms << " R atom " << j << " sigma ij  " << R_sigmaij << " epsilon ij " <<   R_epsilonij << " dist R " << distR << " LJ energy " << LJ_R << endl;
           double Aj_P = 4*params.P_epsilon[j]* ( pow( params.P_sigma[j],12) ) ;
           double Bj_P = 4*params.P_epsilon[j]* ( pow( params.P_sigma[j],6) ) ;
	   double Aij_P = sqrt( Ai * Aj_P ) ;
	   double Bij_P = sqrt( Bi * Bj_P ) ;

           LJ_P += (Aij_P/pow(distP,12)) - (Bij_P/pow(distP,6))  ;

           double Aj_TS = 4*params.TS_epsilon[j]* ( pow( params.TS_sigma[j],12) ) ;
           double Bj_TS = 4*params.TS_epsilon[j]* ( pow( params.TS_sigma[j],6) ) ;
	   double Aij_TS = sqrt( Ai * Aj_TS ) ;
	   double Bij_TS = sqrt( Bi * Bj_TS ) ;

           LJ_TS += (Aij_TS/pow(distTS,12)) - (Bij_TS/pow(distTS,6))  ;
           cout << " LJ defalt " << endl;    
         } 

       }
    }
    cout << "Step " << step << endl;
    cout << " CoulombR " << Coulomb_R << " LJ R " << LJ_R << endl;
    cout << "CoulombP " << Coulomb_P << "LJ P " << LJ_P <<  endl;
    cout << "CoulombTS " << Coulomb_TS <<  " LJ TS " << LJ_TS << endl;
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

