/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "tools/PDB.h"

#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <iostream>

// Kabsch algorithm implementation
#include "kabsch.h"

typedef double real;

using namespace std;

namespace PLMD
{
namespace colvar
{

vector<double> set_bs_values(vector<Vector> grid_positions,
			     vector<Vector> site_positions,
			     double theta, double Bsmin, double deltaBS);

vector<vector<int> > init_grid_neighbors(vector<Vector> grid_positions,
					 double GP1_min, double deltaGP1,
					 double GP2_min, double deltaGP2);

double s_on(double k,double v,double v_min,double delta);
double s_off(double k,double v,double v_min,double delta);
double ds_off_dm(double k,double v,double v_min,double delta);
double ds_off_dk(double k,double v,double v_min,double delta);
double ds_on_dm(double k,double v,double v_min,double delta);
double ds_on_dk(double k,double v,double v_min,double delta);

//+PLUMEDOC COLVAR JEDI
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

 JEDI ATOM=<atom selection> APOLAR=<atom selection> POLAR=<atom selection> COM_x=     COM_y=     COM_z=     B_grid=     N_grid=     CUTOFF_close=0.22  CUTOFF_far=0.14  CUTOFF_enclosure=3  CUTOFF_hull=3  CUTOFF_surface=0.21  CUTOFF_contact=1  CUTOFF_hydro=0.4  CUTOFF_con=5  dump=     SIGMA=0.05 

*/
//+ENDPLUMEDOC

// Internal class used to store parameter set
class jediparameters
{
public:
  jediparameters();
  bool readParams(string &parameters_file);
  double alpha;
  double beta;
  double gamma;
  double theta;
  double CC_mind;
  double deltaCC;
  double Emin;
  double deltaE;
  double BSmin;
  double deltaBS;
  double CC2_min;
  double deltaCC2;
  double GP1_min;
  double deltaGP1;
  double GP2_min;
  double deltaGP2;
  double r_hydro;
  double deltar_hydro;
  double V_max;
  double deltaV_max;
  double V_min;
  double deltaV_min;
};

// Constructor
jediparameters::jediparameters()
{
  alpha = 0.0;
  beta = 0.0;
  gamma = 0.0;
  theta = 0.0;
  CC_mind = 0.0;
  deltaCC = 0.0;
  Emin = 0.0;
  deltaE = 0.0;
  BSmin = 0.0;
  deltaBS = 0.0;
  CC2_min = 0.0;
  deltaCC2 = 0.0;
  GP1_min = 0.0;
  deltaGP1 = 0.0;
  GP2_min = 0.0;
  deltaGP2 = 0.0;
  r_hydro = 0.0;
  deltar_hydro = 0.0;
  V_max = 0.0;
  deltaV_max = 0.0;
  V_min = 0.0;
  deltaV_min = 0.0;
}

bool jediparameters::readParams(string &parameters_file)
{
  FILE* fp=fopen(parameters_file.c_str(),"r");
  if (!fp) return false;
  string line;
  while (Tools::getline(fp, line))
    {
      if (line[0] == '#')
	continue;
      //cout << line << endl;
      istringstream iss(line);
      istream_iterator<string> beg(iss), end;
      vector<string> tokens(beg, end); // done!
      if (tokens.size() < 2)
	{
	  cout << "ABORT ! JEDI parameters file appear malformed check the line below !" << endl; 
	  cout << line << endl;
	  exit(-1);
	}
      string key = tokens[0];
      double item = atof(tokens[2].c_str());
      //cout << " key " << key << " " << " item " << item << endl;
      // set jedi parameters datastruct
      // FIXME: Check that all parameters have been set
      if ( key == string("alpha") )
	alpha = item;
      else if ( key == string("beta") )
	beta = item;
      else if ( key == string("gamma") )
	gamma = item;
      else if ( key == string("theta") )
	theta = item;
      else if ( key == string("CC_mind") )
	CC_mind = item;
      else if ( key == string("deltaCC") )
	deltaCC = item;
      else if ( key == string("Emin") )
	Emin = item;
      else if ( key == string("deltaE") )
	deltaE = item;
      else if ( key == string("BSmin") )
	BSmin = item;
      else if ( key == string("deltaBS") )
	deltaBS = item;
      else if ( key == string("CC2_min") )
	CC2_min = item;
      else if ( key == string("deltaCC2") )
	deltaCC2 = item;
      else if ( key == string("GP1_min") )
	GP1_min = item;
      else if ( key == string("deltaGP1") )
	deltaGP1 = item;
      else if ( key == string("GP2_min") )
	GP2_min = item;
      else if ( key == string("deltaGP2") )
	deltaGP2 = item;
      else if ( key == string("r_hydro") )
	r_hydro = item;
      else if ( key == string("deltar_hydro") )
	deltar_hydro = item;
      else if ( key == string("V_max") )
	V_max = item;
      else if ( key == string("deltaV_max") )
	deltaV_max = item;
      else if ( key == string("V_min") )
	V_min = item;
      else if ( key == string("deltaV_min") )
	deltaV_min = item;
    }
  fclose(fp);

  cout << "*** Values of the JEDI parameter set loaded in memory: ***" << endl;
  cout << "alpha = " << alpha << endl;
  cout << "beta  = " << beta << endl;
  cout << "gamma = " << gamma << endl;
  cout << "theta = " << gamma << endl;  
  cout << "CC_mind  = " << CC_mind << endl;
  cout << "deltaCC  = " << deltaCC << endl;
  cout << "Emin  = " << Emin << endl;
  cout << "deltaE  = " << deltaE << endl;
  cout << "BSmin  = " << BSmin << endl;
  cout << "deltaBS  = " << deltaBS << endl;
  cout << "CC2_min  = " << CC2_min << endl;
  cout << "deltaCC2  = " << deltaCC2 << endl;
  cout << "GP1_min  = " << GP1_min << endl;
  cout << "deltaGP1  = " << deltaGP1 << endl;
  cout << "GP2_min  = " << GP2_min << endl;
  cout << "deltaGP2  = " << deltaGP2 << endl;
  cout << "r_hydro  = " << r_hydro << endl;
  cout << "deltar_hydro  = " << deltar_hydro << endl;
  cout << "V_max  = " << V_max << endl;
  cout << "deltaV_max  = " << deltaV_max << endl;
  cout << "V_min  = " << V_min << endl;
  cout << "deltaV_min  = " << deltaV_min << endl;

  return true;
}



class jedi : public Colvar
{
private:
  bool pbc;
  //for JEDI
  vector<AtomNumber> alignmentatoms;//list of atoms used for alignments
  vector<AtomNumber> apolaratoms;//list of apolar atoms used for CV
  vector<AtomNumber> polaratoms;//list of polar atoms used for CV
  vector<Vector> ref_pos;// coordinates reference structure for alignment.
  double ref_com[3];// coordinates of the center of mass of the reference structure for alignments
  vector<Vector> grid_positions;//coordinates of the reference grid for alignment
  vector<double> grid_s_off_bsi;//binding site score of grid point (eq 5 term 1 Cuchillo et al. JCTC 2015)
  jediparameters params;// parameters druggability estimator
  vector<vector<int> > neighbors;//list of grid indices that are neighbors of a given grid point
  string summary_file;//path to output file
  int stride;//frequency of output (in timesteps) to summary file;
  string gridstats_folder;//path to output grid folder;
  int gridstride;//frequency of output (in timestep) of a grid file;
  //DEPRECATED
  vector<vector<int> > rays;
  //vector<vector<int> > neighbors;//DEPRECATED NOT USED??
  vector<vector<double> > grid_pos;//DEPRECATED replaced by grid_positions
  vector<double> pocket;//DEPRECATED, replaced by grid_s_off_bsi
  //float score[9];//rotation matrix elements
  vector<AtomNumber> Apolar;//for plumed.dat file DEPRECATED replaced by apolaratoms
  vector<AtomNumber> Polar;// for plumed.dat file DEPRECATED replaced by polaratoms

  double site_com[3];//reference coordinates of the center of mass of the binding site region
  double COM_X, COM_Y, COM_Z;// coordinates of the center of mass of the binding site region

  double b_grid;//atom number of the first grid point
  double n_grid;// total number of grid points (double)
  double cutoff_close;//cutoff close contact
  double cutoff_far;//cutoff distant contact
  double cutoff_enclosure;//cutoff NP for distant contact
  double cutoff_hull;//cutoff hull
  double cutoff_surface;//cutoff surface
  double cutoff_contact;//cutoff contact for surface
  double cutoff_hydro;//cutoff hydrophobicity
  double cutoff_con;//cutoff connectivity
  double dump_matrix;//paramater to control the output of jedi (just to check, could be set to 0)
  double delta;
  string line;
public:
 jedi(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(jedi,"JEDI")

void jedi::registerKeywords(Keywords& keys)
{
  Colvar::registerKeywords(keys);
  //  keys.addFlag("JEDI_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
  //  keys.addFlag("JEDI_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  keys.add("compulsory","SIGMA","0.05","Gaussian width for metadynamics calculations");
  keys.add("compulsory","REFERENCE","a file in pdb format containing the atoms to use for computing rotation matrices.");
  keys.add("compulsory","APOLAR","a file in pdb format containing the apolar protein atoms to use for the CV.");
  keys.add("compulsory","POLAR","a file in pdb format containing the polar protein atoms to use for the CV.");
  keys.add("compulsory","GRID","a file in pdb format containing the grid points to use for the CV.");
  keys.add("compulsory","PARAMETERS","a file listing the parameters of the JEDI estimator.");
  keys.add("compulsory", "SITE","a file listing coordinates of atoms used to define a binding site region.");
  keys.add("compulsory","STRIDE","100","frequency of output to jedi summary file.");
  keys.add("compulsory","SUMMARY","jedi_stats.dat","summary file jedi descriptor.");
  keys.add("compulsory","GRIDSTRIDE","100","frequency of output of jedi grid.");
  keys.add("compulsory","GRIDFOLDER","jedi-grids", "folder where jedi grids will be output.");
}

jedi::jedi(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  parse("SIGMA", delta);//FIXME: Where is this set//used?? Plumed convention??;
  string reference_file;
  parse("REFERENCE",reference_file);
  string apolar_file;
  parse("APOLAR", apolar_file);
  string polar_file;
  parse("POLAR", polar_file);
  string grid_file;
  parse("GRID", grid_file);
  string parameters_file;
  parse("PARAMETERS", parameters_file);
  //FIXME: Parse error if site no provided
  string site_file;
  parse("SITE", site_file);
  if (site_file.length() == 0)
    site_file = "null";
  string stride_string;
  parse("STRIDE",stride_string);
  int stride=atoi(stride_string.c_str());
  string summary_file;
  parse("SUMMARY",summary_file);
  string gridstride_string;
  parse("GRIDSTRIDE",gridstride_string);
  int gridstride=atoi(gridstride_string.c_str());
  string gridstats_folder;
  parse("GRIDFOLDER",gridstats_folder);


  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc)
    log.printf("  using periodic boundary conditions\n");
  else 
    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();//FIXME WHAT TO DO PBC?
 
  cout << "*** Initialisation of JEDI collective variable ***" << endl;

  //Apolar
  vector<AtomNumber> Apolar;
  //cout << " Apolar has ? elements " << Apolar.size() << endl;

  PDB apolar_pdb;
  if( !apolar_pdb.read(apolar_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + apolar_file );

  apolaratoms = apolar_pdb.getAtomNumbers();
  cout << " apolaratoms has ? elements " << apolaratoms.size() << endl;  

  //Polar 
  vector<AtomNumber> Polar;
  //cout << " Polar has ? elements " << Polar.size() << endl;

  PDB polar_pdb;
  if( !polar_pdb.read(polar_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + polar_file );
  polaratoms = polar_pdb.getAtomNumbers();
  cout << " polaratoms has ? elements " << polaratoms.size() << endl;

  // Load up alignment file
  PDB reference_pdb;
  // read everything in ang and transform to nm if we are not in natural units
  if( !reference_pdb.read(reference_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + reference_file );
  // Save in memory the reference coordinates of the atoms to use for future alignments
  const std::vector<Vector> alignment_positions = reference_pdb.getPositions();
  // Masses are taken from the occupancy array
  const std::vector<double> alignment_masses = reference_pdb.getOccupancy();
  // Also compute the com of reference coordinates and save for future calcs
  double ref_mass_tot=0.0;
  ref_com[0] = 0.0;
  ref_com[1] = 0.0;
  ref_com[2] = 0.0;
  for (unsigned i=0; i < alignment_positions.size() ; ++i)
    {
      ref_com[0] += alignment_masses[i] * alignment_positions[i][0];
      ref_com[1] += alignment_masses[i] * alignment_positions[i][1];
      ref_com[2] += alignment_masses[i] * alignment_positions[i][2];
      ref_pos.push_back( Vector(alignment_positions[i][0], alignment_positions[i][1], alignment_positions[i][2]) );
      ref_mass_tot += alignment_masses[i];
    }
  ref_com[0] /= ref_mass_tot;
  ref_com[1] /= ref_mass_tot;
  ref_com[2] /= ref_mass_tot;

  // Add the alignment atom numbers to the list of atoms to request from the CV
  alignmentatoms = reference_pdb.getAtomNumbers();
  cout << " alignmentatoms has ? elements " << alignmentatoms.size() << endl;    
  //  for(unsigned i=0; i < alignmentatoms.size();++i)
  //  cout << " i " << i << " alignmentatoms number " << alignmentatoms[i].serial() << endl;

  //vector<AtomNumber> allatoms( Apolar.size() + Polar.size() + alignmentatoms.size() );
  vector<AtomNumber> allatoms( apolaratoms.size() + polaratoms.size() + alignmentatoms.size() );

  for ( unsigned i = 0; i < apolaratoms.size() ; ++i)
    allatoms[i]=apolaratoms[i];

   for ( unsigned i = 0; i < polaratoms.size() ; ++i)
    allatoms[apolaratoms.size()+i] = polaratoms[i];

  for ( unsigned i=0; i < alignmentatoms.size() ; ++i)
    allatoms[ apolaratoms.size() + polaratoms.size() + i ] = alignmentatoms[i];

  requestAtoms(allatoms);

  // FIXME set reference COM of binding site region
  // THIS CAN BE DONE BY DOING COM CALCULATION OVER DIFFERENT SET
  // OF ATOMS (POLAR+APOLAR)

  //READ jedi.parameters here
  params.readParams(parameters_file);

  // (Optional) If specified, load up site file
  PDB site_pdb;
  std::vector<Vector> site_positions;
  if ( site_file != string("null") )
    {
      if( !site_pdb.read(site_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
	error("missing input file " + site_file );
      site_positions = site_pdb.getPositions();
    }
  cout << " site_positions has ? elements " << site_positions.size() << endl;

  // Load up grid file
  PDB grid_pdb;
  // read everything in ang and transform to nm if we are not in natural units
  if( !grid_pdb.read(grid_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + grid_file );
  // Save in memory the reference coordinates of the grid points for future alignments
  //const std::vector<Vector>
  grid_positions = grid_pdb.getPositions();
  cout << " grid_positions has ? elements " << grid_positions.size() << endl;
  // Set maximum activity of grid points
  grid_s_off_bsi = set_bs_values(grid_positions, site_positions, params.theta, params.BSmin, params.deltaBS);

  //const std::vector<AtomNumber> gridnumbers = grid_pdb.getAtomNumbers();
  //cout << " grid has ? elements " << gridnumbers.size() << endl;
  neighbors = init_grid_neighbors( grid_positions,
				   params.GP1_min, params.deltaGP1,
				   params.GP2_min, params.deltaGP2);
  //Setup output
  // JM TODO: Check behavior CV init upon job restart.
  ofstream wfile;
  wfile.open(summary_file.c_str());
  wfile << "#step \t JEDI \t Va \t Ha \t Va/Vmax \t COM_x \t COM-y \t COM_z \t RotMat[0][0].Rotmat[0][1]....Rotmat[2][2]" << endl;
  wfile.close();

  //TODO CHECK IF FOLDER GRIDSTATS_FOLDER EXISTS
  // IF YES DELETE
  // CREATE NEW EMPTY FOLDER

  cout << "*** Completed initialisation JEDI collective variable" << endl;
  exit(0);
}

vector<double> set_bs_values( vector<Vector> grid_pos,
			      vector<Vector> site_pos,
			      double theta, double BSmin, double deltaBS)
{
  vector<double> grid_s_off_bsi;

  for (int i=0; i < grid_pos.size(); i++)
    {
      grid_s_off_bsi.push_back(1.0);
    }
  
  if (site_pos.size() > 0)
    {
      for (int i=0; i < grid_pos.size(); i++)
	{
	  double sum=0.0;
	  for (int j=0; j < site_pos.size(); j++)
	    {
	      // Precompute first term eq5 Cuchillo et al. JCTC 2015
	      double dij2 = pow(grid_pos[i][0] - site_pos[j][0],2) +
		pow(grid_pos[i][1] - site_pos[j][1],2) +
		pow(grid_pos[i][2] - site_pos[j][2],2);
	      double dij = sqrt(dij2);
	      sum = sum + exp(theta/dij);
	    };
	  double bsi = theta/log(sum);
	  //grid_bsi.push_back(bsi);
	  double ai = s_off(1.0, bsi, BSmin, deltaBS);
	  //cout << " i " << i << " bsi " << bsi << " ai " << ai << endl;
	  grid_s_off_bsi[i] = ai;
	}
    }
  return grid_s_off_bsi;
}

vector<vector<int> > init_grid_neighbors(vector<Vector> grid_pos,
					 double GP1_min, double deltaGP1,
					 double GP2_min, double deltaGP2)
{
  size_t grid_size=grid_pos.size();
  vector<vector<int> > neighbors;
  neighbors.reserve(grid_size);

  for (int i=0; i < grid_pos.size(); i++)
    {
      vector<int> list;
      neighbors.push_back(list);
    }

  double dmin = pow(GP1_min,2);
  double dmax = pow(GP2_min+deltaGP2,2);

  for (int i=0; i < grid_pos.size(); i++)
    {
      for (int j=i+1; j < grid_pos.size(); j++)
	{
	  //Prepare calculation of eq 7 Cuchillo et atl. JCTC 2015
	  double dij2 = pow(grid_pos[i][0] - grid_pos[j][0],2) +
	    pow(grid_pos[i][1] - grid_pos[j][1],2) +
	    pow(grid_pos[i][2] - grid_pos[j][2],2);
	  if (dij2 > dmin and dij2 < dmax)
	    {
	      //	      cout << " i " << i << " j " << j << " dij2 " << dij2 << " dmin " << dmin << " dmax " << dmax << endl;
	      neighbors[i].push_back(j);
	      neighbors[j].push_back(i);
	    }
	}
      //cout << i << " : ";
      //for (int k=0;k< neighbors[0].size();k++)
      //	cout << neighbors[i][k] << " ";
      //cout << endl;
      //exit(0);
    }
  // JM comment 09/15. Seem to get at most 38 neighbors but paper said 44. Parameter set?
  //cout << " Done populating neighbors list" << endl;
  //exit(0);
  //for (int i=0; i < grid_pos.size(); i++)
  //  {
  //    cout << " i " << i << " : ";
  //    cout << " has ? elements " << neighbors[i].size() << " : " ;
  //    for (int j=0; j < neighbors[i].size(); j++)
  //    	{
  //    	  cout << neighbors[i][j] << " ";
  //     	}
  //    cout << endl;
  //  }

  return neighbors;
}

///////////////////////////////////////////////////////////////////////////////////
//////// This is the routine where JEDI and its derivatives are calculated/////////
///////////////////////////////////////////////////////////////////////////////////

// subroutines to calculate the switching functions and their derivatives
 double s_on(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
    s=0.;
    }
  else if (m > 1.)
    {
    s=k;
    }
  else
    {
      s=k*(1.-(pow((1.-pow(m,2)),2))*((1.+2*(pow(m,2)))));
    }
  return s;
 }

 double s_off(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
    s=k;
    }
  else if (m > 1.)
    {
    s=0.;
    }
  else
    {
      s=k*(pow((1.-pow(m,2)),2)*((1.+2*(pow(m,2)))));
    }
  return s;
 }

 double ds_off_dm(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
      s = 0.;
    }
  else if (m > 1.)
    {
      s = 0.;
    }
  else
    {
      s = 4*k*m*(pow((1.-pow(m,2)),2)) - 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
    }
  return s;
 }

 double ds_off_dk(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
      s = 1.;
    }
  else if (m > 1.)
    {
      s = 0.;
    }
  else
    {
      s = (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2))));
    }
  return s;
 }

 double ds_on_dm(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
      s = 0.;
    }
  else if (m > 1.)
    {
      s = 0.;
    }
  else
    {
      s = -4*k*m*(pow((1.-pow(m,2)),2)) + 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
    }
  return s;
 }

 double ds_on_dk(double k,double v,double v_min,double delta)
 {
  double s, m;
  m = (v-v_min)/(delta);
  if (m < 0.)
    {
      s = 0.;
    }
  else if (m > 1.)
    {
      s = 1.;
    }
  else
    {
      s = 1. - ( (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2)))) );
    }
  return s;
 }

// calculator
void jedi::calculate(){

  int i, j, jj, k, l, gridpoint;
  //int N, protein, gridpoint1, gridpoint2;
  double Jedi, volume, a, b, hydrophobicity, constant, grd_x, grd_y, grd_z;
  //double enclosure,c,d, connectivity;
  double resolution, Vg, mod_rij, mod_rjk, mod_rik, pe, sum, min_dist;// penalty_close
  double penalty_far;
  double ncoord, D_far;
  double D_close;// D_enclosure;
  double s1, s2, current, time1, time2, Va, Vmax, Vmin, D_Vmax, D_Vmin;//, prev_snap;
  //double s3, s4, dump_time;
  //double connectivity_tot; 
  //double ncontact
  double apolar, polar, D_hydro, hydrophobicity_tot;
  //double surface, hull, D_hull, D_surface, D_contact;
  //double ncontact_hull, ncontact_surface;
  
  int size_grid;
  size_grid =  n_grid; // total number of grid points !! for loops and array !!(integer)
  
  //unsigned int size_protein;
  unsigned int size_apolar;
  unsigned int size_polar;
  //unsigned int size_connect_map;
   // total number of grid points !! for loops and array !!(integer)
  vector<int> active_grid;// array with the index of active grid point (0 if inactive)
  double penalty_list[size_grid];
  double enclosure_list[size_grid];// array with the penalty of distant contact
  double mind_list[size_grid];// array with the penalty of close contact
  vector<double> s3_m;// array with the derivative of distant contact
  vector<double> s4_m; // array with the derivative of distant contact
  double sum_list[size_grid];// array with the sum of the distances between each grid point and all atoms in the CV
  double min_dist_list[size_grid];// array with the minimum distance between each grid point and all atoms in the CV
  //double NP[size_grid];// array with the number of protein atoms within around each grid points according to s_on
  double new_x[size_grid];// array with update of x coordinates of each grid points according to translation/rotation
  double new_y[size_grid];// array with update of y coordinates of each grid points according to translation/rotation
  double new_z[size_grid];// array with update of z coordinates of each grid points according to translation/rotation
  vector<double> connectivity_list; // array with the number of active grid points around each grid points according to s_on
  vector<double> hull_list;// array with the hull score of each grid points
  vector<double> surface_list;// array with the surface score of each grid points
  vector<double> ncontact_surface_list;// array number of protein atoms in contact with each grid point (for derivatives)
  double hydrophobicity_list[size_grid];
  double apolarity[size_grid];// array number of apolar protein atoms around each grid point
  double polarity[size_grid];// array number of polar protein atoms around each grid point
  double ray[size_grid];
  
  
  volume      = 0.;// score volume : total number of active and partially active grid points
  //N           = 0;// to modify
  D_close     = 0.05;// delta close contact
  D_far       = 0.05;// delta distant contact
  //D_enclosure = 5.0;// delta NP contact ( for distant contact)
  resolution  = 0.15;// grid spacing
  s1          = 0.0;// minimum volume according our drug-like dataset
  s2          = 0.0;// maximum volume according our drug-like dataset
  Va          = 0.0;// penalty non drug-like pockets ( s1 * s2 )
  Vmin        = cutoff_hull;// minimum volume if volume < Vmin --> s1 = 0
  D_Vmin      = 10.0;// delta minimum volume volume > Vmin + D_Vmin --> s1 = 1
  Vmax        = 90.0;// minimum volume if volume < Vmax --> s2 = 1
  D_Vmax      = 50.0;// delta maximum volume volume > Vmax + D_Vmin --> s1 = 0
  //enclosure   = 0.;// score enclosure : surface / hull
  //hull        = 0.;// score hull
  //surface     = 0.;// score surface
  //D_hull      = 3;// delta hull
  //D_surface   = 0.05;// delta surface
  //D_contact   = 3;// delta contact for surface
  Jedi        = 0.;// jedi score :-)
  a           = 0.059*Vmax;//just to normalise the volume according the trainning dataset Vmax
  b           = 24.29;//coefficient hydrophobicity derived from PLS
  //  c           =  0.8032;//coefficient connectivity derived from PLS
  constant    = -13.39;//constant derived from PLS
  Vg          = resolution*resolution*resolution; // grid point volume nm^3
  time1       = getTimeStep();// time step from gromacs
  time2       = getStep();// step i from gromacs
  current     = time1*time2;// current snapshot (in ps)

  //size_protein     =  Apolar.size() + Polar.size();//total number of protein atoms in the CV
  size_apolar      =  Apolar.size();//total number of apolar atoms in the CV
  size_polar       =  Polar.size();//total number of polar atoms in the CV
  //size_connect_map   = 27;//maximum number of grid point surrounding a grid point ( used to speed up derivatives calculations)
  
  //connectivity_tot   = 0.;//sum of the connectivity scores of each grid point
  //connectivity       = 0.;//score connectivity
  hydrophobicity     = 0.;//score hydrophobicity
  hydrophobicity_tot = 0.;//sum of the hydrophobicity scores of each grid point
  D_hydro            = 0.05;

  double beta;
  //double connectivity_map[size_grid][27];
  //int short_list_grid[size_grid][27];

  beta = 5.0;

  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP1  :  Update grid coordinates
  //
  //-----------------------------------------
  //-----------------------------------------

  //int start;
  int max;
  //start = b_grid;
  max   = n_grid;
  //int nbr[max];
  
  // FIXME: This seem to do nothing useful?
  for (i=0; i<max; i++)
   {
     j = b_grid+i-1;
     //nbr[i] = j;
   }
  //cout << "nbr[i] is " << nbr[i] << " j is " << j << endl;

  // Check the translation of the center of mass of the binding site region

  // calculate COM from current snapshot
  double COM_x=0, sum_x=0;
  double COM_y=0, sum_y=0;
  double COM_z=0, sum_z=0;
  double mass_tot=0;
  // cout << "Checking translation of the center of mass" << endl; 
  
  //cout << " init sum_x  " << sum_x << " " << sum_y << " " << sum_z << " " << mass_tot << endl;
  for( j = 0; j < Apolar.size() + Polar.size() ; j++) 
    {
      sum_x += ( getMass(j) ) * ( getPosition(j)[0] );
      sum_y += ( getMass(j) ) * ( getPosition(j)[1] );
      sum_z += ( getMass(j) ) * ( getPosition(j)[2] );
      mass_tot += getMass(j);
      //cout << " j " << j;
      //cout << " jx " << getPosition(j)[0] << " jy " << getPosition(j)[1] << " jz " << getPosition(j)[2] << endl;
      //cout << "sum_x " << sum_x << " sum_y " << sum_y << " sum_z " << sum_z << " mass_tot " << mass_tot << endl;
    }

  COM_x = sum_x/mass_tot;
  COM_y = sum_y/mass_tot;
  COM_z = sum_z/mass_tot;
  //cout << "calculating COM over xyz " << COM_x << " " << COM_y << " " << COM_z << endl;    
  
  // JM OPT) replace system calls which subroutine to compute rotation matrix
  // generate a pdb file of the last snapshot of the trajectory (prev_snap)
  // to compute the rotation matrix.
  
  // Code below commented out while working out how to do this without system call
  // TODO) MAKE CODE BELOW WORK SO CAN COMPARE RESULTS OF ROTATION MATRIX OPERATION USING BOTH IMPLEMENTATION TO MAKE SURE 
  // RESULTS ARE CONSISTENT
  /*
  if (current > 0.0)
    {
      prev_snap = current-dump_matrix;
      //system(("echo 1 | trjconv -f md.xtc -o last.gro -b " + std::to_string(prev_snap) + " -e " + std::to_string(prev_snap) + " -s md.tpr > /dev/null 2>&1").c_str());
      string snap_string = static_cast<ostringstream*>( &(ostringstream() << prev_snap) )->str();
      string sys_string;
      sys_string = "echo 1 | trjconv -f md.xtc -o last.gro -b " + snap_string + " -e " + snap_string + " -s md.tpr > /dev/null 2>&1";;
      system(sys_string.c_str());
      //string str;
      //str = "echo 1 | trjconv -f md.xtc -o last.gro -b " + string("%s", prev_snap) + " -e " + string("%s", prev_snap) + " -s md.tpr > /dev/null 2>&1";
      cout << sys_string << endl;
      //cout << prev_snap << "@@" << string("%s", prev_snap);
      //THERE IS SOMETHING NOT RIGHT WITH THE SYNTAX OF THIS COMMAND HERE
    }
  // compute the rotation matrix using gromacs between the first and the last snapshots of the trajectory
  // using the heavy atom of the entire protein (group 3 in the default index file of gromacs)
  system("echo 2 | g_rotmat -f last.gro -s md.tpr > /dev/null 2>&1");
  remove("rotmat.xvg");
  remove("coord.xvg");
  remove("#last.gro.1#");

  // Read rotmat-jediv2 file ///////////
  ifstream file4("rotmat-jedi");// This output file is produced by the custom code change of gmx_rotmat.c in 4.5.5.
  while (getline(file4 , line))
    {
      istringstream ss(line);
      float value;
      int score_pos=0;
      while (ss >> value)
	{
	  score[score_pos] = value;
	  score_pos = score_pos +1;
	}
    }
  //   cout << "Reading rotmat-jedi file" << endl;
  
  //ifstream file4("rotmat-jediv2");// This output file is produced by the custom code change of gmx_rotmat.c in 4.5.5.
  //if(!file4)
  //  {
  //float 
  */
  // double ref_xlist[][3] : two dimensional array of coordinates for reference conformation   --> check PLUMED RMSD code to see how to store this data
  // double mov_xlist[][3] : two dimensional array of coordinates for current conformation     --> check PLUMED RMSD code to see how to access this data
  // int n_list            : the number of atoms used for the alignment                        --> easy
  // double mov_com[3]     : the centre of mass of the move_list (check if cog or com)         --> easy
  // double mov_to_ref[3]  : the vector between the com of mov and ref                         --> easy
  // double rotmat[3][3]   : the rotation matrix for least-squares fit                         --> the desired output
  // double* rmsd          : will contain the RMSD of the fit (but not needed for my purposes) --> calc for debugging purposes. remove if bottleneck.
  
  // Initialise data to pass
  int n_list = alignmentatoms.size();  
  double ref_xlist[n_list][3];//get this one directly from object
  double mov_xlist[n_list][3];
  double mov_com[3] = {0.0,0.0,0.0};
  double mov_to_ref[3];
  double rotmat[3][3];
  double rmsd = 0.0;

  double mov_mass_tot=0.0;
  
  for (unsigned i=0; i < n_list ; ++i)
    {
      ref_xlist[i][0] = ref_pos[i][0];//FIXME change calculate_rotation_rmsd args to take directly vector in
      ref_xlist[i][1] = ref_pos[i][1];
      ref_xlist[i][2] = ref_pos[i][2];
      Vector i_pos = getPosition( Apolar.size() + Polar.size() + i );// Not sure this is giving expected behavior
      mov_xlist[i][0] = i_pos[0];
      mov_xlist[i][1] = i_pos[1];
      mov_xlist[i][2] = i_pos[2];
      //cout << " i " << i << " mov_xlist " << mov_xlist[i][0] << " " << mov_xlist[i][1] << " " << mov_xlist[i][2] << endl;
      //cout << " i " << i << " ref_xlist " << ref_xlist[i][0] << " " << ref_xlist[i][1] << " " << ref_xlist[i][2] << endl;      
      // Also get data to compute mov_com
      double i_mass = getMass( Apolar.size() + Polar.size() + i );//FIXME this is a constant so could be cached
      mov_mass_tot += i_mass;
      mov_com[0] += i_mass * i_pos[0];
      mov_com[1] += i_mass * i_pos[1];
      mov_com[2] += i_mass * i_pos[2];
    }
  // Set mov_com and mov_to_ref
  mov_com[0] /= mov_mass_tot;
  mov_com[1] /= mov_mass_tot;
  mov_com[2] /= mov_mass_tot;

  mov_to_ref[0] = ref_com[0] - mov_com[0];//ref_com defined during CV init
  mov_to_ref[1] = ref_com[1] - mov_com[1];
  mov_to_ref[2] = ref_com[2] - mov_com[2];
  
  cout << "mov_com " << mov_com[0] << " " << mov_com[1] << " " << mov_com[2] << endl;
  cout << "ref_com " << ref_com[0] << " " << ref_com[1] << " " << ref_com[2] << endl;

  calculate_rotation_rmsd( ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, rotmat, &rmsd  );

  // NOTE: THE COM IS CALCULATED USING THE COM OF THE BINDING SITE REGION
  // HOWEVER THE ROTATION IS DONE A DIFFERENT SET OF ATOMS. USING A LARGER ALIGNMENT REGION (E.G PROTEIN BACKBONE ATOMS) HAS BEEN 
  // FOUND TO GIVE MORE STABLE GRID ROTATIONS

  cout << "Just called calculated_rotation_rmsd with dummy arguments " << endl;
  cout << "rotmat elements :" << endl;
  cout << rotmat[0][0] << " " << rotmat[0][1] << " " << rotmat[0][2] << endl;
  cout << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << endl;
  cout << rotmat[2][0] << " " << rotmat[2][1] << " " << rotmat[2][2] << endl; 

  //float score[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  //  }  
  // score contains the matrix elements of the 3x3 rotation matrix 
  //cout << "rotmat the old way " << endl;
  //cout << score[0] << " " << score[1] << " " << score[2] << endl;
  //cout << score[3] << " " << score[4] << " " << score[5] << endl;
  //cout << score[6] << " " << score[7] << " " << score[8] << endl;  
  //exit(0);

  for(i=0;i<max;i++) 
    {
      gridpoint = i;
      //new_x[i]  = ( score[0]*(grid_pos[gridpoint][0]) + score[3]*(grid_pos[gridpoint][1]) + score[6]*(grid_pos[gridpoint][2]) ) + (COM_x-COM_X) ;
      //new_y[i]  = ( score[1]*(grid_pos[gridpoint][0]) + score[4]*(grid_pos[gridpoint][1]) + score[7]*(grid_pos[gridpoint][2]) ) + (COM_y-COM_Y) ;
      //new_z[i]  = ( score[2]*(grid_pos[gridpoint][0]) + score[5]*(grid_pos[gridpoint][1]) + score[8]*(grid_pos[gridpoint][2]) ) + (COM_z-COM_Z) ;

      new_x[i]  = ( rotmat[0][0]*(grid_pos[gridpoint][0]) + rotmat[1][0]*(grid_pos[gridpoint][1]) + rotmat[2][0]*(grid_pos[gridpoint][2]) ) + (COM_x-COM_X) ;
      new_y[i]  = ( rotmat[0][1]*(grid_pos[gridpoint][0]) + rotmat[1][1]*(grid_pos[gridpoint][1]) + rotmat[2][1]*(grid_pos[gridpoint][2]) ) + (COM_y-COM_Y) ;
      new_z[i]  = ( rotmat[0][2]*(grid_pos[gridpoint][0]) + rotmat[1][2]*(grid_pos[gridpoint][1]) + rotmat[2][2]*(grid_pos[gridpoint][2]) ) + (COM_z-COM_Z) ;

      //cout << "grid_pos[gridpoint][0] is " << grid_pos[gridpoint][0] << endl;
      //cout << "grid_pos[gridpoint][1] is " << grid_pos[gridpoint][1] << endl;
      //cout << "grid_pos[gridpoint][2] is " << grid_pos[gridpoint][2] << endl;
      //cout << "new_x , new_y , new_z" << new_x[i] << "," << new_y[i] << "," << new_z[i] << endl; 
    }
  
  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP2  :  JEDI score  
  //
  //-----------------------------------------
  //-----------------------------------------
  //cout << " Starting Step 2 where JEDI score is calculated " << endl;
  double ray_score = .0;
  double rij[3];

  //----------> Compute activity of grid points and also VOLUME
  for(i = 0; i < max ; i++) 
    {
      ncoord = 0.;
      sum = 0.;
      //penalty_close = 0.;
      //penalty_far = 0.;
      min_dist = 0.;
      grd_x = new_x[i];
      grd_y = new_y[i];
      grd_z = new_z[i];
      k = 0;

      //-----> FIXME) Restructure to reuse code for loops over apolar/polar atom lists
      for( j = 0; j < Apolar.size(); j++) 
	{
	  rij[0] = grd_x - getPosition(j)[0];
	  rij[1] = grd_y - getPosition(j)[1];
	  rij[2] = grd_z - getPosition(j)[2];
	  // FIXME No PBC check !!  Code may need changes for PBC 
	  mod_rij   = sqrt( rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < 0.03)//FIXME: Constants
	    {
	      mod_rij = 0.03;
	    }
	  //cout << "grd_x " << grd_x << " y " << grd_y << " z " << grd_z << " jx " << getPosition(j)[0] << " jy " << getPosition(j)[1] << " jz " << getPosition(j)[2] << endl;
//      cout << " i " << i << " j " << j << " mod_rij " << mod_rij << endl;
	  // calculate mindist between grid points and protein atoms (EQUATION X.XX)
	  // also get number of neighbors in the same pass ...(EQUATION X.XX)
	  ncoord += s_off( 1.0, mod_rij, cutoff_far, D_far);
	  pe = exp(beta/mod_rij);
	  sum += pe;
	}
      //-----> polar check
      for( j = 0; j < Polar.size(); j++) 
	{
	  rij[0] = grd_x - getPosition( Apolar.size() + j )[0];
	  rij[1] = grd_y - getPosition( Apolar.size() + j)[1];
	  rij[2] = grd_z - getPosition( Apolar.size() + j)[2];
	  mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < 0.03)//FIXME Constant
	    {
	      mod_rij = 0.03;
	    }
	  //      cout << "grd_x " << grd_x << " y " << grd_y << " z " << grd_z << " jx " << getPosition(Apolar.size()+j)[0] << " jy " << getPosition(Apolar.size()+j)[1] << " jz " << getPosition(Apolar.size()+j)[2] << endl;
	  //      cout << " i " << i << " j " << j << " mod_rij " << mod_rij << endl;
	  // calculate mindist between grid points and protein atoms //
	  ncoord += s_off( 1.0, mod_rij, cutoff_far, D_far);
	  pe = exp(beta/mod_rij);
	  sum += pe;
	}
      sum_list[i] = sum;
      min_dist = beta/std::log(sum);
      min_dist_list[i] = min_dist;
      //NP[i] = ncoord;
      mind_list[i]=s_on( 1.0, min_dist, cutoff_close, D_close);
      //cout << "mind_dist " << min_dist << " " << cutoff_close << " " << D_close << endl;
      ncoord = 0.;
    }

  //FIXME: Comment what is going on here
  for( i = 0; i < max; i++) 
    {
      sum = .0;
      ray_score = 0.0;
      grd_x = new_x[i];
      grd_y = new_y[i];
      grd_z = new_z[i];
      for(j = 0; j < 44; j++)//FIXME number of neighbors is hardcoded 
	{
	  k = rays[i][j] - b_grid;
	  if (rays[i][j] > 0)
	    {
	      rij[0] = grd_x - new_x[k];
	      rij[1] = grd_y - new_y[k];
	      rij[2] = grd_z - new_z[k];
	      // FIXME: No PBC check !!  Code may need changes for PBC 
	      mod_rij   = sqrt( rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	      if (mod_rij < 0.03)// FIXME: Distance hardcoded
		{
		  mod_rij = 0.03;
		}
	      ray_score += s_off(1.0,min_dist_list[k],0.15,cutoff_far) * s_on(1.0,mod_rij,0.25,0.05) * s_off(1.0,mod_rij,0.3,0.05);// FIXME: Constants
	    }
	  else
	    {
	      break;
	    }
	}
      ray[i] = ray_score;
      penalty_far = ( s_on( 1.0, ray_score, 10.0, cutoff_enclosure));//FIXME Constants, documentation
      enclosure_list[i]  = ( s_on( 1.0, ray_score, 10.0, cutoff_enclosure));//FIX ME Constants, documentation
      
      penalty_list[i] = mind_list[i] * enclosure_list[i] * pocket[i];//FIXME variable naming

      if (0.0 < penalty_list[i])
	{
	  active_grid.push_back(i);
	}
      volume += penalty_list[i];
      //cout << "i " << i << " penalty_list[i] " << penalty_list[i] << " mind_list[i] " << mind_list[i] << " enclosure_list[i] "<< enclosure_list[i] << " pocket[i] " << pocket[i]  << " volume " <<  volume << endl;
    }
  
  //----------> "Drug-like" volume
  s1 = ( s_off( 1.0, volume, Vmax, D_Vmax));
  s2 = (s_on( 1.0, volume, Vmin, D_Vmin));
  Va = s1 * s2;

  //----------> Hydrophobicity
  for( i = 0; i < active_grid.size() ; i++) 
    {
      apolar = 0.;
      polar  = 0.;
      // FIX_ME loops over apolar/polar lists
      for( j = 0; j < Apolar.size() ; j++) 
	{
	  rij[0] = new_x[active_grid[i]] - getPosition(j)[0];
	  rij[1] = new_y[active_grid[i]] - getPosition(j)[1];
	  rij[2] = new_z[active_grid[i]] - getPosition(j)[2];
	  mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  apolar += s_off( 1.0, mod_rij, cutoff_hydro, D_hydro);//FIXME equation number
	}

      apolarity[active_grid[i]] = apolar;
      for( j = 0; j < Polar.size(); j++) 
	{
	  rij[0] = new_x[active_grid[i]] - getPosition(Apolar.size()+j)[0];
	  rij[1] = new_y[active_grid[i]] - getPosition(Apolar.size()+j)[1];
	  rij[2] = new_z[active_grid[i]] - getPosition(Apolar.size()+j)[2];
	  mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  polar += s_off(1.0,mod_rij,cutoff_hydro,D_hydro);//FIXME equation number
	}
      polarity[active_grid[i]] = polar;
      
      if( polar + apolar > 0.)
	{
	  hydrophobicity_list[active_grid[i]]=apolar/(apolar+polar);   
	}
      else
	{
	  hydrophobicity_list[active_grid[i]]=0.;   
	}
      hydrophobicity_tot+=hydrophobicity_list[active_grid[i]]*penalty_list[active_grid[i]];
    }
  
  if(volume > 0.)
    {
      hydrophobicity = hydrophobicity_tot/volume;   
    }
  else
    {
      hydrophobicity = 0.0;
    }
   //cout << "Va a volume Vmax b hydrophobicity constant " << endl;
  //cout << Va << " " << a << " " << volume << " " << Vmax << " " << b << " " << hydrophobicity << " " << constant << endl;

  //----------> JEDI SCORE
  Jedi=Va*(a*volume/Vmax+b*hydrophobicity+constant);//JEDI score without connectivity
  setValue(Jedi);
  // cout.precision(9);
  // cout << "current, Jedi, Va, hydrophobicity, volume/Vmax, COM_x, COM_y, COM_z, score[0], score[1], score[2], score[3], score[4], score[5], score[6], score[7], score[8]" << endl;
  // cout << current << " " << Jedi << " " << Va << " " << hydrophobicity << " " << volume/Vmax
  // << " " << COM_x << " " << COM_y << " " << COM_z << " " << score[0] << " " << score[1] << " "
  // << score[2] << " " << score[3] << " " << score[4] << " " << score[5] << " " << score[6] << " " << score[7] << " " << score[8] <<endl;
  

  cout << "current is " << current << " dump_matrix is " << dump_matrix << endl;
  cout << "fmod is " << fmod(current,dump_matrix) << endl;
  if( fmod(current,dump_matrix) < 0.00001 )//FIXME dodgy mathematics
    {
      ofstream wfile;
      wfile.open("jedi_output.dat",std::ios_base::app); // This command allows you to append data to a file that already exists
      //wfile << current << " " << Jedi << " " << Va << " " << hydrophobicity << " " << volume/Vmax
      //	    << " " << COM_x << " " << COM_y << " " << COM_z << " " << score[0] << " " << score[1] << " "
      //    << score[2] << " " << score[3] << " " << score[4] << " " << score[5] << " " << score[6] << " " 
      //    << score[7] << " " << score[8] <<endl;
      wfile << current << " " << Jedi << " " << Va << " " << hydrophobicity << " " << volume/Vmax
      	    << " " << COM_x << " " << COM_y << " " << COM_z << " " << rotmat[0][0] << " " << rotmat[0][1] << " "
          << rotmat[0][2] << " " << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << " " << rotmat[2][0] << " " 
          << rotmat[2][1] << " " << rotmat[2][2] << endl;
      wfile.close();
    }
 
  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP 3  :  Derivatives  
  //
  //-----------------------------------------
  //-----------------------------------------

  //FIXME: Compute derivatives at the same time as potential for speed 

  //double dV_dx;
  double ds1_dx, ds1_dy, ds1_dz;//derivatives Vmin
  double ds2_dx, ds2_dy, ds2_dz;//derivatives Vmax
  double da_dx, da_dy, da_dz;//derivatives activity
  double ds2_dm, dAp_dm, dAp_dx, dP_dx;


  //FIXME: Initialise arrays with '0' 
  double dVa_dx_apolar[Apolar.size()];
  double dVa_dx_polar[Polar.size()];

  double dV_dx_apolar[Apolar.size()];
  double dV_dx_polar[Polar.size()];

  double dVa_dy_apolar[Apolar.size()];
  double dVa_dy_polar[Polar.size()];

  double dV_dy_apolar[Apolar.size()];
  double dV_dy_polar[Polar.size()];

  double dVa_dz_apolar[Apolar.size()];
  double dVa_dz_polar[Polar.size()];

  double dV_dz_apolar[Apolar.size()];
  double dV_dz_polar[Polar.size()];

  double dH_dx_apolar[Apolar.size()];
  double dH_dx_polar[Polar.size()];

  double dH_dy_apolar[Apolar.size()];
  double dH_dy_polar[Polar.size()];

  double dH_dz_apolar[Apolar.size()];
  double dH_dz_polar[Polar.size()];

  double dHa_dx, dHa_dy, dHa_dz, dP_dm, ds1_dm;

  //FIXME: Magic numbers
  double da_dx_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV
  double da_dy_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV
  double da_dz_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV

  double da_dx_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV
  double da_dy_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV
  double da_dz_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV

  double dSenclosure_dx, dSmind_dx, dSenclosure_dy, dSmind_dy, dSenclosure_dz, dSmind_dz;//dSmin: just for CC pemalty
  double dsmind_dx, dsmind_dy, dsmind_dz;//dsmin: just for solvent exposed pemalty
  double dmind_dx, dmind_dy, dmind_dz, dmind_dr, dr_dx, dr_dy, dr_dz,dSmind_dm,dSenclosure_dm,s1_ray,s2_ray;

  double dsum_enclosure_x, dsum_enclosure_y, dsum_enclosure_z;// ds_off_dx, ds_on_dx, ds_off_dy, ds_on_dy, ds_off_dz, ds_on_dz;

  //----------> da_dx, da_dy, da_dz ( activity )

  //  cout << " loop i over active_grid.size() is " << active_grid.size() << endl; 
  for( i = 0; i < active_grid.size() ; i++) 
    {
      grd_x = new_x[active_grid[i]];
      grd_y = new_y[active_grid[i]];
      grd_z = new_z[active_grid[i]];
      // da_dx_apolar
      //      cout << " loop j over apolar.size() is " << Apolar.size() << endl;
      for( j = 0; j < Apolar.size() ; j++) 
	{
	  rij[0]    = new_x[active_grid[i]] - getPosition(j)[0];
	  rij[1]    = new_y[active_grid[i]] - getPosition(j)[1];
	  rij[2]    = new_z[active_grid[i]] - getPosition(j)[2];
	  mod_rij   = sqrt( rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < 0.03)// FIXME: Constant
	    {
	    mod_rij = 0.03;
	    }
	  dsum_enclosure_x = 0.0;

	  dmind_dr = ((beta*beta)*(exp(beta/mod_rij))) / ((mod_rij*mod_rij) * (sum_list[active_grid[i]]) * pow((std::log(sum_list[active_grid[i]])),2));
	  dSmind_dm = ds_on_dm(1.0,min_dist_list[active_grid[i]],cutoff_close,D_close) * 1.0/D_close;

	  dr_dx  = - (new_x[active_grid[i]]-getPosition(j)[0]) / (mod_rij);
	  dmind_dx  = dmind_dr*(dr_dx);
	  dSmind_dx = dSmind_dm * dmind_dx;
	  dsum_enclosure_y = 0.0;

	  dr_dy  = - (new_y[active_grid[i]]-getPosition(j)[1]) / (mod_rij);
	  dmind_dy  = dmind_dr*(dr_dy);
	  dSmind_dy = dSmind_dm * dmind_dy;
	  dsum_enclosure_z = 0.0;

	  dr_dz  = - (new_z[active_grid[i]]-getPosition(j)[2]) / (mod_rij);
	  dmind_dz  = dmind_dr*(dr_dz);
	  dSmind_dz = dSmind_dm * dmind_dz;

	  for( l = 0; l < 44 ; l++)//FIXME Magic number
	    {
	      if (rays[active_grid[i]][l] > 0)
		{
		  double rjk[3];
		  double rik[3];
		  k = rays[active_grid[i]][l] - b_grid;
		  //cout << " k is " << k << endl;
		  //cout << "j is " << j << endl;
		  //cout << "active_grid[ i ] " << active_grid[i] << " i " << i << " l " << l << " b_grid " << b_grid << endl;
		  //cout << " getPosition(0)[0] " << getPosition(j)[0] << endl;
		  //cout << "new_x[0] " << new_x[0] << endl;    
		  rjk[0]    = new_x[k] - getPosition(j)[0];//protein j and grid point k
		  //cout << " rjk[0] " << rjk[0] << endl;
		  rjk[1]    = new_y[k] - getPosition(j)[1];//protein j and grid point k
		  rjk[2]    = new_z[k] - getPosition(j)[2];//protein j and grid point k
		  rik[0]    = new_x[k] - grd_x;
		  rik[1]    = new_y[k] - grd_y;
		  rik[2]    = new_z[k] - grd_z;
		  // FIXME No PBC check !!  Code may need changes for PBC 
		  mod_rjk   = sqrt(rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]);//protein j and grid point k
		  mod_rik   = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);//grid point i and grid point k
		  //FIXME Magic numbers
		  if (mod_rjk < 0.03)
		    {
		      mod_rjk = 0.03;
		    }
		  if (mod_rik < 0.03)
		    {
		      mod_rik = 0.03;
		    }
		  s1_ray    = s_on( 1.0, mod_rik, 0.25, 0.05);//FIXME Constants
		  s2_ray    = s_off(1.0, mod_rik, 0.3, 0.05);// FIXME Constants
		  //cout << "rjk[0] " << rjk[0] << " rik[0] " << rik[0] << endl;
		  dmind_dr  = ((beta*beta)*(exp(beta/mod_rjk))) / ((mod_rjk*mod_rjk) * (sum_list[k]) * pow((std::log(sum_list[k])),2));
		  dSmind_dm = ds_off_dm(1.0,min_dist_list[k],0.15,cutoff_far) * 1.0/cutoff_far ;
		  //...with respect to x
		  dr_dx     = - (new_x[k]-getPosition(j)[0]) / (mod_rjk);
		  //cout << "dr_dx " << dr_dx << endl;
		  dmind_dx  = dmind_dr*(dr_dx);
		  dsmind_dx = dSmind_dm * dmind_dx;
		  dsum_enclosure_x += (dsmind_dx * ( s1_ray * s2_ray ));
		  //...with respect to y
		  dr_dy     = - (new_y[k]-getPosition(j)[1]) / (mod_rjk);
		  dmind_dy  = dmind_dr*(dr_dy);
		  dsmind_dy = dSmind_dm * dmind_dy;
		  dsum_enclosure_y += (dsmind_dy * ( s1_ray * s2_ray ));
		  //...with respect to z
		  dr_dz     = - (new_z[k]-getPosition(j)[2]) / (mod_rjk);
		  dmind_dz  = dmind_dr*(dr_dz);
		  dsmind_dz = dSmind_dm * dmind_dz;
		  dsum_enclosure_z += (dsmind_dz * ( s1_ray * s2_ray ));
		}
	      else
		{
		  break;
		}
	    }
	  dSenclosure_dm     = ds_on_dm( 1.0, ray[active_grid[i]], 10.0, cutoff_enclosure);
	  dSenclosure_dx     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_x;
	  da_dx_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) * pocket[active_grid[i]];
	  
	  dSenclosure_dy     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_y;
	  da_dy_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) *pocket[active_grid[i]];

	  dSenclosure_dz     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_z;
	  da_dz_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) *pocket[active_grid[i]];
	  
	  dsum_enclosure_x   = 0.;
	  dsum_enclosure_y   = 0.;
	  dsum_enclosure_z   = 0.;
	}
      // da_dx_polar
      //cout << " loop j over Polar.size() is " << Polar.size() << endl;
      for( j = 0; j < Polar.size() ; j++) 
	{
	  //jj = j - Apolar.size();
	  jj = j;//FIXME: Clarify need for 'jj'
	  rij[0]    = new_x[active_grid[i]] - getPosition(Apolar.size() + j)[0];
	  rij[1]    = new_y[active_grid[i]] - getPosition(Apolar.size() + j)[1];
	  rij[2]    = new_z[active_grid[i]] - getPosition(Apolar.size() + j)[2];
	  mod_rij   = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < 0.03)// FIXME constant
	    {
	      mod_rij = 0.03;
	    }
	  dmind_dr = ((beta*beta)*(exp(beta/mod_rij))) / ((mod_rij*mod_rij) * (sum_list[active_grid[i]]) * pow((std::log(sum_list[active_grid[i]])),2));
	  dSmind_dm = ds_on_dm(1.0,min_dist_list[active_grid[i]],cutoff_close,D_close) * 1.0/D_close;
	  dsum_enclosure_x = 0.0;
	  dr_dx  = - (new_x[active_grid[i]]-getPosition(Apolar.size()+j)[0]) / (mod_rij);
	  dmind_dx  = dmind_dr*(dr_dx);
	  dSmind_dx = dSmind_dm * dmind_dx;
	  
	  dsum_enclosure_y = 0.0;
	  dr_dy  = - (new_y[active_grid[i]]-getPosition(Apolar.size()+j)[1]) / (mod_rij);
	  dmind_dy  = dmind_dr*(dr_dy);
	  dSmind_dy = dSmind_dm * dmind_dy;

	  dsum_enclosure_z = 0.0;
	  dr_dz     = - (new_z[active_grid[i]]-getPosition(Apolar.size()+j)[2]) / (mod_rij);
	  dmind_dz  = dmind_dr*(dr_dz);
	  dSmind_dz = dSmind_dm * dmind_dz;
	  // cout << " loop l over rays "<< endl;
	  for( l = 0; l < 44 ; l++)
	    {
	      k = rays[active_grid[i]][l] - b_grid;
	      if (rays[active_grid[i]][l] > 0)
		{
		  double rjk[3];
		  double rik[3];
		  rjk[0]    = new_x[k]-getPosition(Apolar.size()+j)[0];//protein j and grid point k
		  rjk[1]    = new_y[k]-getPosition(Apolar.size()+j)[1];//protein j and grid point k
		  rjk[2]    = new_z[k]-getPosition(Apolar.size()+j)[2];//protein j and grid point k
		  rik[0]    = new_x[k]-grd_x;
		  rik[1]    = new_y[k]-grd_y;
		  rik[2]    = new_z[k]-grd_z;
		  // FIXME No PBC check !!  Code may need changes for PBC 
		  mod_rjk   = sqrt(rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]);//protein j and grid point k
		  mod_rik   = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);//grid point i and grid point k
		  if (mod_rjk < 0.03)
		    {
		      mod_rjk = 0.03;
		    }
		  s1_ray    = s_on(1.0,mod_rik,0.25,0.05);//FIXME Constants
		  s2_ray    = s_off(1.0,mod_rik,0.3,0.05);//FIXME Constants
		  dmind_dr  = ((beta*beta)*(exp(beta/mod_rjk))) / ((mod_rjk*mod_rjk) * (sum_list[k]) * pow((std::log(sum_list[k])),2));
		  dSmind_dm = ds_off_dm(1.0, min_dist_list[k], 0.15, cutoff_far)* 1.0/cutoff_far;
		  //...with respect to x
		  dr_dx     = - (new_x[k]-getPosition(Apolar.size()+j)[0]) / (mod_rjk);
		  dmind_dx  = dmind_dr*(dr_dx);
		  dsmind_dx = dSmind_dm * dmind_dx;
		  dsum_enclosure_x += (dsmind_dx * ( s1_ray * s2_ray ));
		  //...with respect to y
		  dr_dy     = - (new_y[k]-getPosition(Apolar.size()+j)[1]) / (mod_rjk);
		  dmind_dy  = dmind_dr*(dr_dy);
		  dsmind_dy = dSmind_dm * dmind_dy;
		  dsum_enclosure_y += (dsmind_dy * ( s1_ray * s2_ray ));
		  //...with respect to z
		  dr_dz     = - (new_z[k]-getPosition(Apolar.size()+j)[2]) / (mod_rjk);
		  dmind_dz  = dmind_dr*(dr_dz);
		  dsmind_dz = dSmind_dm * dmind_dz;
		  dsum_enclosure_z += (dsmind_dz * ( s1_ray * s2_ray ));
		}
	      else
		{
		  break;
		}
	    }
	  dSenclosure_dm     = ds_on_dm(1.0,ray[active_grid[i]],10.0,cutoff_enclosure) * 1.0/cutoff_enclosure;
	  dSenclosure_dx     = dSenclosure_dm * dsum_enclosure_x;
	  da_dx_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) *pocket[active_grid[i]];
	  dSenclosure_dy     = dSenclosure_dm * dsum_enclosure_y;
	  da_dy_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) *pocket[active_grid[i]];
	  dSenclosure_dz     = dSenclosure_dm * dsum_enclosure_z;
	  da_dz_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) *pocket[active_grid[i]];
	  //cout << " i " << i << " jj " << jj << " da_dz_polar " << da_dz_polar[i][jj] << " enclosure_list[active_grid[i]] " << enclosure_list[active_grid[i]] << " dSmind_dz " << dSmind_dz << " mind_list[active_grid[i] " << mind_list[active_grid[i]] << " dSenclosure_dz " << dSenclosure_dz << " pocket[active_grid[i]] " << pocket[active_grid[i]] << endl;
	  dsum_enclosure_x   = 0.;
	  dsum_enclosure_y   = 0.;
	  dsum_enclosure_z   = 0.;
      }
    }

  double vx,vy,vz;
  // apolar atoms first
  //  cout << " summing partial derivatives " << endl;
  for( j = 0; j < Apolar.size(); j++) 
    {
      jj      = j ;
      dHa_dx  = 0.;
      dHa_dy  = 0.;
      dHa_dz  = 0.;
      da_dx   = 0.;
      da_dy   = 0.;
      da_dz   = 0.;
      for( i = 0; i < active_grid.size() ; i++) 
	{
	  gridpoint = active_grid[i];
	  da_dx += da_dx_apolar[i][jj];
	  da_dy += da_dy_apolar[i][jj];
	  da_dz += da_dz_apolar[i][jj];

	  // Does this ever happen??
	  if( (polarity[gridpoint] + apolarity[gridpoint]) > 0.)
	    {
	      rij[0]    = new_x[gridpoint] - getPosition(j)[0];
	      rij[1]    = new_y[gridpoint] - getPosition(j)[1];
	      rij[2]    = new_z[gridpoint] - getPosition(j)[2];
	      mod_rij   = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);

	      dAp_dm  = ds_off_dm(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro);//same value for x,y,z
	      dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[0])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_apolar[i][jj]);//!!! dynamic allocation
	      dHa_dx += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));
	      dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[1])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dy_apolar[i][jj]);//!!! dynamic allocation
	      dHa_dy += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));
	      dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[2])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dz_apolar[i][jj]);//!!! dynamic allocation
	      dHa_dz += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));
	    }
	}

      dH_dx_apolar[jj] = dHa_dx;
      dH_dy_apolar[jj] = dHa_dy;
      dH_dz_apolar[jj] = dHa_dz;

      ds1_dm =ds_off_dm(1.0,volume,Vmax,D_Vmax) * (1.0/D_Vmax) * (Vg) ;

      dV_dx_apolar[jj] = (Vg) * (da_dx);
      dV_dy_apolar[jj] = (Vg) * (da_dy);
      dV_dz_apolar[jj] = (Vg) * (da_dz);
      
      ds1_dx = (ds1_dm) * (da_dx);
      ds1_dy = (ds1_dm) * (da_dy);
      ds1_dz = (ds1_dm) * (da_dz);

      ds2_dm = ds_on_dm(1.0,volume,Vmin,D_Vmin) * (1.0/D_Vmin) * (Vg) ;

      ds2_dx = (ds2_dm) * (da_dx);
      ds2_dy = (ds2_dm) * (da_dy);
      ds2_dz = (ds2_dm) * (da_dz);

      dVa_dx_apolar[jj] = s1 * ds2_dx + s2 * ds1_dx;
      dVa_dy_apolar[jj] = s1 * ds2_dy + s2 * ds1_dy;
      dVa_dz_apolar[jj] = s1 * ds2_dz + s2 * ds1_dz;

      //-----  derivatives 
      vx= Jedi * ((1.0/Va)*(dVa_dx_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dx_apolar[jj]) + b*(dH_dx_apolar[jj]) ));
      vy= Jedi * ((1.0/Va)*(dVa_dy_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dy_apolar[jj]) + b*(dH_dy_apolar[jj]) ));
      vz= Jedi * ((1.0/Va)*(dVa_dz_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dz_apolar[jj]) + b*(dH_dz_apolar[jj]) ));
      //cout << "Jedi * ((1.0/Va)*(dVa_dx_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dx_apolar[jj]) + b*(dH_dx_apolar[jj])))"<< endl;
      //cout << Jedi << " " << Va << " " << dVa_dx_apolar[jj] << " " << a << " " << volume << " " << b << " " << hydrophobicity << " " << constant << " " << dV_dx_apolar[jj] << " " << dH_dx_apolar[jj] << endl;
      //cout << "j " << j << " vx, vy, vz are " << vx << "," << vy << "," << vz << endl;     
      setAtomsDerivatives(j,Vector(vx,vy,vz));
      //cout << "APOLAR j " << j << " vx, vy, vz are " << vx << "," << vy << "," << vz << endl; 
    }
  
  // now polar atoms
  for( j = 0; j < Polar.size() ; j++) 
    {
      //jj = j - Apolar.size();
      jj = j;
      dHa_dx = 0.;
      dHa_dy = 0.;
      dHa_dz = 0.;
      da_dx = 0.;
      da_dy = 0.;
      da_dz = 0.;
      for( i = 0; i < active_grid.size() ; i++) 
	{
	  gridpoint = active_grid[i];
	  da_dx += da_dx_polar[i][jj];
	  da_dy += da_dy_polar[i][jj];
	  da_dz += da_dz_polar[i][jj];

        if ( (polarity[gridpoint] + apolarity[gridpoint]) > 0.)
	  {
	    rij[0]    = new_x[gridpoint]-getPosition(Apolar.size()+j)[0];
	    rij[1]    = new_y[gridpoint]-getPosition(Apolar.size()+j)[1];
	    rij[2]    = new_z[gridpoint]-getPosition(Apolar.size()+j)[2];
	    mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
	    std::max(mod_rij, 0.05);
	    dP_dm   = ds_off_dm(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro);
	    dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[0])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
	    dHa_dx += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
	    dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[1])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
	    dHa_dy += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
	    dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[2])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
	    dHa_dz += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
	  }
	}

      dH_dx_polar[jj] = dHa_dx;
      dH_dy_polar[jj] = dHa_dy;
      dH_dz_polar[jj] = dHa_dz;

      jj =j ;
      ds1_dm =ds_off_dm(1.0,volume,Vmax,D_Vmax) * (1.0/D_Vmax) * (Vg);

      dV_dx_polar[jj] = (Vg) * (da_dx);
      dV_dy_polar[jj] = (Vg) * (da_dy);
      dV_dz_polar[jj] = (Vg) * (da_dz);

      ds1_dx = (ds1_dm) * (da_dx);
      ds1_dy = (ds1_dm) * (da_dy);
      ds1_dz = (ds1_dm) * (da_dz);

      ds2_dm = ds_on_dm(1.0,volume,Vmin,D_Vmin) * (1.0/D_Vmin) * (Vg);

      ds2_dx = (ds2_dm) * (da_dx);
      ds2_dy = (ds2_dm) * (da_dy);
      ds2_dz = (ds2_dm) * (da_dz);

      dVa_dx_polar[jj] = s1 * ds2_dx + s2 * ds1_dx;
      dVa_dy_polar[jj] = s1 * ds2_dy + s2 * ds1_dy;
      dVa_dz_polar[jj] = s1 * ds2_dz + s2 * ds1_dz;

      //-----  derivatives without
      vx = Jedi * ((1.0/Va)*(dVa_dx_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dx_polar[jj]) + b*(dH_dx_polar[jj]) ));
      vy = Jedi * ((1.0/Va)*(dVa_dy_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dy_polar[jj]) + b*(dH_dy_polar[jj]) ));
      vz = Jedi * ((1.0/Va)*(dVa_dz_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dz_polar[jj]) + b*(dH_dz_polar[jj]) ));
      
      setAtomsDerivatives(j+Apolar.size(),Vector(vx,vy,vz));
      //cout << " jj " << jj << " dVa_dz_polar " << dVa_dz_polar[jj] << " dV_dz_polar " << dV_dz_polar[jj] << " dH_dz_polar[jj] " << dH_dz_polar[jj] << " da_dz " << da_dz << " ds1_dm " << ds1_dm << " ds2_dm " << ds2_dm << endl; 
      //cout << "POLAR jcount " << j+Apolar.size() << " vx, vy, vz are " << vx << "," << vy << "," << vz << endl;     
    }

  //TODO: Flag to output infrequently detailed grid statistics
  // dx files with updated position and 'activitiy'
  //                                and 'hydrophobicity'

}//close jedi::calculate
  
}//close namespace colvar
}//cloase namespace PLMD
