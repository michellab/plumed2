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
#include <iomanip>

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

void center_grid( vector<Vector> &grid_positions, double grid_ref_cog[3]);

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
  double resolution;
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
  //int gridstride;//frequency of output (in timestep) of a grid file;
  //DEPRECATED
  //vector<vector<int> > rays;
  //vector<vector<int> > neighbors;//DEPRECATED NOT USED??
  //vector<vector<double> > grid_pos;//DEPRECATED replaced by grid_positions
  //vector<double> pocket;//DEPRECATED, replaced by grid_s_off_bsi
  //float score[9];//rotation matrix elements
  //vector<AtomNumber> Apolar;//for plumed.dat file DEPRECATED replaced by apolaratoms
  //vector<AtomNumber> Polar;// for plumed.dat file DEPRECATED replaced by polaratoms

  double site_com[3];//reference coordinates of the center of mass of the binding site region
  double grid_ref_cog[3];//reference coordinates of the center of geometry of the grid
  //double COM_X, COM_Y, COM_Z;// coordinates of the center of mass of the binding site region

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
  stride=atoi(stride_string.c_str());
  cout << "STRIDE FROM INPUT IS " << stride << endl;
  //string summary_file;
  parse("SUMMARY",summary_file);
  string gridstride_string;
  parse("GRIDSTRIDE",gridstride_string);
  //int gridstride=atoi(gridstride_string.c_str());
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

  addValueWithDerivatives(); setNotPeriodic();

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

  // Work out grid resolution.
  // FIXME Here we ASSUME that all grid points are evenly spread and that we can infer their 
  // separation by computing the distance between the first TWO grid points. 
  double d2 = pow(grid_positions[1][0] - grid_positions[0][0],2) + pow(grid_positions[1][1] - grid_positions[0][1],2) + pow(grid_positions[1][2] - grid_positions[0][2],2);
  params.resolution = sqrt(d2);

  // Set maximum activity of grid points
  grid_s_off_bsi = set_bs_values(grid_positions, site_positions, params.theta, params.BSmin, params.deltaBS);

  //const std::vector<AtomNumber> gridnumbers = grid_pdb.getAtomNumbers();
  //cout << " grid has ? elements " << gridnumbers.size() << endl;
  neighbors = init_grid_neighbors( grid_positions,
				   params.GP1_min, params.deltaGP1,
				   params.GP2_min, params.deltaGP2);

  //Now center grid on origin by removing COG
  center_grid( grid_positions, grid_ref_cog );
  cout << "Grid ref cog " << grid_ref_cog[0] << " " << grid_ref_cog[1] << " " << grid_ref_cog[2] << endl;

  //Setup output
  // JM TODO: Check behavior CV init upon job restart.
  ofstream wfile;
  wfile.open(summary_file.c_str());
  wfile << "#step \t JEDI \t Vdrug_like \t Va \t Ha \t COM_x \t COM-y \t COM_z \t RotMat[0][0].Rotmat[0][1]....Rotmat[2][2]" << endl;
  wfile.close();

  //TODO CHECK IF FOLDER GRIDSTATS_FOLDER EXISTS
  // IF YES DELETE
  // CREATE NEW EMPTY FOLDER

  cout << "*** Completed initialisation JEDI collective variable" << endl;
  //exit(0);
}

vector<double> set_bs_values( vector<Vector> grid_pos,
			      vector<Vector> site_pos,
			      double theta, double BSmin, double deltaBS)
{
  vector<double> grid_s_off_bsi;

  int n_grid = grid_pos.size();

  for (int i=0; i < n_grid; i++)
    {
      grid_s_off_bsi.push_back(1.0);
    }

  int n_site = site_pos.size();
  if (n_site > 0)
    {
      for (int i=0; i < n_grid; i++)
	{
	  double sum=0.0;
	  for (int j=0; j < n_site; j++)
	    {
	      // Precompute first term eq 5 Cuchillo et al. JCTC 2015
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

  for (unsigned i=0; i < grid_pos.size(); i++)
    {
      vector<int> list;
      neighbors.push_back(list);
    }

  double dmin = pow(GP1_min,2);
  double dmax = pow(GP2_min+deltaGP2,2);

  for (unsigned i=0; i < grid_pos.size(); i++)
    {
      for (unsigned j=i+1; j < grid_pos.size(); j++)
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
	      //FIXME: Also precompute contributions to equation 7
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

  void center_grid( vector<Vector> &grid_positions, double grid_ref_cog[3] )
{
  int n_grid = grid_positions.size();
  double cog_x = 0.0;
  double cog_y = 0.0;
  double cog_z = 0.0;
  for (int i=0 ; i < n_grid ; i++)
    {
      cog_x += grid_positions[i][0];
      cog_y += grid_positions[i][1];
      cog_z += grid_positions[i][2];
    }
  cog_x /= n_grid;
  cog_y /= n_grid;
  cog_z /= n_grid;
  for (int i=0 ; i < n_grid ; i++)
    {
      grid_positions[i][0] = grid_positions[i][0] - cog_x;
      grid_positions[i][1] = grid_positions[i][1] - cog_y;
      grid_positions[i][2] = grid_positions[i][2] - cog_z;
    }
  grid_ref_cog[0] = cog_x;
  grid_ref_cog[1] = cog_y;
  grid_ref_cog[2] = cog_z;
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
      //s = 4*k*m*(pow((1.-pow(m,2)),2)) - 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
      double m2=m*m;
      s = 4*k*m*( (1.-m2)*(1-m2) - (1.-m2)*(1.+2*m2) );
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
      //s = (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2))));
      double m2=m*m;
      s = ( (1-m2)*(1-m2) ) * (1+2*m2);
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
      //s = -4*k*m*(pow((1.-pow(m,2)),2)) + 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
      double m2=m*m;
      s = -4*k*m*( (1.-m2)*(1-m2) - (1.-m2)*(1.+2*m2) );
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
      //s = 1. - ( (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2)))) );
      double m2=m*m;
      s = 1 -(1-m2)*(1-m2)*(1+2*m2);
    }
  return s;
 }

// calculator
void jedi::calculate(){

  //int i, j, jj, k, l, gridpoint;
  //int i,j,k,l,gridpoint;
  //double Jedi, volume, a, b, hydrophobicity, constant, grd_x, grd_y, grd_z;
  //double Vg, mod_rij, mod_rjk, mod_rik, pe, sum, min_dist;// penalty_close
  //double penalty_far;
  //double ncoord;// D_far;
  //double D_close;// D_enclosure;
  //double s1, s2, current, time1, time2, Va, Vmax, Vmin, D_Vmax, D_Vmin;//, prev_snap;
  //double apolar, polar, D_hydro, hydrophobicity_tot;
  //unsigned int size_apolar;
  //unsigned int size_polar;
  //vector<double> s3_m;// array with the derivative of distant contact
  //vector<double> s4_m; // array with the derivative of distant contact
  //vector<double> connectivity_list; // array with the number of active grid points around each grid points according to s_on
  //vector<double> hull_list;// array with the hull score of each grid points
  //vector<double> surface_list;// array with the surface score of each grid points
  //vector<double> ncontact_surface_list;// array number of protein atoms in contact with each grid point (for derivatives)
  //volume      = 0.;// score volume : total number of active and partially active grid points
  //D_close     = 0.05;// delta close contact
  //D_far       = 0.05;// delta distant contact
  //resolution  = 0.15;// grid spacing
  //s1          = 0.0;// minimum volume according our drug-like dataset
  //s2          = 0.0;// maximum volume according our drug-like dataset
  //Va          = 0.0;// penalty non drug-like pockets ( s1 * s2 )
  //Vmin        = cutoff_hull;// minimum volume if volume < Vmin --> s1 = 0
  //D_Vmin      = 10.0;// delta minimum volume volume > Vmin + D_Vmin --> s1 = 1
  //Vmax        = 90.0;// minimum volume if volume < Vmax --> s2 = 1
  //D_Vmax      = 50.0;// delta maximum volume volume > Vmax + D_Vmin --> s1 = 0
  //Jedi        = 0.;// jedi score :-)
  //a           = 0.059*Vmax;//just to normalise the volume according the trainning dataset Vmax
  //b           = 24.29;//coefficient hydrophobicity derived from PLS
  //  c           =  0.8032;//coefficient connectivity derived from PLS
  //constant    = -13.39;//constant derived from PLS
  //time1       = getTimeStep();// time step from gromacs
  //time2       = getStep();// step i from gromacs
  //current     = time1*time2;// current snapshot (in ps)

  //size_apolar      =  Apolar.size();//total number of apolar atoms in the CV
  //size_polar       =  Polar.size();//total number of polar atoms in the CV
  //hydrophobicity     = 0.;//score hydrophobicity
  //hydrophobicity_tot = 0.;//sum of the hydrophobicity scores of each grid point
  //D_hydro            = 0.05;

  int size_grid = grid_positions.size();
  vector<int> active_grid;// array with the index of active grid point (0 if inactive)
  //double activity[size_grid];//array of activity scores for each grid point
  vector<double> activity;
  activity.reserve(size_grid);
  //double s_on_exposure[size_grid];// array of s_on_exposure scores for each grid point
  vector<double> s_on_exposure;
  s_on_exposure.reserve(size_grid);
  //double s_on_mind[size_grid];// array with the penalty of close contact
  vector <double> s_on_mind;
  s_on_mind.reserve(size_grid);
  //double sum_dist[size_grid];// array with the sum of the distances between each grid point and all atoms in the CV
  vector<double> sum_dist;
  sum_dist.reserve(size_grid);
  //double min_dist_list[size_grid];// array with the minimum distance between each grid point and all atoms in the CV
  vector<double> min_dist_list;
  min_dist_list.reserve(size_grid);
  //double new_x[size_grid];// array with update of x coordinates of each grid points according to translation/rotation
  vector<double> new_x;
  new_x.reserve(size_grid);
  //double new_y[size_grid];// array with update of y coordinates of each grid points according to translation/rotation
  vector<double> new_y;
  new_y.reserve(size_grid);
  //double new_z[size_grid];// array with update of z coordinates of each grid points according to translation/rotation
  vector<double> new_z;
  new_z.reserve(size_grid);
  //double hydrophobicity_list[size_grid];
  vector<double> hydrophobicity_list;
  hydrophobicity_list.reserve(size_grid);

  //double apolarity[size_grid];// array number of apolar protein atoms around each grid point
  vector<double> apolarity;
  apolarity.reserve(size_grid);
  //double polarity[size_grid];// array number of polar protein atoms around each grid point
  vector<double> polarity;
  polarity.reserve(size_grid);
  //double exposure[size_grid];
  vector<double> exposure;
  exposure.reserve(size_grid);

  double Vg   = pow(params.resolution,3); // grid point volume nm^3
  double beta = 5.0;
  double min_modr = 0.03;

  //cout << "*** Beginning calculation of JEDI collective variable" << endl;

  if(pbc) makeWhole();


  //exit(0);

  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP1  :  Update grid coordinates
  //
  //-----------------------------------------
  //-----------------------------------------

  //int start;
  //int max;
  //max   = n_grid;
  // FIXME: This seem to do nothing useful?
  //for (i=0; i<max; i++)
  // {
  //   j = b_grid+i-1;
  //   //nbr[i] = j;
  // }
  //cout << "nbr[i] is " << nbr[i] << " j is " << j << endl;

  // Check the translation of the center of mass of the binding site region

  // calculate COM from current snapshot
  //double COM_x=0, sum_x=0;
  //double COM_y=0, sum_y=0;
  //double COM_z=0, sum_z=0;

  //double sum_x=0;
  //double sum_y=0;
  //double sum_z=0;
  //double mass_tot=0;
  // cout << "Checking translation of the center of mass" << endl; 
  //cout << " init sum_x  " << sum_x << " " << sum_y << " " << sum_z << " " << mass_tot << endl;

  //cout << "Apolar has size " << Apolar.size() << endl;
  //cout << " Polar has size " << Polar.size() << endl;

  //for( j = 0; j < Apolar.size() + Polar.size() ; j++)
  //int jstart = apolaratoms.size() + polaratoms.size();
  //for( unsigned j = jstart; j < jstart + alignmentatoms.size()  ; j++)
  //  {
  //    sum_x += ( getMass(j) ) * ( getPosition(j)[0] );
  //    sum_y += ( getMass(j) ) * ( getPosition(j)[1] );
  //    sum_z += ( getMass(j) ) * ( getPosition(j)[2] );
  //    mass_tot += getMass(j);
  //    //cout << " J " << j << " jx " << getPosition(j)[0] << " jy " << getPosition(j)[1] << " jz " << getPosition(j)[2] << endl;
  //    //cout << "sum_x " << sum_x << " sum_y " << sum_y << " sum_z " << sum_z << " mass_tot " << mass_tot << endl;
  //  }
  //
  //double COM_x = sum_x/mass_tot;
  //double COM_y = sum_y/mass_tot;
  //double COM_z = sum_z/mass_tot;
  //cout << "calculating COM over xyz " << COM_x << " " << COM_y << " " << COM_z << endl;    

  // double ref_xlist[][3] : two dimensional array of coordinates for reference conformation   --> check PLUMED RMSD code to see how to store this data
  // double mov_xlist[][3] : two dimensional array of coordinates for current conformation     --> check PLUMED RMSD code to see how to access this data
  // int n_list            : the number of atoms used for the alignment                        --> easy
  // double mov_com[3]     : the centre of mass of the move_list (check if cog or com)         --> easy
  // double mov_to_ref[3]  : the vector between the com of mov and ref                         --> easy
  // double rotmat[3][3]   : the rotation matrix for least-squares fit                         --> the desired output
  // double* rmsd          : will contain the RMSD of the fit (but not needed for my purposes) --> calc for debugging purposes. remove if bottleneck.

  // Initialise data to pass
  int n_align = alignmentatoms.size();
  // JM: LOOK UP BEST WAY TO RESERVE 2D ARRAY IN C++
  double ref_xlist[n_align][3];//get this one directly from object
  double mov_xlist[n_align][3];
  double mov_com[3] = {0.0,0.0,0.0};
  double mov_to_ref[3];
  double rotmat[3][3];
  double rmsd = 0.0;

  double mov_mass_tot=0.0;

  for (int i=0; i < n_align ; ++i)
    {
      ref_xlist[i][0] = ref_pos[i][0];//FIXME change calculate_rotation_rmsd args to take directly vector in
      ref_xlist[i][1] = ref_pos[i][1];
      ref_xlist[i][2] = ref_pos[i][2];
      //Vector i_pos = getPosition( Apolar.size() + Polar.size() + i );// Not sure this is giving expected behavior
      Vector i_pos = getPosition( apolaratoms.size() + polaratoms.size() + i );// MUST CHECK that PBC DO NOT MESS THINGS UP
      mov_xlist[i][0] = i_pos[0];
      mov_xlist[i][1] = i_pos[1];
      mov_xlist[i][2] = i_pos[2];
      //cout << " i " << i << " mov_xlist " << mov_xlist[i][0] << " " << mov_xlist[i][1] << " " << mov_xlist[i][2] << endl;
      //cout << " i " << i << " ref_xlist " << ref_xlist[i][0] << " " << ref_xlist[i][1] << " " << ref_xlist[i][2] << endl;      
      // Also get data to compute mov_com
      double i_mass = getMass( apolaratoms.size() + polaratoms.size() + i );//FIXME this is a constant so could be cached
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

  //  cout << "mov_com " << mov_com[0] << " " << mov_com[1] << " " << mov_com[2] << endl;
  //cout << "ref_com " << ref_com[0] << " " << ref_com[1] << " " << ref_com[2] << endl;
  //cout << " mov_to_ref " << mov_to_ref[0] << " " << mov_to_ref[1] << " " << mov_to_ref[2] << endl;


  // Here for debugging purposes, write oordinates of mov-xlist are indeed those of the alignement atoms
  //ofstream wfile000;
  //wfile000.open("mov-xlist-bef.xyz");
  //wfile000 << n_align << endl;
  //wfile000 << "comment" << endl;
  //for (int i=0; i < n_align; i++)
  //  {
  //    wfile000<< "C " << std::fixed << std::setprecision(5) << mov_xlist[i][0]*10 << " " << mov_xlist[i][1]*10 << " " << mov_xlist[i][2]*10 << endl;
  //  }
  //wfile000.close();

  calculate_rotation_rmsd( ref_xlist, mov_xlist, n_align, mov_com, mov_to_ref, rotmat, &rmsd  );

  // Here mov_xlist should end up centered on origin
  //ofstream wfile001;
  //wfile001.open("mov-xlist.xyz");
  //wfile001 << n_align << endl;
  //wfile001 << "comment" << endl;
  //for (int i=0; i < n_align; i++)
  //  {
  //    wfile001<< "C " << std::fixed << std::setprecision(5) << mov_xlist[i][0]*10 << " " << mov_xlist[i][1]*10 << " " << mov_xlist[i][2]*10 << endl;
  //  }
  //wfile001.close();

  //ofstream wfile002;
  //wfile002.open("ref-xlist.xyz");
  //wfile002 << n_align << endl;
  //wfile002 << "comment" << endl;
  //for (int i=0; i < n_align; i++)
  //  {
  //    wfile002<< "C " << std::fixed << std::setprecision(5) << ref_xlist[i][0]*10 << " " << ref_xlist[i][1]*10 << " " << ref_xlist[i][2]*10 << endl;
  //  }
  //wfile002.close();

  // Normally expect low rmsd value
  //cout << " The rmsd is " << rmsd << endl;

  // OPTIMISATION?
  // For systems that fluctuate a lot in alignement coordinates it may be better 
  // to update the reference coordinates after one iteration, so as to keep
  // rmsd fits as low as possible

  //rotmat[0][0] = 1.0;
  //rotmat[0][1] = 0.0;
  //rotmat[0][2] = 0.0;
  //rotmat[1][0] = 0.0;
  //rotmat[1][1] = 1.0;
  //rotmat[1][2] = 0.0;
  //rotmat[2][0] = 0.0;
  //rotmat[2][1] = 0.0;
  //rotmat[2][2] = 1.0;

  //cout << "Just called calculated_rotation_rmsd" << endl;
  //cout << "rotmat elements :" << endl;
  //cout << rotmat[0][0] << " " << rotmat[0][1] << " " << rotmat[0][2] << endl;
  //cout << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << endl;
  //cout << rotmat[2][0] << " " << rotmat[2][1] << " " << rotmat[2][2] << endl; 

  double new_grid_cog_x, local_grid_cog_x;
  double new_grid_cog_y, local_grid_cog_y;
  double new_grid_cog_z, local_grid_cog_z;
  // Apply same translation to grid_cog
  new_grid_cog_x = grid_ref_cog[0] - mov_to_ref[0];
  new_grid_cog_y = grid_ref_cog[1] - mov_to_ref[1];
  new_grid_cog_z = grid_ref_cog[2] - mov_to_ref[2];
  //cout << " New grid cog " << new_grid_cog_x << " " << new_grid_cog_y << " " << new_grid_cog_z << endl;
  // Now shift origin to mov_com
  local_grid_cog_x = new_grid_cog_x - mov_com[0];
  local_grid_cog_y = new_grid_cog_y - mov_com[1];
  local_grid_cog_z = new_grid_cog_z - mov_com[2];
  //cout << " local grid cog " << local_grid_cog_x << " " << local_grid_cog_y << " " << local_grid_cog_z << endl;
  // Now rotate about origin
  local_grid_cog_x = rotmat[0][0]*local_grid_cog_x + rotmat[1][0]*local_grid_cog_y + rotmat[2][0]*local_grid_cog_z;
  local_grid_cog_y = rotmat[0][1]*local_grid_cog_x + rotmat[1][1]*local_grid_cog_y + rotmat[2][1]*local_grid_cog_z;
  local_grid_cog_z = rotmat[0][2]*local_grid_cog_x + rotmat[1][2]*local_grid_cog_y + rotmat[2][2]*local_grid_cog_z;  
  //cout << " local grid cog " << local_grid_cog_x << " " << local_grid_cog_y << " " << local_grid_cog_z << endl;
  // Now map back to Cartesian
  new_grid_cog_x = local_grid_cog_x + mov_com[0];
  new_grid_cog_y = local_grid_cog_y + mov_com[1];
  new_grid_cog_z = local_grid_cog_z + mov_com[2];

  //new_grid_cog_x = grid_ref_cog[0];
  //new_grid_cog_y = grid_ref_cog[1];
  //new_grid_cog_z = grid_ref_cog[2];

  //new_grid_cog_x = -5.68689;
  //new_grid_cog_y = 1.53507;
  //new_grid_cog_z = -0.213876;

  //rotmat[0][0] = 1.0;
  //rotmat[0][1] = 0.0;
  //rotmat[0][2] = 0.0;
  //rotmat[1][0] = 0.0;
  //rotmat[1][1] = 1.0;
  //rotmat[1][2] = 0.0;
  //rotmat[2][0] = 0.0;
  //rotmat[2][1] = 0.0;
  //rotmat[2][2] = 1.0;

  //cout << " New grid cog " << new_grid_cog_x << " " << new_grid_cog_y << " " << new_grid_cog_z << endl;

  // Now rotate all grid points at origin and then translate to new cog
  for(int i=0;i< size_grid;i++)
    {
      new_x[i]  = ( rotmat[0][0]*(grid_positions[i][0]) + rotmat[1][0]*(grid_positions[i][1]) + rotmat[2][0]*(grid_positions[i][2]) ) + new_grid_cog_x;
      new_y[i]  = ( rotmat[0][1]*(grid_positions[i][0]) + rotmat[1][1]*(grid_positions[i][1]) + rotmat[2][1]*(grid_positions[i][2]) ) + new_grid_cog_y;
      new_z[i]  = ( rotmat[0][2]*(grid_positions[i][0]) + rotmat[1][2]*(grid_positions[i][1]) + rotmat[2][2]*(grid_positions[i][2]) ) + new_grid_cog_z;

      //new_x[i]  = ( rotmat[0][0]*(grid_positions[i][0]) + rotmat[0][1]*(grid_positions[i][1]) + rotmat[0][2]*(grid_positions[i][2]) ) + new_grid_cog_x;
      //new_y[i]  = ( rotmat[1][0]*(grid_positions[i][0]) + rotmat[1][1]*(grid_positions[i][1]) + rotmat[1][2]*(grid_positions[i][2]) ) + new_grid_cog_y;
      //new_z[i]  = ( rotmat[2][0]*(grid_positions[i][0]) + rotmat[2][1]*(grid_positions[i][1]) + rotmat[2][2]*(grid_positions[i][2]) ) + new_grid_cog_z;


      //cout << "grid_positions[gridpoint][0] is " << grid_positions[gridpoint][0] << endl;
      //cout << "grid_positions[gridpoint][1] is " << grid_positions[gridpoint][1] << endl;
      //cout << "grid_positions[gridpoint][2] is " << grid_positions[gridpoint][2] << endl;
      //cout << "new_x , new_y , new_z" << new_x[i] << "," << new_y[i] << "," << new_z[i] << endl; 
    }

  // Here for debugging purposes, write oordinates of updated grid to a XYZ file
  // This can be used to check that the grid has been trans/roted correctly
  //ofstream wfile;
  //wfile.open("grid-step1.xyz");
  //wfile << size_grid << endl;
  //wfile << "comment" << endl;
  //for (int i=0; i < size_grid; i++)
  //  {
  //    wfile << "C " << std::fixed << std::setprecision(5) << new_x[i]*10 << " " << new_y[i]*10 << " " << new_z[i]*10 << endl;
  //  }
  //wfile.close();

  //cout << "*** Getting ready for STEP 2" << endl;

  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP2  :  JEDI score
  //
  //-----------------------------------------
  //-----------------------------------------
  //cout << " Starting Step 2 where JEDI score is calculated " << endl;
  double rij[3];

  //----------> Compute activity of grid points and also VOLUME
  int n_apolarpolar = apolaratoms.size() + polaratoms.size();
  for(int i = 0; i < size_grid ; i++)
    {
      //ncoord = 0.;
      double sum = 0.;
      //penalty_close = 0.;
      //penalty_far = 0.;
      double grd_x = new_x[i];
      double grd_y = new_y[i];
      double grd_z = new_z[i];
      //int k = 0;

      //-----> FIXME) Restructure to reuse code for loops over apolar/polar atom lists
      for( int j = 0; j < n_apolarpolar; j++)
	{
	  rij[0] = grd_x - getPosition(j)[0];
	  rij[1] = grd_y - getPosition(j)[1];
	  rij[2] = grd_z - getPosition(j)[2];
	  // FIXME No PBC check !!  Code may need changes for PBC
	  // FIXME avoid sqrt if possible
	  double mod_rij   = sqrt( rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < min_modr)
	      mod_rij = min_modr;
	  //cout << "grd_x " << grd_x << " y " << grd_y << " z " << grd_z << " jx " << getPosition(j)[0] << " jy " << getPosition(j)[1] << " jz " << getPosition(j)[2] << endl;
//      cout << " i " << i << " j " << j << " mod_rij " << mod_rij << endl;
	  // calculate mindist between grid points and protein atoms (EQUATION 5 2nd term)
	  //I don't get what we do with ncoord also get number of neighbors in the same pass ...(EQUATION 7 1st term)
	  //ncoord += s_off( 1.0, mod_rij, cutoff_far, D_far);
	  //ncoord += s_off( 1.0, mod_rij, params.CC2_min, params.deltaCC2);
	  double pe = exp(beta/mod_rij);
	  sum += pe;
	}
      sum_dist[i] = sum;
      double min_dist = beta/std::log(sum);
      min_dist_list[i] = min_dist;
      //NP[i] = ncoord;
      //s_on_mind[i]=s_on( 1.0, min_dist, cutoff_close, D_close);
      s_on_mind[i] = s_on( 1.0, min_dist, params.CC_mind, params.deltaCC);
      //cout << " i " << i << " mind_dist " << min_dist << " " << params.CC_mind << " " << params.deltaCC << endl;
    }

  //FIXME: Comment what is going on here
  // For each grid point...
  double volume=0;

  for( int i = 0; i < size_grid; i++)
    {
      //double sum = .0;
      double exposure_score = 0.0;
      double grd_x = new_x[i];
      double grd_y = new_y[i];
      double grd_z = new_z[i];
      vector<int> neighbors_i = neighbors[i];
      for(unsigned j = 0; j < neighbors_i.size(); j++)
	{
	  int k = neighbors_i[j];
	  double rij[3];
	  rij[0] = grd_x - new_x[k];
	  rij[1] = grd_y - new_y[k];
	  rij[2] = grd_z - new_z[k];
	  // FIXME: No PBC check !!  Code may need changes for PBC
	  // OPTME. Distances should be constant so mod_rij result pre-computed !
	  double mod_rij   = sqrt( rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  if (mod_rij < min_modr)// FIXME: Distance hardcoded
	    mod_rij = min_modr;
	  // This is equation 7 of the paper so exposure_score is s_on_exposure_i
	  // OPTM term two and three should be constant and precomputed !!
	  //cout << " k " << k << " min_dist_list[k] " << min_dist_list[k] << endl;
	  double term1 = s_off(1.0,min_dist_list[k],params.CC2_min,params.deltaCC2);
	  double term2 = s_on(1.0,mod_rij,params.GP1_min,params.deltaGP1);
	  double term3 = s_off(1.0,mod_rij,params.GP2_min,params.deltaGP2);
	  //cout << " term1 " << term1 << " term2 " << term2 << " term3 " << term3 << endl;
	  exposure_score += term1 * term2 * term3;
	}
      exposure[i] = exposure_score;
      s_on_exposure[i]  = s_on( 1.0, exposure_score, params.Emin, params.deltaE);
      // This now gives the activity value from equation 5
      activity[i] = s_on_mind[i] * s_on_exposure[i] * grid_s_off_bsi[i];

      // Only grid points with activity > 0 are kept for the next calculations
      if (activity[i] > 0.0)
	active_grid.push_back(i);
      volume += activity[i];
      //cout << "i " << i << " activity[i] " << activity[i] << " s_on_mind[i] " << s_on_mind[i] << " s_on_exposure[i] "<< s_on_exposure[i] << " grid_s_off_bsi[i] " << grid_s_off_bsi[i]  << " exposure[i] " << exposure[i] << " volume " <<  volume << endl;
    }

  //cout << "There are " << active_grid.size() << " active grid points " << endl;
  double sum_activity = volume;
  double sum_activity2 = sum_activity*sum_activity;
  volume *= Vg;//Equation 4 of the paper
  // ----------> "Drug-like" volume
  double s_off_V = (s_off( 1.0, volume, params.V_max, params.deltaV_max));
  double s1 = s_off_V;
  double s2 = (s_on( 1.0, volume, params.V_min, params.deltaV_min));
  double s_on_V = s2;
  double Vdrug_like = s_off_V * s_on_V;
  //cout << " volume " << volume << " Vdrug_like " << Vdrug_like << endl;

  //----------> Hydrophobicity
  double hydrophobicity_tot=0.0;
  for( unsigned i = 0; i < active_grid.size() ; i++)
    {
      double apolar = 0.;
      for( unsigned j = 0; j < apolaratoms.size() ; j++)
	{
	  rij[0] = new_x[active_grid[i]] - getPosition(j)[0];
	  rij[1] = new_y[active_grid[i]] - getPosition(j)[1];
	  rij[2] = new_z[active_grid[i]] - getPosition(j)[2];
	  double mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  apolar += s_off( 1.0, mod_rij, params.r_hydro, params.deltar_hydro);//equation 12a paper (but multiplie dby ai later)
	}
      apolarity[active_grid[i]] = apolar;

      double polar  = 0.;
      for( unsigned j = 0; j < polaratoms.size(); j++)
	{
	  rij[0] = new_x[active_grid[i]] - getPosition(apolaratoms.size()+j)[0];
	  rij[1] = new_y[active_grid[i]] - getPosition(apolaratoms.size()+j)[1];
	  rij[2] = new_z[active_grid[i]] - getPosition(apolaratoms.size()+j)[2];
	  double mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
	  polar += s_off(1.0,mod_rij,params.r_hydro,params.deltar_hydro);//equation 12b, but multiplied by ai later
	}
      polarity[active_grid[i]] = polar;

      if( polar + apolar > 0.)
	  hydrophobicity_list[active_grid[i]]=apolar/(apolar+polar);
      else
	  hydrophobicity_list[active_grid[i]]=0.;


      //      cout << " i " << i << " active grid point " << active_grid[i] << " apolar " << apolar << " polar " << polar << endl;
      //cout << " hydrophobicity " << hydrophobicity_list[active_grid[i]] <<  endl;
      hydrophobicity_tot+=hydrophobicity_list[active_grid[i]]*activity[active_grid[i]];
    }

  //  cout << " hydrophobicity_tot " << hydrophobicity_tot << " volume " << volume << endl;

  // JM Feb 2016. Error in paper. Equation 10 numerator needs to be multiplied by Vg to get units
  double Ha;
  if(volume > 0.)
    Ha = hydrophobicity_tot*Vg/volume;//equation 10
  else
    Ha = 0.0;
    //cout << "Va a volume Vmax b hydrophobicity constant " << endl;
  //cout << Va << " " << a << " " << volume << " " << Vmax << " " << b << " " << hydrophobicity << " " << constant << endl;

  //----------> JEDI SCORE

  double Va = volume/params.V_max;
  double Jedi=Vdrug_like*(params.alpha*Va + params.beta*Ha + params.gamma);//JEDI score without connectivity
  setValue(Jedi);

  //cout << "Jedi score is " << Jedi << " alpha " << params.alpha << " Va " << Va << " beta " << params.beta << " Ha " << Ha << " gamma " << params.gamma << endl;
  //exit(0);

  // cout.precision(9);
  // cout << "current, Jedi, Vdrug_like, Ha, volume/Vmax, COM_x, COM_y, COM_z, score[0], score[1], score[2], score[3], score[4], score[5], score[6], score[7], score[8]" << endl;
  // cout << current << " " << Jedi << " " << Vdrug_like << " " << Ha << " " << volume/Vmax
  // << " " << COM_x << " " << COM_y << " " << COM_z << " " << score[0] << " " << score[1] << " "
  // << score[2] << " " << score[3] << " " << score[4] << " " << score[5] << " " << score[6] << " " << score[7] << " " << score[8] <<endl;

  int step= getStep();
  int mod =fmod(step,stride);
  //cout << " we are at step " << step << " we dump ? " << mod << endl;
  if (!mod)
    {
      ofstream wfile;
      //wfile.open("jedi_output.dat",std::ios_base::app); // This command allows you to append data to a file that already exists
      wfile.open(summary_file.c_str(),std::ios_base::app);
      //wfile << current << " " << Jedi << " " << Vdrug_like << " " << Ha << " " << Va
      //	    << " " << COM_x << " " << COM_y << " " << COM_z << " " << score[0] << " " << score[1] << " "
      //    << score[2] << " " << score[3] << " " << score[4] << " " << score[5] << " " << score[6] << " " 
      //    << score[7] << " " << score[8] <<endl;
      wfile << step << " " << Jedi << " " << Vdrug_like << " " << Va << " " << Ha 
      	    << " " << new_grid_cog_x << " " << new_grid_cog_y << " " << new_grid_cog_z << " " << rotmat[0][0] << " " << rotmat[0][1] << " "
          << rotmat[0][2] << " " << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << " " << rotmat[2][0] << " " 
          << rotmat[2][1] << " " << rotmat[2][2] << endl;
      wfile.close();
    }

  //FIXME: Add code to write grid coordinates as a dx file
  // Write with GRIDSTRIDE frequency
  // grid_stepnumber_activity.dx
  // grid_stepnumber_Ha.dx
  // These can then be visualized together with a gromacs trajectory to check
  // how the grid is moving over time, and to analyse how the jedi components
  // fluctuate over time

  //exit(0);

  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP 3  :  Derivatives
  //
  //-----------------------------------------
  //-----------------------------------------

  //FIXME: Compute derivatives at the same time as potential for speed
  //cout << " Doing derivatives " << endl;

  unsigned n_apolar = apolaratoms.size();
  unsigned n_polar = polaratoms.size();
  unsigned n_cv_atoms = n_apolar + n_polar;
  unsigned n_active_grid = active_grid.size();

  vector<double> d_ai_xpj_vec;
  d_ai_xpj_vec.reserve(n_active_grid);
  vector<double> d_ai_ypj_vec;
  d_ai_ypj_vec.reserve(n_active_grid);
  vector<double> d_ai_zpj_vec;
  d_ai_zpj_vec.reserve(n_active_grid);
  vector<double> dij_vec;
  dij_vec.reserve(n_active_grid);

  //n_cv_atoms=0;
  for ( unsigned j=0; j < n_cv_atoms ; j++)
    {
      double xj = getPosition(j)[0];
      double yj = getPosition(j)[1];
      double zj = getPosition(j)[2];
      // Compute activity derivatives (should we include pts with 0 activity?)
      double sum_d_ai_xpj=0.0;
      double sum_d_ai_ypj=0.0;
      double sum_d_ai_zpj=0.0;
      for (unsigned int i=0; i < n_active_grid ; i++)
	{
	  int grid_index = active_grid[i];
	  double xi = new_x[grid_index];
	  double yi = new_y[grid_index];
	  double zi = new_z[grid_index];
	  // d_rij over x/y/z
	  double dij_x = xi-xj;
	  double dij_y = yi-yj;
	  double dij_z = zi-zj;
	  //cout << " grid index " << grid_index << endl;
	  //cout << " xi " << xi << " yi " << yi << " zi " << zi << endl;
	  //cout << " xj " << xj << " yj " << yj << " zj " << zj << endl;
	  double dij2 = dij_x*dij_x + dij_y*dij_y + dij_z*dij_z;
	  double dij = sqrt(dij2);
	  if (dij < min_modr)//JM since also done in JEDI CV calc
	    {
	      dij = min_modr;
	      dij2 = dij*dij;
	    }
	  dij_vec[i] = dij;
	  //cout << " dij " << dij << endl;
	  //exit(0);
	  // COULD USE CUTOFFS TO SPEED UP CALCULATION?
	  double d_rij_xpj = - dij_x / dij;
	  double d_rij_ypj = - dij_y / dij;
	  double d_rij_zpj = - dij_z / dij;
	  // d_mindi over d_rij
	  double num = (beta*beta)*(exp(beta/dij));
	  // sum_dist was saved during calculation of activities in step 2
	  double den = dij2*sum_dist[grid_index]*pow(std::log(sum_dist[grid_index]),2);
	  double d_mindi_rij = num/den;
	  // d_mindi over x/y/z
	  double d_mindi_xpj = d_mindi_rij*d_rij_xpj;
	  double d_mindi_ypj = d_mindi_rij*d_rij_ypj;
	  double d_mindi_zpj = d_mindi_rij*d_rij_zpj;
	  // d_Sonmindi_dm
	  double d_Sonmindi_dm = ds_on_dm(1.0,min_dist_list[grid_index],\
					  params.CC_mind,params.deltaCC)* \
	    (1.0/params.deltaCC);
	  // d_Sonmindi over x/y/z
	  double d_Sonmindi_xpj = d_Sonmindi_dm*d_mindi_xpj;
	  //cout << " mindi " << beta/std::log(sum_dist[grid_index]) <<  " d_Sonmindi_dm " <<  d_Sonmindi_dm << " d_mindi_rij " << d_mindi_rij << " d_rij_xpj " << d_rij_xpj << endl;
	  double d_Sonmindi_ypj = d_Sonmindi_dm*d_mindi_ypj;
	  double d_Sonmindi_zpj = d_Sonmindi_dm*d_mindi_zpj;
	  // For exposure calculations we must inspect all neighbors of i
	  vector<int> neighbors_i = neighbors[grid_index];
	  double d_Emin_xpj = 0.0;
	  double d_Emin_ypj = 0.0;
	  double d_Emin_zpj = 0.0;
	  for (unsigned l=0; l < neighbors_i.size(); l++)
	    {
	      int k = neighbors_i[l];
	      double xk = new_x[k];
	      double yk = new_y[k];
	      double zk = new_z[k];
	      double rik_x = xk-xi;
	      double rik_y = yk-yi;
	      double rik_z = zk-zi;
	      double dik2 = rik_x*rik_x + rik_y*rik_y + rik_z*rik_z;
	      double dik = sqrt(dik2);//should be constant for given j,k pair
	      double rjk_x = xk-xj;
	      double rjk_y = yk-yj;
	      double rjk_z = zk-zj;
	      double djk2 = rjk_x*rjk_x + rjk_y*rjk_y + rjk_z*rjk_z;
	      double djk = sqrt(djk2);
	      if (djk < min_modr)
		{
		  djk = min_modr;
		  djk2 = djk*djk;
		}
	      double son_rik = s_on( 1.0, dik, params.GP1_min, params.deltaGP1);//eq 25 SI
	      double soff_rik = s_off( 1.0, dik, params.GP2_min, params.deltaGP2);
	      // Could opt by doing 1 pass to compute all d_mindk_rik and then
	      // look up instead of recomputing as we go through each neighbors
	      double num = (beta*beta)*(exp(beta/djk));
	      double den = djk2*sum_dist[k]*pow(std::log(sum_dist[k]),2);
	      double d_mindk_rjk = num/den;
	      //cout << " d_mindk_rjk " << d_mindk_rjk << " min_dist_list " << min_dist_list[k] << endl;
	      double d_rjk_xpj =  - rjk_x / djk;
	      double d_rjk_ypj =  - rjk_y / djk;
	      double d_rjk_zpj =  - rjk_z / djk;
	      double d_mindk_xpj = d_mindk_rjk * d_rjk_xpj;
	      double d_mindk_ypj = d_mindk_rjk * d_rjk_ypj;
	      double d_mindk_zpj = d_mindk_rjk * d_rjk_zpj;
	      double d_Soffmindk_m = ds_off_dm(1.0,min_dist_list[k],\
					       params.CC2_min,params.deltaCC2);
	      double allterms_x = d_Soffmindk_m*d_mindk_xpj*son_rik*soff_rik;
	      double allterms_y = d_Soffmindk_m*d_mindk_ypj*son_rik*soff_rik;
	      double allterms_z = d_Soffmindk_m*d_mindk_zpj*son_rik*soff_rik;
	      //cout << " l " << l << " djk " << djk << " d_Soffmindk_m " << d_Soffmindk_m << " d_mindk_xpj " << d_mindk_xpj << " son_rik " << son_rik << " soff_rik " << soff_rik << endl;
	      //cout << " allterms_x " << allterms_x << endl;
	      //exit(0);
	      d_Emin_xpj += allterms_x;
	      d_Emin_ypj += allterms_y;
	      d_Emin_zpj += allterms_z;
	    }
	  d_Emin_xpj *= (1.0/params.deltaCC2);
	  d_Emin_ypj *= (1.0/params.deltaCC2);
	  d_Emin_zpj *= (1.0/params.deltaCC2);
	  // d_Sonexposurei_m
	  double d_Sonexposurei_m = ds_on_dm(1.0, exposure[grid_index],\
					     params.Emin, params.deltaE);
	  double d_Sonexposurei_xpj = d_Sonexposurei_m*(1.0/params.deltaE)\
	    *d_Emin_xpj;
	  //cout << "  d_Emin_xpj" <<  d_Emin_xpj << endl;
	  double d_Sonexposurei_ypj = d_Sonexposurei_m*(1.0/params.deltaE)\
	    *d_Emin_ypj;
	  double d_Sonexposurei_zpj = d_Sonexposurei_m*(1.0/params.deltaE)\
	    *d_Emin_zpj;;
	  // d_ai_xpj
	  double term1 = s_on_exposure[grid_index] * d_Sonmindi_xpj;
	  double term2 = s_on_mind[grid_index] * d_Sonexposurei_xpj;
	  //cout << "  d_Sonexposurei_xpj" <<  d_Sonexposurei_xpj << endl;
	  double d_ai_xpj = grid_s_off_bsi[grid_index] * term1 * term2;
	  //cout << " grid_s_off_bsi[grid_index] " << grid_s_off_bsi[grid_index] << " term1 " << term1 << " term2 " << term2 << endl;
	  // d_ai_ypj
	  term1 = s_on_exposure[grid_index] * d_Sonmindi_ypj;
	  term2 = s_on_mind[grid_index] * d_Sonexposurei_ypj;
	  double d_ai_ypj = grid_s_off_bsi[grid_index] * term1 * term2;
	  // d_ai_zpj
	  term1 = s_on_exposure[grid_index] * d_Sonmindi_zpj;
	  term2 = s_on_mind[grid_index] * d_Sonexposurei_zpj;
	  double d_ai_zpj = grid_s_off_bsi[grid_index] * term1 * term2;
	  //cout << " j atom " << j << " active grid i " << i << " activity " << activity[grid_index] << " d_ij " << dij << " d_Sonmindi_xpj " <<  d_Sonmindi_xpj << " y " << d_Sonmindi_ypj << " z " << d_Sonmindi_zpj << " d_Emin_xpj " << d_Emin_xpj << " d_Sonexposurei_xpj  " << d_Sonexposurei_xpj << endl;
	  /*cout << " j atom " << j << " i grid " << i << " activity "	\
	       << activity[grid_index] << " dij " << dij << " d_ai_xpj " \
	       << d_ai_xpj << " d_ai_ypj " << d_ai_ypj		\
	       << " d_ai_zpj " << d_ai_zpj				\
	       << " d_Sonexposurei_m " << d_Sonexposurei_m		\
	       << " Exposure " << exposure[grid_index]		\
	       << " d_Emin_xpj " << d_Emin_xpj \
	       << " s_on_exposure " << s_on_exposure[grid_index] \
	       << " d_Sonexposurei_xpj " << d_Sonexposurei_xpj \
	       << " s_on_mindi " << s_on_mind[grid_index] \
	       << " d_Sonmindi_xpj " << d_Sonmindi_xpj \
	       << endl;*/
	  //exit(0);
	  sum_d_ai_xpj+= d_ai_xpj;
	  sum_d_ai_ypj+= d_ai_ypj;
	  sum_d_ai_zpj+= d_ai_zpj;
	  d_ai_xpj_vec[i] = d_ai_xpj;
	  d_ai_ypj_vec[i] = d_ai_ypj;
	  d_ai_zpj_vec[i] = d_ai_zpj;
	  // Also accumulate now some partial derivatives we need for Ha later
	  //double d_Soffrij_m = ds_off_dm(activity[grid_index], dij, params.r_hydro, params.deltar_hydro);
	  //double d_Soffrij_ai = ds_off_dk(activity[grid_index],dij,params.r_hydro,params.deltar_hydro);
	  //double d_apolarpolari_xpj = d_Soffrij_m*(1.0/params.r_hydro)*d_rij_xpj+ \
	  //      d_Soffrij_ai*d_ai_xpj;
	  //double d_apolarpolari_ypj = d_Soffrij_m*(1.0/params.r_hydro)*d_rij_ypj+ \
	  //   d_Soffrij_ai*d_ai_ypj;
	  //double d_apolarpolari_zpj = d_Soffrij_m*(1.0/params.r_hydro)*d_rij_zpj+ \
	  //  d_Soffrij_ai*d_ai_zpj;
	  //cout << " d_apolarpolari_xpj " << d_apolarpolari_xpj << " d_apolarpolari_ypj " << d_apolarpolari_ypj << " d_apolarpolari_zpj " << d_apolarpolari_zpj << endl;
	  //double num_x;
	  //double num_y;
	  //double num_z;
	  //if (j < n_apolar)
	  //  {
	  //    num_x=d_apolarpolari_xpj*(apolarity[grid_index]+polarity[grid_index]) \
		  //	-apolarity[grid_index]*d_apolarpolari_xpj;
	  //num_y=d_apolarpolari_ypj*(apolarity[grid_index]+polarity[grid_index]) \
	      //	-apolarity[grid_index]*d_apolarpolari_ypj;
	      //num_z=d_apolarpolari_zpj*(apolarity[grid_index]+polarity[grid_index]) \
		  //-apolarity[grid_index]*d_apolarpolari_zpj;
	  //}
	  //else
	  //  {
	  //  num_x=-d_apolarpolari_xpj*(apolarity[grid_index]);
	  //  num_y=-d_apolarpolari_ypj*(apolarity[grid_index]);
	  //  num_z=-d_apolarpolari_zpj*(apolarity[grid_index]);
	  //}
	  //double den2=(apolarity[grid_index]+polarity[grid_index])*(apolarity[grid_index]+polarity[grid_index]);
	  //double d_Hi_xpj=num_x/den2;
	  //double d_Hi_ypj=num_y/den2;
	  //double d_Hi_zpj=num_z/den2;
	  //cout << " d_Hi_xpj " << d_Hi_xpj << " d_Hi_ypj " << d_Hi_ypj << " d_Hi_zpj " << d_Hi_zpj << endl;
	  //sum_d_Hi_xpj += d_Hi_xpj;
	  //sum_d_Hi_ypj += d_Hi_ypj;
	  //sum_d_Hi_zpj += d_Hi_zpj;
	}
      // Compute Vdrug_like derivatives
      double d_V_xpj = sum_d_ai_xpj*Vg;
      double d_V_ypj = sum_d_ai_ypj*Vg;
      double d_V_zpj = sum_d_ai_zpj*Vg;
      double d_sonV_m = ds_on_dm(1.0,volume,params.V_min,params.deltaV_min);
      double d_sonV_xpj = d_sonV_m*(1.0/params.deltaV_min)*d_V_xpj;
      double d_sonV_ypj = d_sonV_m*(1.0/params.deltaV_min)*d_V_ypj;
      double d_sonV_zpj = d_sonV_m*(1.0/params.deltaV_min)*d_V_zpj;
      double d_soffV_m = ds_off_dm(1.0,volume,params.V_max,params.deltaV_max);
      double d_soffV_xpj = d_soffV_m*(1.0/params.deltaV_max)*d_V_xpj;
      double d_soffV_ypj = d_soffV_m*(1.0/params.deltaV_max)*d_V_ypj;
      double d_soffV_zpj = d_soffV_m*(1.0/params.deltaV_max)*d_V_zpj;

      double d_Vdruglike_xpj = s_off_V * d_sonV_xpj + s_on_V * d_soffV_xpj;
      double d_Vdruglike_ypj = s_off_V * d_sonV_ypj + s_on_V * d_soffV_ypj;
      double d_Vdruglike_zpj = s_off_V * d_sonV_zpj + s_on_V * d_soffV_zpj;
      // Compute Va derivatives
      double d_Va_xpj = d_V_xpj/params.V_max;
      double d_Va_ypj = d_V_ypj/params.V_max;
      double d_Va_zpj = d_V_zpj/params.V_max;
      // Compute Ha derivatives
      double sum_d_Hderiv_xpj=0.0;
      double sum_d_Hderiv_ypj=0.0;
      double sum_d_Hderiv_zpj=0.0;
      for (unsigned i=0; i < n_active_grid; i++)
	{
	  int grid_index = active_grid[i];

	  double d_apolari_xpj;
	  double d_polari_xpj;
	  double d_apolari_ypj;
	  double d_polari_ypj;
	  double d_apolari_zpj;
	  double d_polari_zpj;

	  // d_Hi_xpj
	  double dij=dij_vec[i];
	  double d_rij_xpj=-(new_x[grid_index]-xj)/dij;//Or store during previous loop?	  
	  double d_rij_ypj=-(new_y[grid_index]-yj)/dij;
	  double d_rij_zpj=-(new_z[grid_index]-zj)/dij;
	  double ai=activity[grid_index];
	  double d_Soffrij_m=ds_off_dm(ai,dij,params.r_hydro,params.deltar_hydro);
	  double d_Soffrij_ai=ds_off_dk(ai,dij,params.r_hydro,params.deltar_hydro);
	  double deriv_x=d_Soffrij_m*(1/params.r_hydro)*d_rij_xpj+d_Soffrij_ai \
	    *d_ai_xpj_vec[i];
	  double deriv_y=d_Soffrij_m*(1/params.r_hydro)*d_rij_ypj+d_Soffrij_ai\
	    *d_ai_ypj_vec[i];
	  double deriv_z=d_Soffrij_m*(1/params.r_hydro)*d_rij_zpj+d_Soffrij_ai\
	    *d_ai_zpj_vec[i];

	  if (j < n_apolar)
	    {
	      d_polari_xpj=0.0;
	      d_polari_ypj=0.0;
	      d_polari_zpj=0.0;
	      d_apolari_xpj=deriv_x;
	      d_apolari_ypj=deriv_y;
	      d_apolari_zpj=deriv_z;
	    }
	  else
	    {
	      d_polari_xpj=deriv_x;
	      d_polari_ypj=deriv_y;
	      d_polari_zpj=deriv_z;
	      d_apolari_xpj=0.0;
	      d_apolari_ypj=0.0;
	      d_apolari_zpj=0.0;
	    }

	  double apolar=apolarity[grid_index];
	  double term1=(apolar+polarity[grid_index]);
	  double num_x=term1*d_apolari_xpj-\
	    apolar*(d_apolari_xpj+d_polari_xpj);
	  double num_y=term1*d_apolari_ypj-\
	    apolar*(d_apolari_ypj+d_polari_ypj);
	  double num_z=term1*d_apolari_zpj-\
	    apolar*(d_apolari_zpj+d_polari_zpj);
	  double den=term1*term1;
	  double d_Hi_xpj=num_x/den;
	  double d_Hi_ypj=num_y/den;
	  double d_Hi_zpj=num_z/den;
	  //cout << " d_Hi_xpj " << d_Hi_xpj << " d_Hi_ypj " << d_Hi_ypj << " d_Hi_zpj " << d_Hi_zpj << endl;
	  // d_Hderiv_xpj
	  double Hi=hydrophobicity_list[grid_index];
	  double term2_x=Hi*d_ai_xpj_vec[i];
	  double term2_y=Hi*d_ai_ypj_vec[i];
	  double term2_z=Hi*d_ai_zpj_vec[i];
	  double term3_x=ai*d_Hi_xpj;
	  double term3_y=ai*d_Hi_ypj;
	  double term3_z=ai*d_Hi_zpj;
	  double term4=Hi*ai;
	  double num2_x=sum_activity*(term2_x+term3_x)-term4*sum_d_ai_xpj;
	  double num2_y=sum_activity*(term2_y+term3_y)-term4*sum_d_ai_ypj;
	  double num2_z=sum_activity*(term2_z+term3_z)-term4*sum_d_ai_zpj;
	  double d_Hderiv_xpj=num2_x/sum_activity2;
	  double d_Hderiv_ypj=num2_y/sum_activity2;
	  double d_Hderiv_zpj=num2_z/sum_activity2;
	  //cout << " j " << j << " i " << i << " d_Hderiv_xpj " << d_Hderiv_xpj << endl;
	  sum_d_Hderiv_xpj += d_Hderiv_xpj;
	  sum_d_Hderiv_ypj += d_Hderiv_ypj;
	  sum_d_Hderiv_zpj += d_Hderiv_zpj;
	}
      //cout << " sum_d_Hi_xpj " << sum_d_Hi_xpj << " sum_d_Hi_ypj " << sum_d_Hi_ypj << " sum_d_Hi_zpj " << sum_d_Hi_zpj << endl;
      //double d_Ha_xpj = sum_d_Hi_xpj/params.V_max;
      //double d_Ha_ypj = sum_d_Hi_ypj/params.V_max;
      //double d_Ha_zpj = sum_d_Hi_zpj/params.V_max;
      double d_Ha_xpj = sum_d_Hderiv_xpj;
      double d_Ha_ypj = sum_d_Hderiv_ypj;
      double d_Ha_zpj = sum_d_Hderiv_zpj;
      //exit(0);
      // Compute JEDI derivatives
      double term1_x=0.0;
      double term1_y=0.0;
      double term1_z=0.0;
      if (Vdrug_like>0)
	{
	  term1_x=(1/Vdrug_like)*d_Vdruglike_xpj;
	  term1_y=(1/Vdrug_like)*d_Vdruglike_ypj;
	  term1_z=(1/Vdrug_like)*d_Vdruglike_zpj;
	}
      double term2=(1/(params.alpha*Va+params.beta*Ha+params.gamma));
      double term3_x=(params.alpha*d_Va_xpj+params.beta*d_Ha_xpj);
      double term3_y=(params.alpha*d_Va_xpj+params.beta*d_Ha_ypj);
      double term3_z=(params.alpha*d_Va_xpj+params.beta*d_Ha_zpj);
      double d_Jedi_xpj=Jedi*(term1_x+term2*term3_x);
      double d_Jedi_ypj=Jedi*(term1_y+term2*term3_y);
      double d_Jedi_zpj=Jedi*(term1_z+term2*term3_z);
      //cout << " *** atom j *** " << j << endl;
      //cout << " d_V_xpj " << d_V_xpj << " d_V_ypj " << d_V_ypj << " d_V_zpj " << d_V_zpj <<endl;
      //cout << " d_Vdruglike_xpj " << d_Vdruglike_xpj << " d_Vdruglike_ypj " << d_Vdruglike_ypj << " d_Vdruglike_zpj " << d_Vdruglike_zpj <<endl; 
      //cout << " d_Va_xpj " << d_Va_xpj << " d_Va_ypj " << d_Va_ypj << " d_Va_zpj " << d_Va_zpj <<endl; 
      //cout << " d_Ha_xpj " << d_Ha_xpj << " d_Ha_ypj " << d_Ha_ypj << " d_Ha_zpj " << d_Ha_zpj <<endl; 
      //cout << "atom j " << j << " d_Jedi_xpj " << d_Jedi_xpj << " d_Jedi_ypj " << d_Jedi_ypj << " d_Jedi_zpj " << d_Jedi_zpj << endl;
      setAtomsDerivatives(j,Vector(d_Jedi_xpj,d_Jedi_ypj,d_Jedi_zpj));
      //exit(0);
    }

  //TODO: Flag to output infrequently detailed grid statistics
  // dx files with updated position and 'activitiy'
  //                                and 'Ha'

  //exit(0);

}//close jedi::calculate 
}//close namespace colvar
}//cloase namespace PLMD
