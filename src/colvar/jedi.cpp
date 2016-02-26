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
// Lapack needed for l2-mininum norm solution
#include "../tools/lapack/lapack.h"

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
  int grid_extent[3];//the number of grid points along each axis
  int grid_origin_idx;//the index of the grid point at the origin of the grid
  vector<double> grid_s_off_bsi;//binding site score of grid point (eq 5 term 1 Cuchillo et al. JCTC 2015)
  jediparameters params;// parameters druggability estimator
  vector<vector<int> > neighbors;//list of grid indices that are neighbors of a given grid point
  string summary_file;//path to output file
  int stride;//frequency of output (in timesteps) to summary file;
  string gridstats_folder;//path to output grid folder;
  int gridstride;//frequency of output (in timestep) of a grid file;

  double refsite_com[3];//reference coordinates of the center of mass of the binding site region
  double grid_ref_cog[3];//reference coordinates of the center of geometry of the grid

  //double b_grid;//atom number of the first grid point
  double n_grid;// total number of grid points (double)
  //double cutoff_close;//cutoff close contact
  //double cutoff_far;//cutoff distant contact
  //double cutoff_enclosure;//cutoff NP for distant contact
  //double cutoff_hull;//cutoff hull
  //double cutoff_surface;//cutoff surface
  //double cutoff_contact;//cutoff contact for surface
  //double cutoff_hydro;//cutoff hydrophobicity
  //double cutoff_con;//cutoff connectivity
  //double dump_matrix;//paramater to control the output of jedi (just to check, could be set to 0)
  double delta;
  //string line;
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
  keys.add("optional", "SITE","a file listing coordinates of atoms used to define a binding site region.");
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
  gridstride=atoi(gridstride_string.c_str());
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
  //vector<AtomNumber> Apolar;
  //cout << " Apolar has ? elements " << Apolar.size() << endl;

  PDB apolar_pdb;
  if( !apolar_pdb.read(apolar_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + apolar_file );

  apolaratoms = apolar_pdb.getAtomNumbers();
  const std::vector<Vector> apolar_positions = apolar_pdb.getPositions();
  const std::vector<double> apolar_masses = apolar_pdb.getOccupancy();
  cout << " apolaratoms has ? elements " << apolaratoms.size() << endl;  

  //Polar 
  //vector<AtomNumber> Polar;
  //cout << " Polar has ? elements " << Polar.size() << endl;

  PDB polar_pdb;
  if( !polar_pdb.read(polar_file,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + polar_file );
  polaratoms = polar_pdb.getAtomNumbers();
  const std::vector<Vector> polar_positions = polar_pdb.getPositions();
  const std::vector<double> polar_masses = polar_pdb.getOccupancy();
  cout << " polaratoms has ? elements " << polaratoms.size() << endl;

  // Also compute the com of reference coordinates and save for future calcs
  double site_mass_tot=0.0;
  refsite_com[0] = 0.0;
  refsite_com[1] = 0.0;
  refsite_com[2] = 0.0;
  for (unsigned i=0; i < apolaratoms.size() ; ++i)
    {
      refsite_com[0] += apolar_masses[i] * apolar_positions[i][0];
      refsite_com[1] += apolar_masses[i] * apolar_positions[i][1];
      refsite_com[2] += apolar_masses[i] * apolar_positions[i][2];
      ref_pos.push_back( Vector(apolar_positions[i][0], apolar_positions[i][1], apolar_positions[i][2]) );
      site_mass_tot += apolar_masses[i];
    }
  for (unsigned i=0; i < polaratoms.size() ; ++i)
    {
      refsite_com[0] += polar_masses[i] * polar_positions[i][0];
      refsite_com[1] += polar_masses[i] * polar_positions[i][1];
      refsite_com[2] += polar_masses[i] * polar_positions[i][2];
      ref_pos.push_back( Vector(polar_positions[i][0], polar_positions[i][1], polar_positions[i][2]) );
      site_mass_tot += polar_masses[i];
    }
  refsite_com[0] /= site_mass_tot;
  refsite_com[1] /= site_mass_tot;
  refsite_com[2] /= site_mass_tot;

  cout << " refsite_com " << refsite_com[0] << " " << refsite_com[1] << " " << refsite_com[2] << endl;
  //exit(0);

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
      //ref_pos.push_back( Vector(alignment_positions[i][0], alignment_positions[i][1], alignment_positions[i][2]) );
      ref_mass_tot += alignment_masses[i];
    }
  ref_com[0] /= ref_mass_tot;
  ref_com[1] /= ref_mass_tot;
  ref_com[2] /= ref_mass_tot;

  cout << " ref_com " << ref_com[0] << " " << ref_com[1] << " " << ref_com[2] << endl;


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
  double gmin_x=999999.0;
  double gmax_x=-99999.0;
  double gmin_y=999999.0;
  double gmax_y=-99999.0;
  double gmin_z=999999.0;
  double gmax_z=-99999.0;  
  grid_origin_idx=0;//Assumed to be the first grid atom. Is that always right?
  for (unsigned i=0; i < grid_positions.size() ; i++)
    {
      double gx = grid_positions[i][0];
      double gy = grid_positions[i][1];
      double gz = grid_positions[i][2];
      if (gx < gmin_x)
	gmin_x = gx;
      if (gx > gmax_x)
	gmax_x = gx;
      if (gy < gmin_y)
	gmin_y = gy;
      if (gy > gmax_y)
	gmax_y = gy;
      if (gz < gmin_z)
	gmin_z = gz;
      if (gz > gmax_z)
	gmax_z = gz;
    }
  int n_gx = int ( (gmax_x-gmin_x)/params.resolution ) + 1;
  int n_gy = int ( (gmax_y-gmin_y)/params.resolution ) + 1;
  int n_gz = int ( (gmax_z-gmin_z)/params.resolution ) + 1;
  cout << "n_g is " << n_gx << " " << n_gy << " " << n_gz << endl;
  grid_extent[0] = n_gx;
  grid_extent[1] = n_gy;
  grid_extent[2] = n_gz;
  // Work out number of grid points along x/y/z components
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
  wfile << "#step \t JEDI \t Vdrug_like \t Va \t Ha \t MaxDerivIdx_raw \t max_deriv_x_raw \t max_deriv_y_raw \t max_deriv_z_raw \t MaxDerivIdx \t max_deriv_x \t max_deriv_y \t max_deriv_z \t rmsd" << endl;
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
      //s=k*(1.-(pow((1.-pow(m,2)),2))*((1.+2*(pow(m,2)))));
      double m2=m*m;
      double m4=m2*m2;
      double m6=m4*m2;
      s=k*(3*m4-2*m6);
      //double m3 = m*m*m;
      //double m4 = m3*m;
      //double m5 = m4*m;
      //s=k*(6*m5-15*m4+10*m3);
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
      //s=k*(pow((1.-pow(m,2)),2)*((1.+2*(pow(m,2)))));
      double m2=m*m;
      double m4=m2*m2;
      double m6=m4*m2;
      s=k*(1-3*m4+2*m6);
      //double m3 = m*m*m;
      //double m4 = m3*m;
      //double m5 = m4*m;
      //s=k*(1-6*m5+15*m4-10*m3);
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
      //double m2=m*m;
      //s = 4*k*m*( (1.-m2)*(1.-m2) - (1.-m2)*(1.+2*m2) );
      double m2=m*m;
      double m3=m2*m;
      double m5=m3*m2;
      s = 12*k*(m5-m3);
      //double m2=m*m;
      //double m3=m2*m;
      //double m4=m2*m2;
      //s = k*(-30*m4+60*m3-30*m2);
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
      //double m2=m*m;
      //s = ( (1-m2)*(1-m2) ) * (1+2*m2);
      double m2=m*m;
      double m4=m2*m2;
      double m6=m4*m2;
      s=1-3*m4+2*m6;
      //double m3 = m*m*m;
      //double m4 = m3*m;
      //double m5 = m4*m;
      //s=1-6*m5+15*m4-10*m3;
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
      //double m2=m*m;
      //s = -4*k*m*( (1.-m2)*(1-m2) - (1.-m2)*(1.+2*m2) );
      double m2=m*m;
      double m3=m2*m;
      double m5=m3*m2;
      s = 12*k*(m3-m5);
      //double m2=m*m;
      //double m3=m2*m;
      //double m4=m2*m2;
      //s = k*(30*m4-60*m3+30*m2);
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
      //double m2=m*m;
      //s = 1 -(1-m2)*(1-m2)*(1+2*m2);
      double m2=m*m;
      double m4=m2*m2;
      double m6=m4*m2;
      s=(3*m4-2*m6);
      //double m3 = m*m*m;
      //double m4 = m3*m;
      //double m5 = m4*m;
      //s=(6*m5-15*m4+10*m3);
    }
  return s;
 }

// calculator
void jedi::calculate(){

  unsigned size_grid = grid_positions.size();
  //vector<int> active_grid;// array with the index of active grid point (0 if inactive)
  //double activity[size_grid];//array of activity scores for each grid point
  vector<double> activity;
  activity.reserve(size_grid);
  //double s_on_exposure[size_grid];// array of s_on_exposure scores for each grid point
  vector<double> s_on_exposure;
  s_on_exposure.reserve(size_grid);
  //double s_on_mind[size_grid];// array with the penalty of close contact
  vector <double> s_on_mind;
  s_on_mind.reserve(size_grid);
  vector <double> s_on_contacts;
  s_on_contacts.reserve(size_grid);
  //double sum_dist[size_grid];// array with the sum of the distances between each grid point and all atoms in the CV
  vector<double> sum_dist;
  sum_dist.reserve(size_grid);
  //double min_dist_list[size_grid];// array with the minimum distance between each grid point and all atoms in the CV
  vector<double> min_dist_list;
  min_dist_list.reserve(size_grid);
  vector<double> num_contacts;
  num_contacts.reserve(size_grid);  
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
  double min_modr = 0.10;// that's 0.5 Angstrom

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

  // Check the translation of the center of mass of the binding site region

  unsigned n_apolar = apolaratoms.size();
  unsigned n_apolarpolar = n_apolar + polaratoms.size();
  double site_com[3] = {0.0, 0.0 ,0.0};
  double site_mass=0.0;
  for (unsigned j =0; j < n_apolarpolar ; j++)
    {
      double j_mass = getMass(j);
      Vector j_pos = getPosition(j);
      site_com[0] += j_mass * j_pos[0];
      site_com[1] += j_mass * j_pos[1];
      site_com[2] += j_mass * j_pos[2];
      site_mass += j_mass;
    }
  site_com[0] /= site_mass;
  site_com[1] /= site_mass;
  site_com[2] /= site_mass;

  //cout << " site_com " << site_com[0] << " " << site_com[1] << " " << site_com[2] << endl;
  double delta_com[3];
  delta_com[0] = site_com[0] - refsite_com[0];
  delta_com[1] = site_com[1] - refsite_com[1];
  delta_com[2] = site_com[2] - refsite_com[2];
  //cout << " delta_com " << delta_com[0] << " " << delta_com[1] << " " << delta_com[2] << endl;
  //exit(0);

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

  //for (int i=0; i < n_align ; ++i)
  for (int i=0; i < n_apolarpolar ; ++i)
    {
      ref_xlist[i][0] = ref_pos[i][0];//FIXME change calculate_rotation_rmsd args to take directly vector in
      ref_xlist[i][1] = ref_pos[i][1];
      ref_xlist[i][2] = ref_pos[i][2];
      //Vector i_pos = getPosition( n_apolarpolar + i );// MUST CHECK that PBC DO NOT MESS THINGS UP
      Vector i_pos = getPosition( i );
      mov_xlist[i][0] = i_pos[0];
      mov_xlist[i][1] = i_pos[1];
      mov_xlist[i][2] = i_pos[2];
      //cout << " i " << i << " mov_xlist " << mov_xlist[i][0] << " " << mov_xlist[i][1] << " " << mov_xlist[i][2] << endl;
      //cout << " i " << i << " ref_xlist " << ref_xlist[i][0] << " " << ref_xlist[i][1] << " " << ref_xlist[i][2] << endl;      
      // Also get data to compute mov_com
      //double i_mass = getMass( n_apolarpolar + i );//FIXME this is a constant so could be cached
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

  mov_to_ref[0] = ref_com[0] - mov_com[0];//ref_com defined during CV init
  mov_to_ref[1] = ref_com[1] - mov_com[1];
  mov_to_ref[2] = ref_com[2] - mov_com[2];

  //cout << "mov_com " << mov_com[0] << " " << mov_com[1] << " " << mov_com[2] << endl;
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

  rotmat[0][0] = 1.0;
  rotmat[0][1] = 0.0;
  rotmat[0][2] = 0.0;
  rotmat[1][0] = 0.0;
  rotmat[1][1] = 1.0;
  rotmat[1][2] = 0.0;
  rotmat[2][0] = 0.0;
  rotmat[2][1] = 0.0;
  rotmat[2][2] = 1.0;

  //calculate_rotation_rmsd( ref_xlist, mov_xlist, n_align, mov_com, mov_to_ref, rotmat, &rmsd  );
  calculate_rotation_rmsd( ref_xlist, mov_xlist, n_apolarpolar, mov_com, mov_to_ref, rotmat, &rmsd  );

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

  //cout << "Just called calculated_rotation_rmsd" << endl;
  //cout << "rotmat elements :" << endl;
  //cout << rotmat[0][0] << " " << rotmat[0][1] << " " << rotmat[0][2] << endl;
  //cout << rotmat[1][0] << " " << rotmat[1][1] << " " << rotmat[1][2] << endl;
  //cout << rotmat[2][0] << " " << rotmat[2][1] << " " << rotmat[2][2] << endl; 

  //exit(0);
  /*
  double new_grid_cog_x, local_grid_cog_x;
  double new_grid_cog_y, local_grid_cog_y;
  double new_grid_cog_z, local_grid_cog_z;
  // Apply same translation to grid_cog
  new_grid_cog_x = grid_ref_cog[0] - mov_to_ref[0];
  new_grid_cog_y = grid_ref_cog[1] - mov_to_ref[1];
  new_grid_cog_z = grid_ref_cog[2] - mov_to_ref[2];
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
  */



  // Translate grid_cog by delta_com
  double new_grid_cog_x;
  double new_grid_cog_y;
  double new_grid_cog_z;
  //new_grid_cog_x = grid_ref_cog[0] + delta_com[0];
  //new_grid_cog_y = grid_ref_cog[1] + delta_com[1];
  //new_grid_cog_z = grid_ref_cog[2] + delta_com[2];
  // THIS DOESN'T SOUND RIGHT
  // WHAT ABOUT FITTING GRID COM ON REFSITE_COM ?
  // I THINK THIS DOES NOT WORK BECAUSE IT DO NOT ROTATE ABOUT SAME CENTER
  new_grid_cog_x = site_com[0];
  new_grid_cog_y = site_com[1];
  new_grid_cog_z = site_com[2];

  //cout << " New grid cog " << new_grid_cog_x << " " << new_grid_cog_y << " " << new_grid_cog_z << endl;


  // Now rotate all grid points at origin and then translate to new cog
  //double grid_min[3] = {99999.0, 99999.0, 99999.0};//minimum grid coordinates
  for(unsigned i=0;i< size_grid;i++)
    {
      new_x[i]  = ( rotmat[0][0]*(grid_positions[i][0]) + rotmat[1][0]*(grid_positions[i][1]) + rotmat[2][0]*(grid_positions[i][2]) ) + new_grid_cog_x;
      new_y[i]  = ( rotmat[0][1]*(grid_positions[i][0]) + rotmat[1][1]*(grid_positions[i][1]) + rotmat[2][1]*(grid_positions[i][2]) ) + new_grid_cog_y;
      new_z[i]  = ( rotmat[0][2]*(grid_positions[i][0]) + rotmat[1][2]*(grid_positions[i][1]) + rotmat[2][2]*(grid_positions[i][2]) ) + new_grid_cog_z;
      //if (new_x[i] < grid_min[0])
      //	grid_min[0] = new_x[i];
      //if (new_y[i] < grid_min[1])
      //	grid_min[1] = new_y[i];
      //if (new_z[i] < grid_min[2])
      //	grid_min[2] = new_z[i];
      //new_x[i]  = ( rotmat[0][0]*(grid_positions[i][0]) + rotmat[0][1]*(grid_positions[i][1]) + rotmat[0][2]*(grid_positions[i][2]) ) + new_grid_cog_x;
      //new_y[i]  = ( rotmat[1][0]*(grid_positions[i][0]) + rotmat[1][1]*(grid_positions[i][1]) + rotmat[1][2]*(grid_positions[i][2]) ) + new_grid_cog_y;
      //new_z[i]  = ( rotmat[2][0]*(grid_positions[i][0]) + rotmat[2][1]*(grid_positions[i][1]) + rotmat[2][2]*(grid_positions[i][2]) ) + new_grid_cog_z;


      //cout << "grid_positions[gridpoint][0] is " << grid_positions[gridpoint][0] << endl;
      //cout << "grid_positions[gridpoint][1] is " << grid_positions[gridpoint][1] << endl;
      //cout << "grid_positions[gridpoint][2] is " << grid_positions[gridpoint][2] << endl;
      //cout << "new_x , new_y , new_z" << new_x[i] << "," << new_y[i] << "," << new_z[i] << endl; 
    }


  // Periodically update ref_pos and refsite_com with current values
  // Note: do not update too frequently, otherwise fitting algorithm doesn't
  // detect noticeable rotation??
  int step= getStep();
  //int mod = fmod(step,gridstride);
  int mod = fmod(step,1);
  bool iszero = mod;
  if (!iszero)
    {
      //cout << "step " << step << " Updating reference grid positions " << endl;
      //if (step >0)
      //	exit(0);
      for (int i=0; i < n_apolarpolar ; ++i)
	{
	  Vector i_pos = getPosition( i );
	  ref_pos[i][0] = i_pos[0];
	  ref_pos[i][1] = i_pos[1];
	  ref_pos[i][2] = i_pos[2];
	}
      refsite_com[0] = site_com[0];
      refsite_com[1] = site_com[1];
      refsite_com[2] = site_com[2];
      grid_ref_cog[0] = new_grid_cog_x;
      grid_ref_cog[1] = new_grid_cog_y;
      grid_ref_cog[2] = new_grid_cog_z;
    }
  //exit(0);
  //cout << "*** Getting ready for STEP 2" << endl;

  //-----------------------------------------
  //-----------------------------------------
  //
  //       STEP2  :  JEDI score
  //
  //-----------------------------------------
  //-----------------------------------------
  //cout << " Starting Step 2 where JEDI score is calculated " << endl;
  // Hard cutoff
  double cutoff=0.6;
  double cutoff2=cutoff*cutoff;
  double volume=0.0;
  double hydrophobicity_tot=0.0;
  //----------> Compute activity of grid points and also VOLUME
  //unsigned n_apolar = apolaratoms.size();
  //unsigned n_apolarpolar = n_apolar + polaratoms.size();
  for(unsigned i = 0; i < size_grid ; i++)
    {
      double sum = 0.;
      double grd_x = new_x[i];
      double grd_y = new_y[i];
      double grd_z = new_z[i];
      double apolar=0.0;
      double polar=0.0;

      for( unsigned j = 0; j < n_apolarpolar; j++)
	{
	  double rij[3];
	  rij[0] = grd_x - getPosition(j)[0];
	  rij[1] = grd_y - getPosition(j)[1];
	  rij[2] = grd_z - getPosition(j)[2];
	  // FIXME No PBC check !!  Code may need changes for PBC
	  // FIXME avoid sqrt if possible
	  double mod_rij2   = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
	  if (mod_rij2 > cutoff2)
	    continue;
	  double mod_rij= sqrt(mod_rij2);
	  if (mod_rij < min_modr)
	      mod_rij = min_modr;
	  //cout << "grd_x " << grd_x << " y " << grd_y << " z " << grd_z << " jx " << getPosition(j)[0] << " jy " << getPosition(j)[1] << " jz " << getPosition(j)[2] << endl;
//      cout << " i " << i << " j " << j << " mod_rij " << mod_rij << endl;
	  // calculate mindist between grid points and protein atoms (EQUATION 5 2nd term)
	  //I don't get what we< do with ncoord also get number of neighbors in the same pass ...(EQUATION 7 1st term)
	  //ncoord += s_off( 1.0, mod_rij, cutoff_far, D_far);
	  double contact = s_off( 1.0, mod_rij, params.r_hydro, params.deltar_hydro);
	  if (j < n_apolar)
	    apolar += contact;
	  else
	    polar += contact;
	  sum += exp(beta/mod_rij);
	}
      if (sum < 0.0001)
	{
	  sum=exp(beta/cutoff);
	}
      sum_dist[i] = sum;
      double min_dist = beta/std::log(sum);
      min_dist_list[i] = min_dist;
      s_on_mind[i] = s_on( 1.0, min_dist, params.CC_mind, params.deltaCC);
      apolarity[i] = apolar;
      polarity[i] = polar;
      double ncontacts = apolar+polar;
      num_contacts[i] = ncontacts;
      //s_on_contacts[i] = s_on( 1.0, ncontacts, params.Emin, params.deltaE);
      //cout << " i " << i << " mind_dist " << min_dist << " " << params.CC_mind << " " << params.deltaCC << endl;
      //activity[i] = s_on_mind[i] * s_on_contacts[i] * grid_s_off_bsi[i];
      if( ncontacts > 0.)
	hydrophobicity_list[i]=apolar/(ncontacts);
      else
	hydrophobicity_list[i]=0.;
      // Only grid points with activity > 0 are kept for the next calculations
      //if (activity[i] > 0.0)
      //	active_grid.push_back(i);
      //volume += activity[i];
      //hydrophobicity_tot+=hydrophobicity_list[i]*activity[i];
      //cout << "i " << i << " activity[i] " << activity[i] << " mindi " << min_dist << " s_on_mind[i] " << s_on_mind[i] << " ncontacts " << ncontacts << " s_on_contacts[i] "<< s_on_contact<s[i] << " grid_s_off_bsi[i] " << grid_s_off_bsi[i]  << " apolar " << apolar << " polar " << polar << " hydrophobicity " << apolar/ncontacts << endl;
    }
  //cout << " volume " << volume << " hydrophobicity " << hydrophobicity_tot;
  //FIXME: Comment what is going on here
  // For each grid point...
  for( int i = 0; i < size_grid; i++)
    {
      double exposure_score = 0.0;
      vector<int> neighbors_i = neighbors[i];
      for(unsigned j = 0; j < neighbors_i.size(); j++)
  	{
  	  int k = neighbors_i[j];
 	  // This is equation 7 of the paper so exposure_score is s_on_exposure_i
	  // OPTM term two and three should be constant and precomputed !!
	  //cout << " k " << k << " min_dist_list[k] " << min_dist_list[k] << endl;
	  double term1 = s_off(1.0,min_dist_list[k],params.CC2_min,params.deltaCC2);
	  //double term2 = s_on(1.0,mod_rij,params.GP1_min,params.deltaGP1);
	  //double term3 = s_off(1.0,mod_rij,params.GP2_min,params.deltaGP2);
	  //cout << " k " << k << " term1 " << term1 << " term2 " << term2 << " term3 " << term3 << endl;
	  //double term2 = 1.0;
	  //double term3 = 1.0;
	  exposure_score += term1;
	}
      exposure[i] = exposure_score;
      s_on_exposure[i]  = s_on( 1.0, exposure_score, params.Emin, params.deltaE);
      //cout << " i " << i << " neighbors " << neighbors_i.size() << " exposure " \
      // 	   << exposure_score << " s_on_exposure " << s_on_exposure[i] << endl;
      //s_on_exposure[i] = 1.0;
      // This now gives the activity value from equation 5
      activity[i] = s_on_mind[i] * s_on_exposure[i] * grid_s_off_bsi[i];

      // Only grid points with activity > 0 are kept for the next calculations
      //if (activity[i] > 0.0)
      //	active_grid.push_back(i);
      volume += activity[i];
      hydrophobicity_tot+=hydrophobicity_list[i]*activity[i];
      //cout << "i " << i << " activity[i] " << activity[i] << " s_on_mind[i] " << s_on_mind[i] << " s_on_exposure[i] "<< s_on_exposure[i] << " grid_s_off_bsi[i] " << grid_s_off_bsi[i]  << " exposure[i] " << exposure[i] << " volume " <<  volume << endl;
    }
  
  //cout << "There are " << active_grid.size() << " active grid points " << endl;
  double sum_activity = volume;
  double sum_activity2 = sum_activity*sum_activity;
  volume *= Vg;//Equation 4 of the paper
  // ----------> "Drug-like" volume
  double s_off_V = (s_off( 1.0, volume, params.V_max, params.deltaV_max));
  //double s1 = s_off_V;
  //double s2 = (s_on( 1.0, volume, params.V_min, params.deltaV_min));
  double s_on_V = s_on( 1.0, volume, params.V_min, params.deltaV_min);
  double Vdrug_like = s_off_V * s_on_V;
  //cout << " volume " << volume << " Vdrug_like " << Vdrug_like << endl;

  //----------> Hydrophobicity
  /*double hydrophobicity_tot=0.0;
  for( unsigned i = 0; i < active_grid.size() ; i++)
    {
      double apolar = 0.;
      for( unsigned j = 0; j < apolaratoms.size() ; j++)
	{
	  double rij[3];
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
	  double rij[3];
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
      }*/

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
  double Jedi=Vdrug_like*(params.alpha*Va + params.beta*Ha + params.gamma);
  //double Jedi=Vdrug_like*(params.alpha*Va);
  //double Jedi=(params.alpha*Va);
  setValue(Jedi);

  cout << "Jedi score is " << Jedi <<  " sum_ai " << sum_activity << " Va " << Va << " Ha " << Ha << endl;
  //exit(0);

  // cout.precision(9);
  // cout << "current, Jedi, Vdrug_like, Ha, volume/Vmax, COM_x, COM_y, COM_z, score[0], score[1], score[2], score[3], score[4], score[5], score[6], score[7], score[8]" << endl;
  // cout << current << " " << Jedi << " " << Vdrug_like << " " << Ha << " " << volume/Vmax
  // << " " << COM_x << " " << COM_y << " " << COM_z << " " << score[0] << " " << score[1] << " "
  // << score[2] << " " << score[3] << " " << score[4] << " " << score[5] << " " << score[6] << " " << score[7] << " " << score[8] <<endl;

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
  
  //cout << "@@@ Doing derivatives step " << step << endl;

  //unsigned n_apolar = apolaratoms.size();
  unsigned n_polar = polaratoms.size();
  unsigned n_cv_atoms = n_apolar + n_polar;
  //unsigned n_active_grid = active_grid.size();

  vector<double> d_ai_xpj_vec;
  //d_ai_xpj_vec.reserve(n_active_grid);
  d_ai_xpj_vec.reserve(size_grid);
  vector<double> d_ai_ypj_vec;
  //d_ai_ypj_vec.reserve(n_active_grid);
  d_ai_ypj_vec.reserve(size_grid);
  vector<double> d_ai_zpj_vec;
  //d_ai_zpj_vec.reserve(n_active_grid);
  d_ai_zpj_vec.reserve(size_grid);
  vector<double> dij_vec;
  //dij_vec.reserve(n_active_grid);
  dij_vec.reserve(size_grid);

  vector<double> d_Jedi_xpj_vec;
  d_Jedi_xpj_vec.reserve(n_cv_atoms);
  vector<double> d_Jedi_ypj_vec;
  d_Jedi_ypj_vec.reserve(n_cv_atoms);
  vector<double> d_Jedi_zpj_vec;
  d_Jedi_zpj_vec.reserve(n_cv_atoms);

  double sum_d_Jedi_xpj =0.0;
  double sum_d_Jedi_ypj =0.0;
  double sum_d_Jedi_zpj =0.0;
  double sum_d_Jedi_torque_xpj=0.0;
  double sum_d_Jedi_torque_ypj=0.0;
  double sum_d_Jedi_torque_zpj=0.0;

  // Update derivative every x steps
  mod = step % 1;
  iszero = mod;

  for ( unsigned j=0; j < n_cv_atoms ; j++)
    {
      double xj = getPosition(j)[0];
      double yj = getPosition(j)[1];
      double zj = getPosition(j)[2];
      // Compute activity derivatives (should we include pts with 0 activity?)
      double sum_d_ai_xpj=0.0;
      double sum_d_ai_ypj=0.0;
      double sum_d_ai_zpj=0.0;
      //for (unsigned int i=0; i < n_active_grid ; i++)
      for (unsigned i=0; i < size_grid ; i++)
	{
	  //int grid_index = active_grid[i];
	  int grid_index=i;
	  d_ai_xpj_vec[i] = 0.0;
	  d_ai_ypj_vec[i] = 0.0;
	  d_ai_zpj_vec[i] = 0.0;
	  double xi = new_x[grid_index];
	  double yi = new_y[grid_index];
	  double zi = new_z[grid_index];
	  // d_rij over x/y/z
 	  //JM 25/02 direction of vector?
	  //double dij_x = xi-xj;
	  //double dij_y = yi-yj;
	  //double dij_z = zi-zj;
	  double dij_x = xj-xi;
	  double dij_y = yj-yi;
	  double dij_z = zj-zi;
	  //cout << " dij_x " << dij_x << endl;
	  //cout << " dij_x " << dij_x << " " << dij_y << " " << dij_z << endl; 
	  //cout << " grid index " << grid_index << endl;
	  //cout << " xi " << xi << " yi " << yi << " zi " << zi << endl;
	  //cout << " xj " << xj << " yj " << yj << " zj " << zj << endl;
	  double dij2 = dij_x*dij_x + dij_y*dij_y + dij_z*dij_z;
	  //cutoff
	  if (dij2 > cutoff2)
	    {
	      dij_vec[i] = cutoff;
	      continue;
	    }
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
	  double d_rij_xpj = dij_x / dij;
	  double d_rij_ypj = dij_y / dij;
	  double d_rij_zpj = dij_z / dij;
	  //cout << " d_rij_xpj " << d_rij_xpj << endl;
	  // d_mindi over d_rij
	  double num = (beta*beta)*(exp(beta/dij));
	  // sum_dist was saved during calculation of activities in step 2
	  double den = dij2*sum_dist[grid_index]*pow(std::log(sum_dist[grid_index]),2);
	  double d_mindi_rij = num/den;
	  //cout << "sum_dist " << sum_dist[grid_index];
	  //cout << " num " << num;
	  //cout << " den " << den;
	  //cout << " d_mindi_rij " << d_mindi_rij << endl;
	  // d_mindi over x/y/z
	  double d_mindi_xpj = d_mindi_rij*d_rij_xpj;
	  //cout << " d_mindi_xpj " << d_mindi_xpj << endl;
	  double d_mindi_ypj = d_mindi_rij*d_rij_ypj;
	  double d_mindi_zpj = d_mindi_rij*d_rij_zpj;
	  // d_Sonmindi_dm
	  double d_Sonmindi_dm = ds_on_dm(1.0,min_dist_list[grid_index],\
					  params.CC_mind,params.deltaCC)* \
	    (1.0/params.deltaCC);
	  //cout << "d_Sonmindi_dm*(1/deltaCC) " << d_Sonmindi_dm << endl;
	  // d_Sonmindi over x/y/z
	  double d_Sonmindi_xpj = d_Sonmindi_dm*d_mindi_xpj;
	  //cout << "d_Sonmindi_xpj " << d_Sonmindi_xpj << endl;
	  //cout << " mindi " << beta/std::log(sum_dist[grid_index]) <<  " d_Sonmindi_dm " <<  d_Sonmindi_dm << " d_mindi_rij " << d_mindi_rij << " d_rij_xpj " << d_rij_xpj << endl;
	  double d_Sonmindi_ypj = d_Sonmindi_dm*d_mindi_ypj;
	  double d_Sonmindi_zpj = d_Sonmindi_dm*d_mindi_zpj;

	  //double d_Soncontactsi_dm = ds_on_dm(1.0,num_contacts[grid_index], \
	  //					  params.Emin,params.deltaE)* \
								      //  (1.0/params.deltaE);
	  //cout << "d_Soncontactsi_dm " << d_Soncontactsi_dm << endl;
	  // d_Soncontactsi over x/y/z
	  //double d_Soncontactsi_xpj = d_Soncontactsi_dm*d_mindi_xpj;
	  //cout << "d_Soncontactsi_xpj " << d_Soncontactsi_xpj << endl;
	  //double d_Soncontactsi_ypj = d_Soncontactsi_dm*d_mindi_ypj;
	  //double d_Soncontactsi_zpj = d_Soncontactsi_dm*d_mindi_zpj;

	  // For exposure calculations we must inspect all neighbors of i
	  vector<int> neighbors_i = neighbors[grid_index];
	  double d_exposure_xpj = 0.0;
	  double d_exposure_ypj = 0.0;
	  double d_exposure_zpj = 0.0;
	  for (unsigned l=0; l < neighbors_i.size(); l++)
	    //for (unsigned l=0; l < 0; l++)
	    {
	      int k = neighbors_i[l];
	      double xk = new_x[k];
	      double yk = new_y[k];
	      double zk = new_z[k];
	      double rjk_x = xj-xk;
	      double rjk_y = yj-yk;
	      double rjk_z = zj-zk;
	      double djk2 = rjk_x*rjk_x + rjk_y*rjk_y + rjk_z*rjk_z;
	      double djk = sqrt(djk2);
	      if (djk < min_modr)
		{
		  //cout << " step " << step << " close contact happened i " << i << " j " << j <<  " k " << k << " djk " << djk << endl;
		  djk = min_modr;
		  djk2 = djk*djk;
		}
	      // Could opt by doing 1 pass to compute all d_mindk_rjk and then
	      // look up instead of recomputing as we go through each neighbors
	      double num = (beta*beta)*(exp(beta/djk));
	      double den = djk2*sum_dist[k]*pow(std::log(sum_dist[k]),2);
	      double d_mindk_rjk = num/den;
	      //cout << " k " << k << " d_mindk_rjk " << d_mindk_rjk << " min_dk " << min_dist_list[k] << endl;
	      double d_rjk_xpj =  rjk_x / djk;
	      double d_rjk_ypj =  rjk_y / djk;
	      double d_rjk_zpj =  rjk_z / djk;
	      double d_mindk_xpj = d_mindk_rjk * d_rjk_xpj;
	      double d_mindk_ypj = d_mindk_rjk * d_rjk_ypj;
	      double d_mindk_zpj = d_mindk_rjk * d_rjk_zpj;
	      double d_Soffmindk_m = ds_off_dm(1.0,min_dist_list[k],\
					       params.CC2_min,params.deltaCC2);
	      double allterms_x = d_Soffmindk_m*d_mindk_xpj;
	      double allterms_y = d_Soffmindk_m*d_mindk_ypj;
	      double allterms_z = d_Soffmindk_m*d_mindk_zpj;
	      //cout << " l " << l << " djk " << djk << " d_Soffmindk_m " << d_Soffmindk_m << " d_mindk_xpj " << d_mindk_xpj << " son_rik " << son_rik << " soff_rik " << soff_rik << endl;
	      //cout << " allterms_x " << allterms_x << endl;
	      //exit(0);
	      d_exposure_xpj += allterms_x;
	      d_exposure_ypj += allterms_y;
	      d_exposure_zpj += allterms_z;
	    }
	  d_exposure_xpj *= (1.0/params.deltaCC2);
	  d_exposure_ypj *= (1.0/params.deltaCC2);
	  d_exposure_zpj *= (1.0/params.deltaCC2);
	  //d_exposure_xpj=0.0;
	  //d_exposure_ypj=0.0;
	  //d_exposure_zpj=0.0;
	  // d_Sonexposurei_m
	  double d_Sonexposurei_m = ds_on_dm(1.0, exposure[grid_index],\
					     params.Emin, params.deltaE);
	  double d_Sonexposurei_xpj = d_Sonexposurei_m*(1.0/params.deltaE)\
	    *d_exposure_xpj;
	  //cout << "  d_exposure_xpj" <<  d_exposure_xpj << endl;
	  double d_Sonexposurei_ypj = d_Sonexposurei_m*(1.0/params.deltaE)\
	    *d_exposure_ypj;
	  double d_Sonexposurei_zpj = d_Sonexposurei_m*(1.0/params.deltaE)\
	  *d_exposure_zpj;
	  // d_ai_xpj
	  double term1_x = s_on_exposure[grid_index] * d_Sonmindi_xpj;
	  double term2_x = s_on_mind[grid_index] * d_Sonexposurei_xpj;
	  //double term1_x = s_on_contacts[grid_index] * d_Sonmindi_xpj;
	  //double term2_x = s_on_mind[grid_index] * d_Soncontactsi_xpj;
	  //cout << " s_on_contacts " << s_on_contacts[grid_index] << " d_Sonmindi_xpj " << d_Sonmindi_xpj << endl;
	  //cout << " s_on_mind " << s_on_mind[grid_index] << " d_Soncontactsi_xpj " << d_Soncontactsi_xpj << endl;
	  //cout << "  d_Sonexposurei_xpj" <<  d_Sonexposurei_xpj << endl;
	  double d_ai_xpj = grid_s_off_bsi[grid_index] * ( term1_x + term2_x );
	  //cout << " grid_s_off_bsi[grid_index] " << grid_s_off_bsi[grid_index] << " term1_x " << term1_x << " term2_x " << term2_x << endl;
	  // d_ai_ypj
	  double term1_y = s_on_exposure[grid_index] * d_Sonmindi_ypj;
	  double term2_y = s_on_mind[grid_index] * d_Sonexposurei_ypj;
	  //double term1_y = s_on_contacts[grid_index] * d_Sonmindi_ypj;
	  //double term2_y = s_on_mind[grid_index] * d_Soncontactsi_ypj;
	  double d_ai_ypj = grid_s_off_bsi[grid_index] * ( term1_y + term2_y );
	  // d_ai_zpj
	  double term1_z = s_on_exposure[grid_index] * d_Sonmindi_zpj;
	  double term2_z = s_on_mind[grid_index] * d_Sonexposurei_zpj;
	  //double term1_z = s_on_contacts[grid_index] * d_Sonmindi_zpj;
	  //double term2_z = s_on_mind[grid_index] * d_Soncontactsi_zpj;
	  double d_ai_zpj = grid_s_off_bsi[grid_index] * ( term1_z + term2_z );
	  //cout << " i " << i << " grid_index " << grid_index << " exposure " << exposure[grid_index] << " d_exposure_xpj  " << d_exposure_xpj << endl;
	  //cout << " j atom " << j << " active grid i " << i << " activity " << activity[grid_index] << " d_ij " << dij << " d_Sonmindi_xpj " <<  d_Sonmindi_xpj << " y " << d_Sonmindi_ypj << " z " << d_Sonmindi_zpj << " d_exposure_xpj " << d_exposure_xpj << " d_Sonexposurei_xpj  " << d_Sonexposurei_xpj << endl;
	  /*cout << " j atom " << j << " i grid " << i << " activity "	\
	       << activity[grid_index] << " dij " << dij << " d_ai_xpj " \
	       << d_ai_xpj << " d_ai_ypj " << d_ai_ypj		\
	       << " d_ai_zpj " << d_ai_zpj << endl;*/
	  /*								\
	       << " d_Sonexposurei_m " << d_Sonexposurei_m		\
	       << " Exposure " << exposure[grid_index]		\
	       << " d_exposure_xpj " << d_exposure_xpj \
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
	}
      //exit(0);
      // Compute Vdrug_like derivatives
      //cout << "sum_d_ai_dj " << sum_d_ai_xpj << " " << sum_d_ai_ypj << " " << sum_d_ai_zpj << endl;
      //exit(0);
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
      //for (unsigned i=0; i < n_active_grid; i++)
      for (unsigned i=0; i < size_grid; i++)
	{
	  //int grid_index = active_grid[i];
	  int grid_index = i;
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
	  double d_Hi_xpj=0.0;
	  double d_Hi_ypj=0.0;
	  double d_Hi_zpj=0.0;
	  if (term1 >0)
	    {
	      double num_x=term1*d_apolari_xpj-		\
		apolar*(d_apolari_xpj+d_polari_xpj);
	      double num_y=term1*d_apolari_ypj-		\
		apolar*(d_apolari_ypj+d_polari_ypj);
	      double num_z=term1*d_apolari_zpj-		\
		apolar*(d_apolari_zpj+d_polari_zpj);
	      double den=term1*term1;
	      d_Hi_xpj=num_x/den;
	      d_Hi_ypj=num_y/den;
	      d_Hi_zpj=num_z/den;
	    }
	  //cout << " d_Hi_xpj " << d_Hi_xpj << " d_Hi_ypj " << d_Hi_ypj << " d_Hi_zpj " << d_Hi_zpj << endl;
	  //exit(0);
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
      double term3_y=(params.alpha*d_Va_ypj+params.beta*d_Ha_ypj);
      double term3_z=(params.alpha*d_Va_zpj+params.beta*d_Ha_zpj);

      double d_Jedi_xpj;
      double d_Jedi_ypj;
      double d_Jedi_zpj;
      if (!iszero)
	{
	  d_Jedi_xpj=Jedi*(term1_x+term2*term3_x);
	  d_Jedi_ypj=Jedi*(term1_y+term2*term3_y);
	  d_Jedi_zpj=Jedi*(term1_z+term2*term3_z);
	}
      else
	{
	  d_Jedi_xpj=0.0;
	  d_Jedi_ypj=0.0;
	  d_Jedi_zpj=0.0;
	}
      //double miniterm1_x=(1/(params.alpha*Va))*(params.alpha*d_Va_xpj);
      //double miniterm2_x=(1/Vdrug_like)*d_Vdruglike_xpj;
      //miniterm2_x=0.0;
      //double miniterm1_y=(1/(params.alpha*Va))*(params.alpha*d_Va_ypj);
      //double miniterm2_y=(1/Vdrug_like)*d_Vdruglike_ypj;
      //miniterm2_y=0.0;
      //double miniterm1_z=(1/(params.alpha*Va))*(params.alpha*d_Va_zpj);
      //double miniterm2_z=(1/Vdrug_like)*d_Vdruglike_zpj;
      //miniterm2_z=0.0;
      //double d_Jedi_xpj=Jedi*(miniterm1_x+miniterm2_x);
      //double d_Jedi_ypj=Jedi*(miniterm1_y+miniterm2_y);
      //double d_Jedi_zpj=Jedi*(miniterm1_z+miniterm2_z);

      //cout << " *** atom j *** " << j << endl;
      //cout << " d_V_xpj " << d_V_xpj << " d_V_ypj " << d_V_ypj << " d_V_zpj " << d_V_zpj <<endl;
      //cout << " d_Vdruglike_xpj " << d_Vdruglike_xpj << " d_Vdruglike_ypj " << d_Vdruglike_ypj << " d_Vdruglike_zpj " << d_Vdruglike_zpj <<endl; 
      //cout << " d_Va_xpj " << d_Va_xpj << " d_Va_ypj " << d_Va_ypj << " d_Va_zpj " << d_Va_zpj <<endl; 
      //cout << " d_Ha_xpj " << d_Ha_xpj << " d_Ha_ypj " << d_Ha_ypj << " d_Ha_zpj " << d_Ha_zpj <<endl; 

      //cout << "***atom j "  << std::fixed << std::setprecision(5) << j << " " <<  d_Jedi_xpj << " " << d_Jedi_ypj << " " << d_Jedi_zpj << " " << d_Va_xpj << " " << d_Va_ypj << " " << d_Va_zpj << " " << " " << d_Ha_xpj << " " << d_Ha_ypj << " " << d_Ha_zpj << " " << d_Vdruglike_xpj << " " << d_Vdruglike_ypj << " " << d_Vdruglike_zpj << endl;
      //double max_dJedi=100.0;
      //double min_dJedi=-100.0;
      //double max_dJedi= 0.1;
      //double min_dJedi=-0.1;
      //d_Jedi_xpj = std::min(d_Jedi_xpj,max_dJedi);
      //d_Jedi_xpj = std::max(d_Jedi_xpj,min_dJedi);
      //d_Jedi_ypj = std::min(d_Jedi_ypj,max_dJedi);
      //d_Jedi_ypj = std::max(d_Jedi_ypj,min_dJedi);
      //d_Jedi_zpj = std::min(d_Jedi_zpj,max_dJedi);
      //d_Jedi_zpj = std::max(d_Jedi_zpj,min_dJedi);

      d_Jedi_xpj_vec[j] = d_Jedi_xpj;
      d_Jedi_ypj_vec[j] = d_Jedi_ypj;
      d_Jedi_zpj_vec[j] = d_Jedi_zpj;
      //exit(0);
      sum_d_Jedi_xpj += d_Jedi_xpj;
      sum_d_Jedi_ypj += d_Jedi_ypj;
      sum_d_Jedi_zpj += d_Jedi_zpj;
      // JM also accumulate torques
      // is this right? r cross gradient
      // tx = ry gradz - rz grady
      // ty = rz gradx - rx gradz
      // tz = rx grady - ry gradx
      double d_Jedi_torque_xpj=(yj)*d_Jedi_zpj-(zj)*d_Jedi_ypj;
      double d_Jedi_torque_ypj=(zj)*d_Jedi_xpj-(xj)*d_Jedi_zpj;
      double d_Jedi_torque_zpj=(xj)*d_Jedi_ypj-(yj)*d_Jedi_xpj;
      // looks like I am doing gradient cross r below...
      //double d_Jedi_torque_xpj=d_Jedi_ypj*(zj-site_com[2])-d_Jedi_zpj*(yj-site_com[1]);
      //double d_Jedi_torque_ypj=d_Jedi_zpj*(xj-site_com[0])-d_Jedi_xpj*(zj-site_com[2]);
      //double d_Jedi_torque_zpj=d_Jedi_xpj*(yj-site_com[1])-d_Jedi_ypj*(xj-site_com[0]);
      sum_d_Jedi_torque_xpj += d_Jedi_torque_xpj;
      sum_d_Jedi_torque_ypj += d_Jedi_torque_ypj;
      sum_d_Jedi_torque_zpj += d_Jedi_torque_zpj;
    }
  //exit(0);

  // Check largest forces we have BEFORE net f/t removal
  double max_d_Jedi_der_raw[3]= {0.0,0.0,0.0};
  double max_norm_raw=0.0;
  int max_der_idx_raw=-1;
  for (int j=0; j < n_cv_atoms;j++)
    {
      double d_Jedi_dx=d_Jedi_xpj_vec[j];
      double d_Jedi_dy=d_Jedi_ypj_vec[j];
      double d_Jedi_dz=d_Jedi_zpj_vec[j];
      double norm2=d_Jedi_dx*d_Jedi_dx+d_Jedi_dy*d_Jedi_dy+d_Jedi_dz*d_Jedi_dz;
      double norm=sqrt(norm2);
      if (norm > max_norm_raw)
	{
	  max_norm_raw=norm;
	  if (j < n_apolar)
	    max_der_idx_raw = apolaratoms[j].index();
	  else
	    max_der_idx_raw = polaratoms[j-n_apolar].index();
	  max_d_Jedi_der_raw[0] = d_Jedi_dx;
	  max_d_Jedi_der_raw[1] = d_Jedi_dy;
	  max_d_Jedi_der_raw[2] = d_Jedi_dz;
	}
    }

  // Now we need to remove net force and torques introduced by JEDI
  //cout << "Jedi max derivative " << max_der_idx_raw << " " << max_d_Jedi_der_raw[0] << " " << max_d_Jedi_der_raw[1] << " " << max_d_Jedi_der_raw[2] << endl;
  //cout << "Jedi gradient " << sum_d_Jedi_xpj << " " << sum_d_Jedi_ypj << " " << sum_d_Jedi_zpj << endl;
  //cout << "Jedi torque " << sum_d_Jedi_torque_xpj << " " << sum_d_Jedi_torque_ypj << " " << sum_d_Jedi_torque_zpj << endl;

  //exit(0);

  // See SI of Ovchinnikov & Karplus J. Phys. Chem. B 2012, 116, 8584−8603
  // and Cline & Plemmons SIAM review Vol 18 January 1 1976
  // We have Ax=y and solution x is undetermined
  // We want the least-squares solution that is given by
  // z = (A+)y
  // Where A+ is the Moore-Penrose pseudo inverse matrix A+ = A*(AA*)-1
  // where A* is the transpose

  // Step 2) Form A
  unsigned nrows=6;
  unsigned ncols=3*n_cv_atoms;
  Matrix<double> A( nrows, ncols );

  for (unsigned j=0; j < ncols; j++)
    {
      if (j < n_cv_atoms)
	{
	  A[0][j] = 1.0;
	  A[1][j] = 0.0;
	  A[2][j] = 0.0;
	  A[3][j] = 0.0;
	  A[4][j] = getPosition(j)[2];//+site_com[2];//Opt by storing once all position vectors?
	  A[5][j] = -getPosition(j)[1];//-site_com[1];
	}
      else if (j < 2*n_cv_atoms)
	{
	  A[0][j] = 0.0;
	  A[1][j] = 1.0;
	  A[2][j] = 0.0;
	  A[3][j] = -getPosition(j-n_cv_atoms)[2];//-site_com[2];
	  A[4][j] = 0.0;
	  A[5][j] = getPosition(j-n_cv_atoms)[0];//+site_com[0];
	}
      else
	{
	  A[0][j] = 0.0;
	  A[1][j] = 0.0;
	  A[2][j] = 1.0;
	  A[3][j] = getPosition(j-2*n_cv_atoms)[1];//+site_com[1];
	  A[4][j] = -getPosition(j-2*n_cv_atoms)[0];//-site_com[0];
	  A[5][j] = 0.0;
	}
    }

  /*cout << "*** A *** " << endl;
  for (int i=0; i < nrows; i++)
    {
      for (int j=0; j < ncols ; j++)
	{
	  //A[i][j] = 1.0;
	  cout << A[i][j] << " ";
	}
      cout << endl;
      }*/
  // Step 3) Compute Moore-Penrose pseudo inverse A+ = A*(AA*)-1
  Matrix<double> Aplus( ncols, nrows);
  // Plumed implementation uses singular value decomposition via plumed_lapack_dgesdd
  // Is this the most numerically stable option?
  pseudoInvert(A, Aplus);

  /*cout << "*** A+ *** " << endl;
  for (int i=0; i < ncols; i++)
    {
      for (int j=0; j < nrows ; j++)
	{
	  cout << Aplus[i][j] << " ";
	}
      cout << endl;
      }*/

  // Multiple iterations
  // Is this wise? adds more random forces,
  // may eventually blow up integrator?
  int max_it = 1;
  double tol= 0.01;
  for (int k=0; k < max_it;k++)
    {
      //cout << "cancelling net force/torque step " << k << " ... " << endl;
      // Step 1) Form y
      vector<double> y;
      y.resize(6);
      y[0] = -sum_d_Jedi_xpj;
      y[1] = -sum_d_Jedi_ypj;
      y[2] = -sum_d_Jedi_zpj;
      y[3] = -sum_d_Jedi_torque_xpj;
      y[4] = -sum_d_Jedi_torque_ypj;
      y[5] = -sum_d_Jedi_torque_zpj;

      /*cout << "*** y" << endl;
	for (int i=0; i < 6; i++)
	{
	cout << y[i] << endl;
	}
      */

      // Step 4) Compute z = (A+)y
      vector<double> z;
      z.resize(3*n_cv_atoms);
      //cout << "y size " << y.size() << endl;
      //cout << "Aplus.ncols " << Aplus.ncols() << endl;
      //cout << "Aplus.nrows " << Aplus.nrows() << endl;
      //cout << " z size " << z.size() << endl;

      mult(Aplus,y,z);
      //cout << "*** z" << endl;
      //for (int i=0;i<3*n_cv_atoms;i++)
      //	{
      // cout << "i " << i << " " << z[i] << endl;
      //	}

      //cout << " *** sum me " << endl;
      double sum_z[3] = {0.0,0.0,0.0};
      double sum_rcrossz[3] = {0.0,0.0,0.0};
      for (int j=0; j < n_cv_atoms; j++)
	{
	  //int j = i/3;
	  double xj = getPosition(j)[0];
	  double yj = getPosition(j)[1];
	  double zj = getPosition(j)[2];
	  double fx = z[j+0*n_cv_atoms];
	  double fy = z[j+1*n_cv_atoms];
	  double fz = z[j+2*n_cv_atoms];
	  //cout << " j " << j << " fx " << fx << " fy " << fy << " fz " << fz << endl;
	  sum_z[0] += fx;
	  sum_z[1] += fy;
	  sum_z[2] += fz;
	  sum_rcrossz[0] += yj*fz-zj*fy;
	  sum_rcrossz[1] += zj*fx-xj*fz;
	  sum_rcrossz[2] += xj*fy-yj*fx;
	}
	  //z[i] = 0.0; z[i+1] =0.0; z[i+2] = 0.0;
      //cout << "sum_z " << sum_z[0] << " " << sum_z[1] << " " << sum_z[2] << endl;
      //cout << "sum_rcrossz " << sum_rcrossz[0] << " " << sum_rcrossz[1] << " " << sum_rcrossz[2] << endl;
      //exit(0);
      // Step 5) Now correct derivatives and verify that net forces and torque are
      // close to zero
      double sum_d_Jedistar_xpj=0.0;
      double sum_d_Jedistar_ypj=0.0;
      double sum_d_Jedistar_zpj=0.0;
      double sum_d_Jedistar_torque_xpj=0.0;
      double sum_d_Jedistar_torque_ypj=0.0;
      double sum_d_Jedistar_torque_zpj=0.0;
      for ( unsigned j=0; j < n_cv_atoms ; j++)
	{
	  double xj = getPosition(j)[0];
	  double yj = getPosition(j)[1];
	  double zj = getPosition(j)[2];
	  double d_Jedistar_xpj = d_Jedi_xpj_vec[j] + z[j+0*n_cv_atoms];
	  double d_Jedistar_ypj = d_Jedi_ypj_vec[j] + z[j+1*n_cv_atoms];
	  double d_Jedistar_zpj = d_Jedi_zpj_vec[j] + z[j+2*n_cv_atoms];
	  double d_Jedistar_torque_xpj = (yj)*d_Jedistar_zpj-(zj)*d_Jedistar_ypj;
	  double d_Jedistar_torque_ypj = (zj)*d_Jedistar_xpj-(xj)*d_Jedistar_zpj;
	  double d_Jedistar_torque_zpj = (xj)*d_Jedistar_ypj-(yj)*d_Jedistar_xpj;
	  sum_d_Jedistar_xpj += d_Jedistar_xpj;
	  sum_d_Jedistar_ypj += d_Jedistar_ypj;
	  sum_d_Jedistar_zpj += d_Jedistar_zpj;
	  sum_d_Jedistar_torque_xpj += d_Jedistar_torque_xpj;
	  sum_d_Jedistar_torque_ypj += d_Jedistar_torque_ypj;
	  sum_d_Jedistar_torque_zpj += d_Jedistar_torque_zpj;

	  d_Jedi_xpj_vec[j] = d_Jedistar_xpj;
	  d_Jedi_ypj_vec[j] = d_Jedistar_ypj;
	  d_Jedi_zpj_vec[j] = d_Jedistar_zpj;
	}

      //cout << "sum_d_Jedistar_der " << sum_d_Jedistar_xpj << " " << sum_d_Jedistar_ypj << " " << sum_d_Jedistar_zpj << endl; 
      double l2_norm2 = sum_d_Jedistar_xpj*sum_d_Jedistar_xpj + sum_d_Jedistar_ypj*sum_d_Jedistar_ypj + sum_d_Jedistar_zpj*sum_d_Jedistar_zpj;
      double l2_norm = sqrt(l2_norm2);
      //cout << " norm d_Jedisar " << l2_norm << endl;
      //cout << "sum_d_Jedistar_torque_der " << sum_d_Jedistar_torque_xpj << " " << sum_d_Jedistar_torque_ypj << " " << sum_d_Jedistar_torque_zpj << endl;

      sum_d_Jedi_xpj= sum_d_Jedistar_xpj;
      sum_d_Jedi_ypj= sum_d_Jedistar_ypj;
      sum_d_Jedi_zpj= sum_d_Jedistar_zpj;
      sum_d_Jedi_torque_xpj= sum_d_Jedistar_torque_xpj;
      sum_d_Jedi_torque_ypj= sum_d_Jedistar_torque_ypj;;
      sum_d_Jedi_torque_zpj= sum_d_Jedistar_torque_zpj;

      if (l2_norm < tol)
	break;

    }

  //cout << "Final Jedi gradient " << sum_d_Jedi_xpj << " " << sum_d_Jedi_ypj << " " << sum_d_Jedi_zpj << endl;
  //cout << "Final Jedi torque " << sum_d_Jedi_torque_xpj << " " << sum_d_Jedi_torque_ypj << " " << sum_d_Jedi_torque_zpj << endl;

  double max_d_Jedi_der[3]= {0.0,0.0,0.0};
  double max_norm=0.0;
  int max_der_idx=-1;
  for (int j=0; j < n_cv_atoms;j++)
    {
      double d_Jedi_dx=d_Jedi_xpj_vec[j];
      double d_Jedi_dy=d_Jedi_ypj_vec[j];
      double d_Jedi_dz=d_Jedi_zpj_vec[j];
      double norm2=d_Jedi_dx*d_Jedi_dx+d_Jedi_dy*d_Jedi_dy+d_Jedi_dz*d_Jedi_dz;
      double norm=sqrt(norm2);
      if (norm > max_norm)
	{
	  max_norm=norm;
	  if (j < n_apolar)
	    max_der_idx = apolaratoms[j].index();
	  else
	    max_der_idx = polaratoms[j-n_apolar].index();
	  max_d_Jedi_der[0] = d_Jedi_dx;
	  max_d_Jedi_der[1] = d_Jedi_dy;
	  max_d_Jedi_der[2] = d_Jedi_dz;
	}
      setAtomsDerivatives(j,Vector(d_Jedi_dx,d_Jedi_dy,d_Jedi_dz));
      //
      unsigned pdb_idx = 0;
      if (j < n_apolar)
	pdb_idx = apolaratoms[j].index();
      else
	pdb_idx = polaratoms[j-n_apolar].index();
      //cout << "***atom j "  << std::fixed << std::setprecision(5) << j << " ( " << pdb_idx << " ) " << d_Jedi_dx << " " << d_Jedi_dy << " " << d_Jedi_dz << endl;
      //
    }

  //exit(0);

  mod = step % stride;
  //cout << " we are at step " << step << " stride is " << stride << " mod is " << mod << endl;
  //<< " " << sum_d_Jedi_xpj << " " << sum_d_Jedi_ypj << " " << sum_d_Jedi_zpj << " " << sum_d_Jedi_torque_xpj << " " << sum_d_Jedi_torque_ypj << " " << sum_d_Jedi_torque_zpj 

  iszero = mod;

  if (!iszero)
    {
      ofstream wfile;
      wfile.open(summary_file.c_str(),std::ios_base::app);
      wfile << std::setprecision(5) << step << " " << Jedi << " " << Vdrug_like << " " << Va << " " << Ha 
	    << " " << max_der_idx_raw << " " << max_d_Jedi_der_raw[0] << " " << max_d_Jedi_der_raw[1] << " " << max_d_Jedi_der_raw[2]  
	    << " " << max_der_idx << " " << max_d_Jedi_der[0] << " " << max_d_Jedi_der[1] << " " << max_d_Jedi_der[2] 
	    << " " << rmsd << endl;
      wfile.close();
    }

  // Occasionally save grids
  mod = step % gridstride;
  iszero = mod;
  //cout << " gridstride is " << gridstride << endl;
  if (!iszero)
    {
      // Here for write oordinates of updated grid to a XYZ file
      // This can be used to check that the grid has been trans/roted correctly
      ofstream wfile;
      string gridfilename = "grid-step-";
      string tail;
      string s;
      stringstream out;
      out << step;
      s = out.str();
      tail.append(s);
      string tail2;
      tail2=tail;
      string gridfilename2;
      gridfilename2 = gridfilename;
      tail.append(".dx");
      gridfilename.append(tail);
      wfile.open(gridfilename.c_str());
      // XYZ
      //wfile << size_grid << endl;
      //wfile << "comment" << endl;
      //for (int i=0; i < size_grid; i++)
      //	{
      //	  wfile << "C " << std::fixed << std::setprecision(5) << new_x[i]*10 << " " << new_y[i]*10 << " " << new_z[i]*10 << endl;
      //}
      // DX
      wfile << "object 1 class gridpositions counts " << grid_extent[0] << " " << grid_extent[1] << " " << grid_extent[2] << endl;

      double origin[3] = {0.0,0.0,0.0};
      double step_z[3] = {0.0,0.0,0.0};
      double delta_z[3] = {0.0,0.0,0.0};
      double step_y[3] = {0.0,0.0,0.0};
      double delta_y[3] = {0.0,0.0,0.0};
      double step_x[3] = {0.0,0.0,0.0};
      double delta_x[3] = {0.0,0.0,0.0};
      origin[0] = new_x[grid_origin_idx];
      origin[1] = new_y[grid_origin_idx];
      origin[2] = new_z[grid_origin_idx];
      step_z[0] = new_x[1];
      step_z[1] = new_y[1];
      step_z[2] = new_z[1];
      delta_z[0] = step_z[0] - origin[0];
      delta_z[1] = step_z[1] - origin[1];
      delta_z[2] = step_z[2] - origin[2];
      step_y[0] = new_x[ grid_extent[2] ];
      step_y[1] = new_y[ grid_extent[2] ];
      step_y[2] = new_z[ grid_extent[2] ];
      delta_y[0] = step_y[0] - origin[0];
      delta_y[1] = step_y[1] - origin[1];
      delta_y[2] = step_y[2] - origin[2];
      step_x[0] = new_x[grid_extent[1]*grid_extent[2] ];
      step_x[1] = new_y[grid_extent[1]*grid_extent[2]];
      step_x[2] = new_z[grid_extent[1]*grid_extent[2]];
      delta_x[0] = step_x[0] - origin[0];
      delta_x[1] = step_x[1] - origin[1];
      delta_x[2] = step_x[2] - origin[2];


      //cout << " origin " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
      //cout << " step_z " << step_z[0] << " " << step_z[1] << " " << step_z[2] << endl;
      //cout << " delta_z " << delta_z[0] << " " << delta_z[1] << " " << delta_z[2] << endl;      
      //cout << " step_y " << step_y[0] << " " << step_y[1] << " " << step_y[2] << endl;
      //cout << " delta_y " << delta_y[0] << " " << delta_y[1] << " " << delta_y[2] << endl; 
      //cout << " step_x " << step_x[0] << " " << step_x[1] << " " << step_x[2] << endl;
      //cout << " delta_x " << delta_x[0] << " " << delta_x[1] << " " << delta_x[2] << endl; 

      //exit(0);

      //wfile << "origin " << grid_min[0]*10.0 << " " << grid_min[1]*10.0 << " " << grid_min[2]*10.0 << endl;
      wfile << "origin " << origin[0]*10.0 << " " << origin[1]*10.0 << " " << origin[2]*10.0 << endl;
      // UPDATE CODE HERE
      // FIND COORDS FIRST GRID POINT ALONG Z --> 1
      // FIND COORDS FIRST GRID POINT ALONG Y --> n_gz +1
      // FIND COORDS FIRST GRID POINT ALONG X --> n_gz*g_y + n_gz + 1
      // CALC DELTA X, DELTA Y, DELTA Z
      wfile << "delta " << delta_x[0]*10.0 << " " << delta_x[1]*10.0 << " " << delta_x[2]*10.0 << endl;
      wfile << "delta " << delta_y[0]*10.0 << " " << delta_y[1]*10.0 << " " << delta_y[2]*10.0 << endl;
      wfile << "delta " << delta_z[0]*10.0 << " " << delta_z[1]*10.0 << " " << delta_z[2]*10.0 << endl;
      wfile << "object 2 class gridconnections counts " << grid_extent[0] << " " << grid_extent[1] << " " << grid_extent[2] << endl;
      wfile << "object 3 class array type double rank 0 items " << size_grid << endl;
      for (int i=0; i < grid_extent[0]; i++)
	{
	  for (int j=0; j < grid_extent[1] ; j++)
	    {
	      for (int k=0 ; k < grid_extent[2] ; k++)
		{
		  // (0,0,0) --> 0
		  // (0,0,n_gz) --> n_gz
		  // (0,1,0) --> n_gz + 1
		  // (0,1,n_gz) --> n_gz + n_gz
		  // (a,b,c) --> a*n_gz*n_gy + b* n_gz + k
		  // (1,1,1) --> 1*n_gz*n_gy + 1*n_gz + 1
		  int index=i*grid_extent[2]*grid_extent[1] + j*grid_extent[2] + k;
		  double ai = activity[index];
		  wfile <<  std::fixed << std::setprecision(5) << ai;
		  wfile << " ";
		}
	      wfile << endl;
	    }
	}
      wfile.close();
      tail2.append(".xyz");
      gridfilename2.append(tail2);
      //cout << gridfilename2 << endl;
      //exit(0);
      wfile.open(gridfilename2.c_str());
      // XYZ
      wfile << size_grid << endl;
      wfile << "comment" << endl;
      for (unsigned i=0; i < size_grid; i++)
      	{
      	  wfile << "C " << std::fixed << std::setprecision(5) << new_x[i]*10 << " " << new_y[i]*10 << " " << new_z[i]*10 << endl;
      }
      wfile.close();
      // output detailed grid statistics in
      // dx files with updated position and 'activitiy'
      //                                and 'Ha'
      //exit(0);
    }

  //exit(0);

}//close jedi::calculate 
}//close namespace colvar
}//cloase namespace PLMD
