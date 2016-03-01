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
#ifndef __PLUMED_vesselbase_ActionWithVessel_h
#define __PLUMED_vesselbase_ActionWithVessel_h

#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "tools/Exception.h"
#include "tools/DynamicList.h"
#include "tools/MultiValue.h"
#include <vector>

namespace PLMD{
class Value;
class Stopwatch;

namespace vesselbase{

class Vessel;
class BridgeVessel;
class StoreDataVessel;

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that are computed by calculating the same function multiple
times.  This is used in PLMD::MultiColvar.
*/

class ActionWithVessel : public virtual Action {
friend class Vessel;
friend class ShortcutVessel;
friend class FunctionVessel;
friend class StoreDataVessel;
friend class BridgeVessel;
friend class ActionWithInputVessel;
private:
/// Do all calculations in serial
  bool serial;
/// Lower memory requirements
  bool lowmem;
/// Are we skipping the calculation of the derivatives
  bool noderiv;
/// This tells plumed that this is used in a bridge
  bool actionIsBridged;
/// The maximum number of derivatives we can use before we need to invoke lowmem
  unsigned maxderivatives;
/// The tolerance on the accumulators 
  double tolerance;
/// Tolerance for quantities being put in neighbor lists
  double nl_tolerance;
/// Pointers to the functions we are using on each value
  std::vector<Vessel*> functions;
/// A pointer to the object that stores data
  StoreDataVessel* mydata;
/// Tempory storage for forces
  std::vector<double> tmpforces;
/// Ths full list of tasks we have to perform
  std::vector<unsigned> fullTaskList;
/// The current number of active tasks
  unsigned nactive_tasks;
/// The indices of the tasks in the full list of tasks
  std::vector<unsigned> indexOfTaskInFullList;
/// The list of currently active tasks
  std::vector<unsigned> partialTaskList;
/// This list is used to update the neighbor list
  std::vector<unsigned> taskFlags;
/// The list of atoms involved in derivatives (we keep a copy here to avoid resizing)
  std::vector<unsigned> der_list;
/// The buffer that we use (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
/// Do we want to output information on the timings of different parts of the calculation
  bool timers;
/// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch;
protected:
/// The terms in the series are locked
  bool contributorsAreUnlocked;
/// Does the weight have derivatives
  bool weightHasDerivatives;
/// This is used for numerical derivatives of bridge variables
  unsigned bridgeVariable;
/// Add a vessel to the list of vessels
  void addVessel( const std::string& name, const std::string& input, const int numlab=0 );
  void addVessel( Vessel* vv );
/// Add a bridging vessel to the list of vessels
  BridgeVessel* addBridgingVessel( ActionWithVessel* tome );
/// Complete the setup of this object (this routine must be called after construction of ActionWithValue)
  void readVesselKeywords();
/// Turn on the derivatives in the vessel
  void needsDerivatives();
/// Return the value of the tolerance
  double getTolerance() const ;
/// Return the value for the neighbor list tolerance
  double getNLTolerance() const ;
/// Get the number of vessels
  unsigned getNumberOfVessels() const;
/// Get a pointer to the ith vessel
   Vessel* getPntrToVessel( const unsigned& i );
/// Calculate the values of all the vessels
  void runAllTasks();
/// Resize all the functions when the number of derivatives change
  void resizeFunctions();
/// This loops over all the vessels calculating them and also 
/// sets all the element derivatives equal to zero
  bool calculateAllVessels( const unsigned& taskCode, MultiValue& myvals, MultiValue& bvals, std::vector<double>& buffer, std::vector<unsigned>& der_list );
/// Retrieve the forces from all the vessels (used in apply)
  bool getForcesFromVessels( std::vector<double>& forcesToApply );
/// Is the calculation being done in serial
  bool serialCalculation() const;
/// Are we using low memory
  bool usingLowMem() const ;
/// Set that we are using low memory
  void setLowMemOption(const bool& );
/// Get the number of tasks that are currently active
  unsigned getCurrentNumberOfActiveTasks() const ;
/// Get the ith of the currently active tasks
  unsigned getActiveTask( const unsigned& ii ) const ;
/// Deactivate all the tasks in the task list
  void deactivateAllTasks();
/// Deactivate all tasks with i in lower \f$\le\f$  i < upper
  void deactivateTasksInRange( const unsigned& lower, const unsigned& upper );
/// Get the size of the buffer
  unsigned getSizeOfBuffer( unsigned& bufsize );
/// Add a task to the full list
  void addTaskToList( const unsigned& taskCode );
public:
  static void registerKeywords(Keywords& keys);
  explicit ActionWithVessel(const ActionOptions&ao);
  ~ActionWithVessel();
  void unlockContributors();
  void lockContributors();
  virtual void finishTaskListUpdate(){};
/// Activate the jth colvar
/// Deactivate the current task in future loops
  virtual void deactivate_task( const unsigned & task_index );
/// Are derivatives required for this quantity
  bool derivativesAreRequired() const ;
/// Finish running all the calculations
  virtual void finishComputations( const std::vector<double>& buffer );
/// Are the base quantities periodic
  virtual bool isPeriodic()=0;
/// What are the domains of the base quantities
  virtual void retrieveDomain( std::string& min, std::string& max);
/// Get the number of derivatives for final calculated quantity 
  virtual unsigned getNumberOfDerivatives()=0;
/// Get the number of quantities that are calculated during each task
  virtual unsigned getNumberOfQuantities();
/// Get the list of indices that have derivatives
//  virtual void getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices );
/// Switch on additional tasks 
  void activateTheseTasks( std::vector<unsigned>& addtionalTasks );
/// Do any jobs that are required before the task list is undertaken
  virtual void doJobsRequiredBeforeTaskList();
/// Get the full size of the taskList dynamic list
  unsigned getFullNumberOfTasks() const ;
/// Get the position of the ith active task in the full list
  unsigned getPositionInFullTaskList( const unsigned& ii ) const ;
/// Get the code for the ii th task in the list
  unsigned getTaskCode( const unsigned& ii ) const ;
/// Calculate one of the functions in the distribution
  virtual void performTask( const unsigned& , const unsigned& , MultiValue& ) const=0;
/// Do the task if we have a bridge
  virtual void transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const;
/// Ensure that data required in other vessels is stored
  StoreDataVessel* buildDataStashes( const bool& allow_wcutoff, const double& wtol );
/// Apply forces from bridge vessel - this is rarely used - currently only in ActionVolume
  virtual void applyBridgeForces( const std::vector<double>& bb ){ plumed_error(); }
/// These are overwritten in MultiColvarFunction
//  virtual void activateIndexes( const unsigned&, const unsigned&, const std::vector<unsigned>& ){}
/// Return a particular named vessel
  Vessel* getVesselWithName( const std::string& mynam );
/// Does the weight have derivatives
  bool weightWithDerivatives() const ;
};

inline
double ActionWithVessel::getTolerance() const {
  return tolerance;
}

inline
double ActionWithVessel::getNLTolerance() const {
  return nl_tolerance;
}

inline
unsigned ActionWithVessel::getNumberOfVessels() const {
  return functions.size();
}

inline
unsigned ActionWithVessel::getNumberOfQuantities(){
  return 2;
}

inline
Vessel* ActionWithVessel::getPntrToVessel( const unsigned& i ){
  plumed_dbg_assert( i<functions.size() );
  return functions[i];
}

inline
unsigned ActionWithVessel::getFullNumberOfTasks() const {
  return fullTaskList.size();
}

inline
unsigned ActionWithVessel::getTaskCode( const unsigned& ii ) const {
  plumed_dbg_assert( ii<fullTaskList.size() );
  return fullTaskList[ii];
}

inline
unsigned ActionWithVessel::getCurrentNumberOfActiveTasks() const {
  return nactive_tasks;
}

inline
unsigned ActionWithVessel::getActiveTask( const unsigned& ii ) const {
  plumed_dbg_assert( ii<nactive_tasks );
  return partialTaskList[ii];
}

inline
unsigned ActionWithVessel::getPositionInFullTaskList( const unsigned& ii ) const {
  plumed_dbg_assert( ii<nactive_tasks );
  return indexOfTaskInFullList[ii];
}

inline
bool ActionWithVessel::serialCalculation() const {
  return serial;
}

inline
bool ActionWithVessel::usingLowMem() const {
  return lowmem;
}

inline
void ActionWithVessel::setLowMemOption(const bool& l){
  lowmem=l;
}

inline
bool ActionWithVessel::derivativesAreRequired() const {
  return !noderiv;
}

inline
bool ActionWithVessel::weightWithDerivatives() const {
  return weightHasDerivatives;
}

} 
}
#endif
