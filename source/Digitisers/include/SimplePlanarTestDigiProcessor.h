#ifndef SimplePlanarTestDigiProcessor_h
#define SimplePlanarTestDigiProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>

#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;


/** ======= SimplePlanarTestDigiProcessor ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The SimTrackerHits should come from a planar detector like VXD, SIT, SET or FTD.
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in u and v direction. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param SimTrkHitRelCollection The name of the TrackerHit SimTrackerHit relation collection <br>
 * (default name VTXTrackerHitRelations) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * <br>
 * 
 */
class SimplePlanarTestDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new SimplePlanarTestDigiProcessor ; }
  
  
  SimplePlanarTestDigiProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
protected:
  
  std::string _inColName ;
  
  std::string _outColName ;
  std::string _outRelColName ;
 
    
  int _nRun ;
  int _nEvt ;
  
  float _resU ;
  float _resV ;
  

  gsl_rng* _rng ;
  
  
} ;

#endif



