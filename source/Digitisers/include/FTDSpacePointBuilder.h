#ifndef FTDSpacePointBuilder_h
#define FTDSpacePointBuilder_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitImpl.h>

using namespace lcio ;
using namespace marlin ;





/** ================= FTD Space Point Builder =================
 * 
 * Builds space points for the FTD. The silicon strip detectors of the FTD consist of sensors on the front
 * and on the back of the petals. 
 * The digitisers create TrackerHitPlanars for the front and the back strips. 
 * In order to get a spacepoint (as is needed by track reconstruction) those strip measurements need to be combined into one space point.
 * 
 * This is done by this processor. 
 * 
 *  <h4>Input - Prerequisites</h4>
 *  
 * The TrackerHitPlanars as created by the Digitisers for the FTD. 
 * This could of course be generalized for other detectors as well, but is specific at the moment
 * as it uses the gear information of the FTD.
 *
 *  <h4>Output</h4> 
 *  
 * A collection of TrackerHits containing all relevant spacepoint hits from the FTD:
 * In there will be the spacepoints combined from the strip measurements and 
 * as well the pixel hits unchanged (as they are already valid spacepoints).
 * The TrackerHits that are created from two strip measurements will store the 
 * strip hits in their rawHits. Also the spacepoints created from the Pixels will have
 * the pixel hits stored in their rawHits.
 * 
 * @param FTDTrackerHitCollection The name of the input collection of TrackerHits on the FTD <br>
 * (default name FTDTrackerHits) <br>
 * 
 * @param SpacePointsCollection The name of the output collection of the created spacepoints on the FTD <br>
 * (default name FTDSpacePoints) <br>
 * 
 * @author Robin Glattauer HEPHY, Vienna
 *
 */
class FTDSpacePointBuilder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FTDSpacePointBuilder ; }
  
  
  FTDSpacePointBuilder() ;
  
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



  /** Input collection name.
  */
  std::string _FTDTrackerHitCollection;

  /** Output collection name.
  */
  std::string _SpacePointsCollection;

  /** Calculates the 2 dimensional crossing point of two lines.
   * Each line is specified by a point (x,y) and a direction vector (ex,ey).
   * 
   * @return 0, if the calculation has been successful, another number if there was a problem.
   */
  static int calculateXingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y );
 
  /** @return a spacepoint (in the form of a TrackerHitImpl* ) created from two TrackerHitPlane* which stand for si-strips */
  TrackerHitImpl* createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b );
 

  int _nRun ;
  int _nEvt ;

  unsigned _nOutOfBoundary;
  unsigned _nStripsTooParallel;


} ;

#endif



