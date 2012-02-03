#include "FTDSpacePointBuilder.h"

#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include "IMPL/LCFlagImpl.h"
#include "UTIL/ILDConf.h"

#include "marlin/VerbosityLevels.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"
#include "marlin/Global.h"
#include <UTIL/ILDConf.h>

#include "MarlinTrk/util/MeasurementSurfaceStore.h"
#include "MarlinTrk/util/MeasurementSurface.h"
#include "MarlinTrk/util/ICoordinateSystem.h"
#include "MarlinTrk/util/CartesianCoordinateSystem.h"

#include <cmath>

using namespace lcio ;
using namespace marlin ;



FTDSpacePointBuilder aFTDSpacePointBuilder ;


FTDSpacePointBuilder::FTDSpacePointBuilder() : Processor("FTDSpacePointBuilder") {

   // modify processor description
   _description = "FTDSpacePointBuilder...." ;


   // register steering parameters: name, description, class-variable, default value
   registerInputCollection(LCIO::TRACKERHIT,
                           "FTDTrackerHitCollection",
                           "FTDTrackerHitCollection",
                           _FTDTrackerHitCollection,
                           std::string("FTDTrackerHits")); 


   registerOutputCollection(LCIO::TRACKERHIT,
                            "SpacePointsCollection",
                            "SpacePointsCollection",
                            _SpacePointsCollection,
                            std::string("FTDSpacePoints"));


  
   
}




void FTDSpacePointBuilder::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

 
  MarlinTrk::GearExtensions::MeasurementSurfaceStore::Instance().initialise(Global::GEAR);

}


void FTDSpacePointBuilder::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void FTDSpacePointBuilder::processEvent( LCEvent * evt ) { 


  LCCollection* col = evt->getCollection( _FTDTrackerHitCollection ) ;


  if( col != NULL ){
    
    
    unsigned createdSpacePointsFromPixel = 0;
    unsigned createdSpacePointsFromStrip = 0;
    unsigned rawStripHits = 0;
    unsigned possibleSpacePointsFromStrips = 0;
    _nOutOfBoundary = 0;
    _nStripsTooParallel = 0;
    
    LCCollectionVec * spCol = new LCCollectionVec(LCIO::TRACKERHIT);
    
    /*LCFlagImpl hitFlag(0) ;
    hitFlag.setBit( LCIO::TRBIT_HITS ) ;
    trkCol->setFlag( hitFlag.getFlag()  ) ;*/ //TODO: do I need to set any flags of the collections?
    
    
    unsigned nHits = col->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << "Number of hits on the FTDs: " << nHits <<"\n";
    
    //store hits in map according to their CellID0
    std::map< int , std::vector< TrackerHitPlane* > > map_cellID0_hits;
    std::map< int , std::vector< TrackerHitPlane* > >::iterator it;
    
    for( unsigned i=0; i<nHits; i++){
      
      TrackerHitPlane* trkHit = dynamic_cast<TrackerHitPlane*>( col->getElementAt( i ) );
      if( trkHit != NULL) map_cellID0_hits[ trkHit->getCellID0() ].push_back( trkHit ); 
      
    }
    
    
    
    
    const gear::FTDParameters& ftdParams = Global::GEAR->getFTDParameters() ;
    const gear::FTDLayerLayout& ftdLayers = ftdParams.getFTDLayerLayout() ;
    
    
    // now loop over all CellID0s
    for( it= map_cellID0_hits.begin(); it!= map_cellID0_hits.end(); it++ ){
     
      
      std::vector< TrackerHitPlane* > hitsFront = it->second;
      
      //find out layer, module, sensor
      
      UTIL::BitField64  cellID( ILDCellID0::encoder_string );
      cellID.setValue( it->first );
      
      int subdet = cellID[ ILDCellID0::subdet ] ;
      int side   = cellID[ ILDCellID0::side ];
      int module = cellID[ ILDCellID0::module ];
      int sensor = cellID[ ILDCellID0::sensor ];
      int layer  = cellID[ ILDCellID0::layer ];
      int cellID0 = cellID.lowWord();
      
      
      //check if this is a double sided petal
      if( ftdLayers.isDoubleSided( layer ) ){ 
        
        rawStripHits += it->second.size();
        
        //check if this sensor is in the front or the back:
        if( sensor <= ftdLayers.getNSensors( layer ) / 2 ){ // in front
          
          
          
          //get the hits on the back:
          UTIL::BitField64  cellIDBack( ILDCellID0::encoder_string );
          cellIDBack.setValue( cellID0 );
          
          int sensorBack = sensor + ftdLayers.getNSensors( layer ) / 2;
          cellIDBack[ ILDCellID0::sensor ] = sensorBack;
          int cellID0Back = cellIDBack.lowWord();
          
          std::vector< TrackerHitPlane* > hitsBack = map_cellID0_hits[ cellID0Back ];
          
          streamlog_out( DEBUG3 ) 
            << "strips: CellID0 " << cellID0  << " (su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensor
            << ")(" << hitsFront.size()
            << " hits) <---> CellID0 " << cellID0Back << " (su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensorBack
            << ")(" << hitsBack.size() << " hits)\n";
          streamlog_out( DEBUG3 )
            << "--> " << hitsFront.size() * hitsBack.size() << " possible combinations\n";
            
          possibleSpacePointsFromStrips += hitsFront.size() * hitsBack.size();
            
          // Now iterate over all combinations and store those that make sense
          for( unsigned i=0; i<hitsFront.size(); i++ ){
            
            TrackerHitPlane* hitFront = hitsFront[i];
            
            for( unsigned j=0; j<hitsBack.size(); j++ ){
              
              
              TrackerHitPlane* hitBack = hitsBack[j];
              
              
              TrackerHitImpl* spacePoint = createSpacePoint( hitFront, hitBack );
              if ( spacePoint == NULL ) continue;
              
              
              CellIDEncoder<TrackerHitImpl> cellid_encoder( ILDCellID0::encoder_string , spCol );
              cellid_encoder[ ILDCellID0::subdet ] = ILDDetID::FTD  ;
              cellid_encoder[ ILDCellID0::side   ] = side;
              cellid_encoder[ ILDCellID0::layer  ] = layer;
              cellid_encoder[ ILDCellID0::module ] = module;
              cellid_encoder[ ILDCellID0::sensor ] = sensor; // use the sensors from the front
              
              cellid_encoder.setCellID( spacePoint ) ;
              
              // store the hits it's composed of:
              spacePoint->rawHits().push_back( hitFront );
              spacePoint->rawHits().push_back( hitBack );
              
              spacePoint->setType( 0 ); // TODO: so what type do we really want to set?
              
              spCol->addElement( spacePoint ) ; 
              
              createdSpacePointsFromStrip++;
              
              
            }
            
          }
          
        }
        
        
      }
      else { // is not double sided
        
        
        for( unsigned i=0; i<hitsFront.size(); i++ ){
          
          TrackerHit* trkHitPlane = hitsFront[i];
          // Just make a TrackerHitImpl version of the hit
          TrackerHitImpl* trackerHit = new TrackerHitImpl();
          
          trackerHit->setPosition( trkHitPlane->getPosition() );
          trackerHit->rawHits().push_back( trkHitPlane );
          
          trackerHit->setCellID0( trkHitPlane->getCellID0() ) ;
          trackerHit->setType( 0 ); // TODO: so what type do we really want to set?
          
          spCol->addElement( trackerHit ) ; 
          
          createdSpacePointsFromPixel++;
          
        }
        
        
      }
      
      
      
      
    }
    
    evt->addCollection(spCol,_SpacePointsCollection.c_str());
    
    streamlog_out( DEBUG4 )<< "\nCreated " << createdSpacePointsFromPixel + createdSpacePointsFromStrip
      << " space points ( " << createdSpacePointsFromPixel << " from Pixel, and " 
      << createdSpacePointsFromStrip << " from overlapping strips (raw strip hits: " << rawStripHits << ") )\n";
    
    streamlog_out( DEBUG3 ) << "  There were " << rawStripHits << " strip hits available, giving " 
      << possibleSpacePointsFromStrips << " possible space points.\n";
    
    streamlog_out( DEBUG3 ) << "  " << _nStripsTooParallel << " space points couldn't be created, cause the strips were too parallel\n";
    streamlog_out( DEBUG3 ) << "  " << _nOutOfBoundary << " space points couldn't be created, cause the result was outside the sensor boundary\n"; 
      
      
    streamlog_out( DEBUG4 ) << "\n";
    
  }


  _nEvt ++ ;
  
}





void FTDSpacePointBuilder::check( LCEvent * evt ) {}


void FTDSpacePointBuilder::end(){
   
   
}

TrackerHitImpl* FTDSpacePointBuilder::createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b ){
  
  const double* p1 = a->getPosition();
  double x1 = p1[0];
  double y1 = p1[1];
  double z1 = p1[2];
  const float* v1 = a->getV();
  float ex1 = cos( v1[1] ) * sin( v1[0] ); 
  float ey1 = sin( v1[1] ) * sin( v1[0] );
  
  const double* p2 = b->getPosition();
  double x2 = p2[0];
  double y2 = p2[1];
  double z2 = p2[2];
  const float* v2 = b->getV();
  float ex2 = cos( v2[1] ) * sin( v2[0] ); 
  float ey2 = sin( v2[1] ) * sin( v2[0] );
  
  streamlog_out( DEBUG2 ) << "\t ( " << x1 << " " << y1 << " " << z1 << " ) <--> ( " << x2 << " " << y2 << " " << z2 << " )\n";
  
  double x=0.;
  double y=0.;
  
  if ( calculateXingPoint( x1, y1, ex1, ey1, x2, y2, ex2, ey2, x, y ) != 0 ){
    
    _nStripsTooParallel++;
    streamlog_out( DEBUG2 ) << "\tStrips too parallel\n\n";
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
  
  }
  
  double z= (z1 + z2)/2.;
  
  streamlog_out( DEBUG2 ) << "\tPosition of space point (global) : ( " << x << " " << y << " " << z << " )\n";
  
  // Check if the new hit is within the boundary
  CLHEP::Hep3Vector globalPoint(x,y,z);
  MarlinTrk::GearExtensions::MeasurementSurface* ms = MarlinTrk::GearExtensions::MeasurementSurfaceStore::Instance().GetMeasurementSurface( a->getCellID0() );
  CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
  localPoint.setZ( 0. ); // we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  if( !ms->isLocalInBoundary( localPoint ) ){
    
    _nOutOfBoundary++;
    streamlog_out( DEBUG2 ) << "\tHit is out of boundary: local coordinates are ( " 
      << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )\n\n";
    
    return NULL;
    
  }
  
  
  //Create the new TrackerHit
  TrackerHitImpl* spacePoint = new TrackerHitImpl();
  
  double pos[3] = {x,y,z};
  spacePoint->setPosition(  pos  ) ;
  
  streamlog_out( DEBUG2 ) << "\tHit accepted\n\n";

  return spacePoint;
  
}

int FTDSpacePointBuilder::calculateXingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ){


  float a = (x1*ey1 - y1*ex1) - (x2*ey1 - y2*ex1);
  float b = ex2*ey1 - ex1*ey2;

  const float epsilon = 0.00001;

  if( fabs(b) < epsilon ) return 1; // if b==0 the two directions e1 and e2 are parallel and there is no crossing!

  float t = a/b;

  x = x2 + t*ex2;
  y = y2 + t*ey2;

  return 0;

  

}
 
