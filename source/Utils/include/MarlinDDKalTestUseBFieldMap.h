#ifndef MarlinDDKalTestUseBFieldMap_h
#define MarlinDDKalTestUseBFieldMap_h 1

#include "marlin/Processor.h"
#include <string>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include "TObject.h"     // from ROOT
#include "TEveTrackPropagator.h"     // from ROOT



/** MarlinDDKalTestUseBFieldMap:
 *  Initialize MarlinDDKalTest to use the (non-homogeneous) B-field map from DD4hep,
 *  
 *  Call before any tracking processor in your steering file - after DD4hep is initialized.
 * 
 * @author F. Gaede, DESY
 * @date June 2018 
 * @version $Id:$
 */

class MarlinDDKalTestUseBFieldMap : public marlin::Processor {
  
protected:
  /// internal helper class for the B-field
  class DD4hepMagField : public TEveMagField {
  public:
    DD4hepMagField(){ 
      fFieldConstant = kFALSE;
    }
    virtual ~DD4hepMagField(){
    }
    
    using   TEveMagField::GetField;
    virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const { 
      double posV[3] = { x*dd4hep::mm, y*dd4hep::mm, z*dd4hep::mm }  ;
      double bfieldV[3] ;
      _det->field().magneticField( posV  , bfieldV  ) ;
      
      //      streamlog_out( DEBUG ) 
//      std::cout << " GetFieldD( " << x*dd4hep::mm <<", "  << y*dd4hep::mm << ", " <<  z*dd4hep::mm  << " ) -> "
//		<<  bfieldV[0]/dd4hep::tesla << ", " << bfieldV[1]/dd4hep::tesla << ", " << bfieldV[2]/dd4hep::tesla << std::endl ;

      
      return TEveVectorD( bfieldV[0]/dd4hep::tesla , bfieldV[1]/dd4hep::tesla, bfieldV[2]/dd4hep::tesla ) ;
    }
    
    dd4hep::Detector* _det = &dd4hep::Detector::getInstance()  ;
  };

 
public:
  
  virtual Processor*  newProcessor() { return new MarlinDDKalTestUseBFieldMap ; }
  
  
  MarlinDDKalTestUseBFieldMap() ;
  
  /** Initiallize the DD4hep geometry.
   */
  virtual void init() ;
  
  /// do nothing
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /// do nothing
  virtual void processEvent( LCEvent * evt ) ; 
  
  /// do nothing
  virtual void check( LCEvent * evt ) ; 
  
  /// do nothing
  virtual void end() ;
  
protected:

  DD4hepMagField* _field=nullptr ;

} ;

#endif



