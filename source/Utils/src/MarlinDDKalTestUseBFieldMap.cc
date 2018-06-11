#include "MarlinDDKalTestUseBFieldMap.h"

#include "TBField.h"

using namespace marlin ;

MarlinDDKalTestUseBFieldMap aMarlinDDKalTestUseBFieldMap ;


MarlinDDKalTestUseBFieldMap::MarlinDDKalTestUseBFieldMap() : Processor("MarlinDDKalTestUseBFieldMap") {

  _description = "MarlinDDKalTestUseBFieldMap configure MarlinDDKalTest to use the DD4hep B-fiel map" ;
}

void MarlinDDKalTestUseBFieldMap::init() { 
  
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  
  streamlog_out( MESSAGE )  << " ---- use B-field map for detector " << theDetector.header().name()  << std::endl ;

  _field =  new DD4hepMagField ;

  TBField::SetBfieldPtr(_field );

  TBField::SetUseUniformBfield(kFALSE);

}


void MarlinDDKalTestUseBFieldMap::processRunHeader( LCRunHeader*) { } 

void MarlinDDKalTestUseBFieldMap::processEvent( LCEvent * ) { }

void MarlinDDKalTestUseBFieldMap::check( LCEvent *  ) { }

void MarlinDDKalTestUseBFieldMap::end(){ 
  delete _field ;
}


// root stuff
//ClassImp(TEveMagField)
