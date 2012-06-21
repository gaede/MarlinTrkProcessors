#include "FullLDCTracking_MarlinTrk.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include "ClusterShapes.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include "gear/FTDLayerLayout.h"
#include "gear/FTDParameters.h"
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include <climits>
#include <cmath>

using namespace lcio ;
using namespace marlin ;
using namespace MarlinTrk ;

FullLDCTracking_MarlinTrk aFullLDCTracking_MarlinTrk ;

FullLDCTracking_MarlinTrk::FullLDCTracking_MarlinTrk() : Processor("FullLDCTracking_MarlinTrk") {  
  _description = "Performs full tracking in ILD detector" ;  
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // Input tracker hit collections
  
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "FTDPixelHitCollectionName",
                          "FTD Pixel Hit Collection Name",
                          _FTDPixelHitCollection,
                          std::string("FTDPixelTrackerHits"));  
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "FTDSpacePointCollectionName",
                          "FTD FTDSpacePoint Collection Name",
                          _FTDSpacePointCollection,
                          std::string("FTDSpacePoints"));  
  
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "VTXHitCollection",
                          "VTX Hit Collection Name",
                          _VTXTrackerHitCollection,
                          std::string("VTXTrackerHits"));  
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "SITHitCollection",
                          "SIT Hit Collection Name",
                          _SITTrackerHitCollection,
                          std::string("SITTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "SETHitCollection",
                          "SET Hit Collection Name",
                          _SETTrackerHitCollection,
                          std::string("SETTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "ETDHitCollection",
                          "ETD Hit Collection Name",
                          _ETDTrackerHitCollection,
                          std::string("ETDTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "TPCHitCollection",
                          "TPC Hit Collection Name",
                          _TPCTrackerHitCollection,
                          std::string("TPCTrackerHits"));
  
  
  // Input track collections
  
  registerInputCollection(LCIO::TRACK,
                          "TPCTracks",
                          "TPC Track Collection",
                          _TPCTrackCollection,
                          std::string("TPCTracks"));
  
  registerInputCollection(LCIO::TRACK,
                          "SiTracks",
                          "Si Track Collection",
                          _SiTrackCollection,
                          std::string("SiTracks"));
  
  // Input relation collections
  
  registerInputCollection(LCIO::LCRELATION,
                          "TPCTracksMCPRelColl",
                          "TPC Track to MCP Relation Collection Name",
                          _TPCTrackMCPCollName,
                          std::string("TPCTracksMCP"));
  
  registerInputCollection(LCIO::LCRELATION,
                          "SiTracksMCPRelColl",
                          "Si Track to Collection",
                          _SiTrackMCPCollName,
                          std::string("SiTracksMCP"));
  
  // Output track collection
  registerOutputCollection(LCIO::TRACK,
                           "LDCTrackCollection",
                           "LDC track collection name",
                           _LDCTrackCollection,
                           std::string("LDCTracks"));
  
  
  
  // steering parameters
  
  registerProcessorParameter("D0CutForMerging",
                             "Cut on D0 difference for merging of Si and TPC segments",
                             _d0CutForMerging,
                             float(500.0));
  
  registerProcessorParameter("Z0CutForMerging",
                             "Cut on Z0 difference for merging of Si and TPC segments",
                             _z0CutForMerging,
                             float(1000.0));
  
  registerProcessorParameter("OmegaCutForMerging",
                             "Cut on Omega difference for merging Si and TPC segments",
                             _dOmegaForMerging,
                             float(0.25));
  
  registerProcessorParameter("AngleCutForMerging",
                             "Cut on Opening Angle for merging Si and TPC segments",
                             _angleForMerging,
                             float(0.10));
  
  registerProcessorParameter("Chi2FitCut",
                             "Cut on fit Chi2",
                             _chi2FitCut,
                             float(100.0));

  registerProcessorParameter("Debug",
                             "Activate debugging?",
                             _debug,
                             int(0));
  
  
  registerProcessorParameter("ForceSiTPCMerging",
                             "Force merging of Si and TPC segments?",
                             _forceMerging,
                             int(0));
  
  registerProcessorParameter("D0CutForForcedMerging",
                             "Cut on D0 difference for forced merging of Si and TPC segments",
                             _d0CutForForcedMerging,
                             float(50.));
  
  registerProcessorParameter("Z0CutForForcedMerging",
                             "Cut on Z0 difference for forced merging of Si and TPC segments",
                             _z0CutForForcedMerging,
                             float(200.));
  
  registerProcessorParameter("OmegaCutForForcedMerging",
                             "Cut on Omega difference for forced merging of Si and TPC segments",
                             _dOmegaForForcedMerging,
                             float(0.15));
  
  registerProcessorParameter("AngleCutForForcedMerging",
                             "Cut on Opening Angle for forced merging of Si and TPC segments",
                             _angleForForcedMerging,
                             float(0.05));
  
  registerProcessorParameter("ForceTPCSegmentsMerging",
                             "Force merging of TPC Segments?",
                             _mergeTPCSegments,
                             int(1));
  
  registerProcessorParameter("D0CutToMergeTPCSegments",
                             "Cut on D0 difference for merging TPC segments",
                             _d0CutToMergeTPC,
                             float(100.));
  
  registerProcessorParameter("Z0CutToMergeTPCSegments",
                             "Cut on Z0 difference for merging TPC segments",
                             _z0CutToMergeTPC,
                             float(5000.0));
  
  registerProcessorParameter("DeltaPCutToMergeTPCSegments",
                             "Cut on dP/P difference for merging TPC segments",
                             _dPCutToMergeTPC,
                             float(0.1));
  
  registerProcessorParameter("PtCutToMergeTPCSegments",
                             "Cut on Pt of tracks for merging TPC segments",
                             _PtCutToMergeTPC,
                             float(1.2));
  
  
  
  registerProcessorParameter("cosThetaCutHighPtMerge",
                             "Cut on cos theta between the two momentum vectors when considering merger of high Pt tracks",
                             _cosThetaCutHighPtMerge,
                             float(0.99));
  
  registerProcessorParameter("cosThetaCutSoftHighPtMerge",
                             "cut on cos theta between the two momentum vectors when considering merger of high Pt tracks for softer dp/p cut",
                             _cosThetaCutSoftHighPtMerge,
                             float(0.998));
  
  registerProcessorParameter("momDiffCutHighPtMerge",
                             "cut on dp/p when considering merger of high Pt tracks",
                             _momDiffCutHighPtMerge,
                             float(0.01));
  
  registerProcessorParameter("momDiffCutSoftHighPtMerge",
                             "softer cut on dp/p when considering merger of high Pt tracks",
                             _momDiffCutSoftHighPtMerge,
                             float(0.25));
  
  registerProcessorParameter("hitDistanceCutHighPtMerge",
                             "cut on 3D distance between hit and helix extrapolation when considering merger of high Pt tracks",
                             _hitDistanceCutHighPtMerge,
                             float(25.0));
  
  registerProcessorParameter("maxHitDistanceCutHighPtMerge",
                             "cut for max 3D distance between any hit and helix extrapolation when considering merger of high Pt tracks",
                             _maxHitDistanceCutHighPtMerge,
                             float(50.0));
  
  registerProcessorParameter("maxFractionOfOutliersCutHighPtMerge",
                             "cut on maximum fraction of outliers when considering merger of high Pt tracks",
                             _maxFractionOfOutliersCutHighPtMerge,
                             float(0.95));
  
  
  
  
  
  registerProcessorParameter("CutOnTPCHits",
                             "Cut on the number of the TPC hits for tracks with no Si hits",
                             _cutOnTPCHits,
                             int(35));
  registerProcessorParameter("CutOnSiHits",
                             "Cut on the number of the Si hits for tracks with no TPC hits",
                             _cutOnSiHits,
                             int(4));
  
  
  registerProcessorParameter("AssignVTXHits",
                             "Assign left over VTX hits",
                             _assignVTXHits,
                             int(1));
  
  registerProcessorParameter("AssignFTDHits",
                             "Assign left over FTD hits",
                             _assignFTDHits,
                             int(1));
  
  registerProcessorParameter("AssignSITHits",
                             "Assign left over SIT hits",
                             _assignSITHits,
                             int(1));
  
  registerProcessorParameter("AssignTPCHits",
                             "Assign left over TPC hits",
                             _assignTPCHits,
                             int(1));
  
  registerProcessorParameter("AssignSETHits",
                             "Assign SET Hits",
                             _assignSETHits,
                             int(1));
  
  
  registerProcessorParameter("AssignETDHits",
                             "Assign ETD Hits",
                             _assignETDHits,
                             int(1));
  
  registerProcessorParameter("NHitsExtrapolation",
                             "number of hits for outward extrapolation",
                             _nHitsExtrapolation,
                             int(35));
    
  
  registerProcessorParameter("VTXHitToTrackDistance",
                             "Cut on distance between track and VTX hits",
                             _distCutForVTXHits,
                             float(1.5));
  
  
  registerProcessorParameter("FTDHitToTrackDistance",
                             "Cut on distance between track and FTD hits",
                             _distCutForFTDHits,
                             float(2.0));
  
  
  registerProcessorParameter("SITHitToTrackDistance",
                             "Cut on distance between track and SIT hits",
                             _distCutForSITHits,
                             float(2.0));
  
  registerProcessorParameter("SETHitToTrackDistance",
                             "Cut on distance between track and SET hits",
                             _distCutForSETHits,
                             float(2.0));
  
  
  registerProcessorParameter("ETDHitToTrackDistance",
                             "Cut on distance between track and ETD hits",
                             _distCutForETDHits,
                             float(10.0));
  
  
  registerProcessorParameter("TPCHitToTrackDistance",
                             "Cut on distance between track and TPC hits",
                             _distCutForTPCHits,
                             float(15.0));
  
  
  registerProcessorParameter("CutOnTrackD0",
                             "Cut on the track parameter D0",
                             _d0TrkCut,
                             float(500.));
  
  
  registerProcessorParameter("CutOnTrackZ0",
                             "Cut on the track parameter Z0",
                             _z0TrkCut,
                             float(500.));
  
  
  registerProcessorParameter("ForbidOverlapInZTPC",
                             "Forbid overlap in Z for the merged TPC segments",
                             _forbidOverlapInZTPC,
                             int(0));
  
  registerProcessorParameter("ForbidOverlapInZComb",
                             "Forbid overlap in Z for combining TPC segments with tracks having Si hits",
                             _forbidOverlapInZComb,
                             int(0));
  
  registerProcessorParameter("MultipleScatteringOn",
                             "Use MultipleScattering in Fit",
                             _MSOn,
                             bool(true));
  
  registerProcessorParameter("EnergyLossOn",
                             "Use Energy Loss in Fit",
                             _ElossOn,
                             bool(true));
  
  registerProcessorParameter("SmoothOn",
                             "Smooth All Mesurement Sites in Fit",
                             _SmoothOn,
                             bool(true));
  
  
  registerProcessorParameter( "InitialTrackErrorD0",
                             "Value used for the initial d0 variance of the trackfit",
                             _initialTrackError_d0,
                             float(1.e6));
  
  registerProcessorParameter( "InitialTrackErrorPhi0",
                             "Value used for the initial phi0 variance of the trackfit",
                             _initialTrackError_phi0,
                             float(1.e2));
  
  registerProcessorParameter( "InitialTrackErrorOmega",
                             "Value used for the initial omega variance of the trackfit",
                             _initialTrackError_omega,
                             float(1.e-4));
  
  registerProcessorParameter( "InitialTrackErrorZ0",
                             "Value used for the initial z0 variance of the trackfit",
                             _initialTrackError_z0,
                             float(1.e6));
  
  registerProcessorParameter( "InitialTrackErrorTanL",
                             "Value used for the initial tanL variance of the trackfit",
                             _initialTrackError_tanL,
                             float(1.e2));
  
  
  registerProcessorParameter( "MaxChi2PerHit",
                             "Maximum Chi-squared value allowed when assigning a hit to a track",
                             _maxChi2PerHit,
                             double(1.e2));
  
  
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  registerOptionalParameter("RunMarlinTrkDiagnostics", "Run MarlinTrk Diagnostics. MarlinTrk must be compiled with MARLINTRK_DIAGNOSTICS_ON defined", _runMarlinTrkDiagnostics, bool(false));
  
  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics", _MarlinTrkDiagnosticsName, std::string("FullLDCTrackingDiagnostics"));    
  
#endif
  
}



void FullLDCTracking_MarlinTrk::init() { 
  
  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;
  
  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  if( _trksystem == 0 ){
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  void * dcv = _trksystem->getDiagnositicsPointer();
  DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
  dc->init(_MarlinTrkDiagnosticsName,_MarlinTrkDiagnosticsName, _runMarlinTrkDiagnostics);
  
#endif
  
  this->setupGearGeom(Global::GEAR);
  
}

void FullLDCTracking_MarlinTrk::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
  _nEvt = 0;
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "FullLDCTracking_MarlinTrk ---> new run : run number = " << _nRun << std::endl;
  
} 

void FullLDCTracking_MarlinTrk::processEvent( LCEvent * evt ) { 
  
  _evt = evt;
  
  
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "FullLDCTracking_MarlinTrk -> run = " << _nRun 
  << "  event = " << _nEvt << std::endl;
  streamlog_out(DEBUG4) << std::endl;
  
  
  prepareVectors( evt );
  streamlog_out(DEBUG3) << "prepareVectors done..." << std::endl;
  streamlog_out(DEBUG3) << "************************************Merge TPC/Si" << std::endl;
  MergeTPCandSiTracks();
  streamlog_out(DEBUG3) << "************************************Merging done..." << std::endl;
  MergeTPCandSiTracksII();
  streamlog_out(DEBUG3) << "************************************Merging II done..." << std::endl;
  Sorting(_allCombinedTracks);
  streamlog_out(DEBUG3) << "************************************Sorting done..." << std::endl;
  SelectCombinedTracks();
  streamlog_out(DEBUG3) << "************************************Selection of combined tracks done..." << std::endl;
  AddNotCombinedTracks( );
  streamlog_out(DEBUG3) << "************************************Not combined tracks added..." << std::endl;
  //CheckTracks( );
  
  AddNotAssignedHits();
  streamlog_out(DEBUG3) << "***********************************Not assigned hits added..." << std::endl;
  AddTrackColToEvt(evt,_trkImplVec,
                   _LDCTrackCollection);
  streamlog_out(DEBUG3) << "Collections added to event..." << std::endl;
  CleanUp();
  streamlog_out(DEBUG3) << "Cleanup is done..." << std::endl;
  _nEvt++;
  //  getchar();
  streamlog_out(DEBUG3) << std::endl;
  streamlog_out(DEBUG3) << std::endl;
  
}

void FullLDCTracking_MarlinTrk::AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
                                                 std::string TrkColName) {
  
  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  colTRK->setFlag( trkFlag.getFlag()  ) ;  
  
  streamlog_out(DEBUG4)<< "Collection " << TrkColName << " is being added to event " << std::endl;
  
  //  LCCollectionVec * colRel = NULL;
  
  
  int nTrkCand = int(trkVec.size());
  
  int nTotTracks = 0;
  float eTot = 0.0;
  float pxTot = 0.0;
  float pyTot = 0.0;
  float pzTot = 0.0;
  
  //SJA:FIXME: So here we are going to do one final refit. This can certainly be optimised, but rather than worry about the mememory management right now lets make it work, and optimise it later ...
  
  
  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    
    TrackExtended * trkCand = trkVec[iTRK];
    TrackerHitExtendedVec& hitVec = trkCand->getTrackerHitExtendedVec();
    
    EVENT::TrackerHitVec trkHits;
    
    int nHits = int(hitVec.size());
    for (int ihit=0;ihit<nHits;++ihit) {
      
      EVENT::TrackerHit* trkHit = hitVec[ihit]->getTrackerHit();
      if(trkHit) { 
        trkHits.push_back(trkHit);   
      }
      else{
        throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: TrackerHit pointer == NULL ")  ) ;
      }
      
    }
    
    
    if( trkHits.size() < 3 ) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Cannot fit less than 3 hits. Number of hits =  " << trkHits.size() << std::endl;
      continue ; 
    }
    
    
    TrackImpl* Track = new TrackImpl ;
    
    // setup initial dummy covariance matrix
    EVENT::FloatVec covMatrix;
    covMatrix.resize(15);
    
    for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }
    
    covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
    
    
    // sort hits in R
    std::vector< std::pair<float, EVENT::TrackerHit*> > r2_values;
    r2_values.reserve(trkHits.size());
    
    for (TrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
      EVENT::TrackerHit* h = *it;
      float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
      r2_values.push_back(std::make_pair(r2, *it));
    }
    
    sort(r2_values.begin(),r2_values.end());
    
    trkHits.clear();
    trkHits.reserve(r2_values.size());
    
    for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
      trkHits.push_back(it->second);
    }
    
    bool fit_backwards = IMarlinTrack::backward;
    
    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
    
    
    int error = 0;
    
    try {
      
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, Track, fit_backwards, covMatrix, _bField, _maxChi2PerHit);                              
      
    } catch (...) {
      
      //      delete Track;
      //      delete marlinTrk;
      
      throw ;
      
    }
    
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    if ( error != IMarlinTrack::success && _runMarlinTrkDiagnostics ) {        
      void * dcv = _trksystem->getDiagnositicsPointer();
      DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
      dc->skip_current_track();
    }        
#endif
    
    
    std::vector<std::pair<EVENT::TrackerHit* , double> > hits_in_fit ;  
    std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
    std::vector<TrackerHit*> all_hits;    
    all_hits.reserve(300);
    
    marlinTrk->getHitsInFit(hits_in_fit);
    
    for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }
    
    UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ; 
    
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);
    
    marlinTrk->getOutliers(outliers);
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }
    
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);
    
    
    delete marlinTrk;
    
    if( error != IMarlinTrack::success ) {       
      
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit failed with error code " << error << " track dropped. Number of hits = "<< trkHits.size() << std::endl;  
      
      delete Track;      
      continue ;
    }
    
    if( Track->getNdf() < 0) {       
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit returns " << Track->getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << std::endl;       
      
      delete Track;
      continue ;
    }
    
    const TrackState* trkStateIP = Track->getTrackState(lcio::TrackState::AtIP);
    
    if (trkStateIP == 0) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit returns " << Track->getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << std::endl;       
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: trkStateIP pointer == NULL ")  ) ;
    }
    
    
    
    GroupTracks * group = trkCand->getGroupTracks();
    
    if (group != NULL) {
      TrackExtendedVec trkVecGrp = group->getTrackExtendedVec();
      int nGrTRK = int(trkVecGrp.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
        TrackExtended * subTrack = trkVecGrp[iGr];
        Track->addTrack(subTrack->getTrack());
      }
    }
    
    float d0TrkCand = trkCand->getD0();
    float z0TrkCand = trkCand->getZ0();
    //    float phi0TrkCand = trkCand->getPhi();
    
    
    int hits_in_vxd = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ];
    int hits_in_ftd = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ];
    int hits_in_sit = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ];
    int hits_in_tpc = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ];
    
    int nHitsSi = hits_in_vxd + hits_in_ftd + hits_in_sit;
    
    bool rejectTrack = (hits_in_tpc < _cutOnTPCHits) && (nHitsSi<=0);
    
    rejectTrack = rejectTrack || ( (hits_in_tpc<=0) && (nHitsSi<_cutOnSiHits) );
    rejectTrack = rejectTrack || ( fabs(d0TrkCand) > _d0TrkCut ) || ( fabs(z0TrkCand) > _z0TrkCut );
    
    if ( rejectTrack ) {
      
      delete Track;
      
    } else {
      
      float omega = trkStateIP->getOmega();
      float tanLambda = trkStateIP->getTanLambda();
      float phi0 = trkStateIP->getPhi();
      float d0 = trkStateIP->getD0();
      float z0 = trkStateIP->getZ0();
      
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      
      float trkPx = helix.getMomentum()[0];
      float trkPy = helix.getMomentum()[1];
      float trkPz = helix.getMomentum()[2];
      float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
      
      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;   
      nTotTracks++;
      
      colTRK->addElement(Track);
      
    }
  }
  
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "Number of accepted " << TrkColName << " = " 
  << nTotTracks << std::endl;
  streamlog_out(DEBUG4) << "Total 4-momentum of " << TrkColName << " : E = " << eTot
  << " Px = " << pxTot
  << " Py = " << pyTot
  << " Pz = " << pzTot << std::endl;
  streamlog_out(DEBUG4) << std::endl;
  
  evt->addCollection(colTRK,TrkColName.c_str());
  
  
}


void FullLDCTracking_MarlinTrk::prepareVectors(LCEvent * event ) {
  
  
  
  
  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allSETHits.clear();
  _allETDHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _allNonCombinedTPCTracks.clear();
  _allNonCombinedSiTracks.clear();
  _trkImplVec.clear();
  _candidateCombinedTracks.clear();
  
  
  std::map <TrackerHit*,TrackerHitExtended*> mapTrackerHits;
  
  // Reading TPC hits
  try {
    
    LCCollection * col = event->getCollection(_TPCTrackerHitCollection.c_str());
    
    int nelem = col->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of TPC hits = " << nelem << std::endl;
    
    for (int ielem=0;ielem<nelem;++ielem) {
      
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
      // For no error in r
      
      double tpcRPhiRes = sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);
      double tpcZRes = sqrt(hit->getCovMatrix()[5]);
      
      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _TPCTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  };
  
  
  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    
    LCCollection * hitCollection = event->getCollection(_FTDPixelHitCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of FTD Pixel hits = " << nelem << std::endl;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      
      double point_res_rphi = sqrt( hit->getdU()*hit->getdU() + hit->getdV()*hit->getdV() );
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));            
      hitExt->setDet(int(INT_MAX));
      
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));
      
      if ( _petalBasedFTDWithOverlaps == true ) {
        
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals 
        if ( petalIndex % 2 == 0 ) {
          layer = 2*layer;
        }
        else {
          layer = 2*layer + 1;
        }
        
      }
      
      if (layer >= _nLayersFTD) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nLayersFTD <<  std::endl;
        exit(1);
      }
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit->getPosition()[i];
      }      
      
      streamlog_out( DEBUG1 ) << " FTD Pixel Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2]  << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << std::endl;  
      
    }
  }
  catch(DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _FTDPixelHitCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  // Reading in FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    
    LCCollection * hitCollection = event->getCollection(_FTDSpacePointCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of FTD SpacePoints hits = " << nelem << std::endl;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      
      TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      
      // SJA:FIXME: fudge for now by a factor of two and ignore covariance
      double point_res_rphi = 2 * sqrt( hit->getCovMatrix()[0] + hit->getCovMatrix()[2] );
      
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
      
      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));
      
      if ( _petalBasedFTDWithOverlaps == true ) {
        
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals 
        if ( petalIndex % 2 == 0 ) {
          layer = 2*layer;
        }
        else {
          layer = 2*layer + 1;
        }
        
      }
      
      if (layer >= _nLayersFTD) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nLayersFTD <<  std::endl;
        exit(1);
      }
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit->getPosition()[i];
      }      
      
      streamlog_out( DEBUG1 ) << " FTD SpacePoint Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2]  << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << std::endl;  
      
      
    }
  }
  catch(DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _FTDSpacePointCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  //  // Reading ETD Hits
  //  try {
  //    LCCollection * col = event->getCollection(_ETDTrackerHitCollection.c_str());
  //    int nelem = col->getNumberOfElements();
  //    for (int ielem=0;ielem<nelem;++ielem) {
  //      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
  //      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
  //      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[0])));
  //      hitExt->setResolutionZ(0.1);
  //      //      hitExt->setResolutionRPhi(_resolutionRPhi_FTD);
  //      //      hitExt->setResolutionZ(_resolutionZ_FTD);
  //      // type and det are no longer used, set to INT_MAX to try and catch any missuse
  //      hitExt->setType(int(INT_MAX));
  //      hitExt->setDet(int(INT_MAX));
  //      _allETDHits.push_back( hitExt );
  //      mapTrackerHits[hit] = hitExt;
  //    }
  //  }
  //  catch( DataNotAvailableException &e ) {
  //      streamlog_out(DEBUG4) << _ETDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  //  }
  
  
  try {
    
    LCCollection *hitCollection = event->getCollection(_SITTrackerHitCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of SIT hits = " << nelem << std::endl;
    
    TrackerHit*          trkhit   = 0;
    TrackerHitPlane*     trkhit_P = 0;
    TrackerHitZCylinder* trkhit_C = 0;
    
    double drphi(NAN);
    double dz(NAN);
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      
      // hit could be of the following type
      // 1) TrackerHit, either ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT or just standard TrackerHit
      // 2) TrackerHitPlane, either 1D or 2D
      // 3) TrackerHitZCylinder, if coming from a simple cylinder design as in the LOI
      
      // Establish which of these it is in the following order of likelyhood
      //    i)   ILDTrkHitTypeBit::ONE_DIMENSIONAL (TrackerHitPlane) Should Never Happen: SpacePoints Must be Used Instead
      //    ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit)
      //    iii) TrackerHitPlane (Two dimentional)
      //    iv)  TrackerHitZCylinder 
      //    v)   Must be standard TrackerHit
      
      trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
      
      int layer = getLayerID(trkhit);
      
      if (layer < 0 || (unsigned)layer >= _nLayersSIT) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in SIT : layer is outside allowed range : " << layer << std::endl;
        exit(1);
      }
      
      // first check that we have not been given 1D hits by mistake, as they won't work here
      if ( BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
        
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
        
      } 
      // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
      else if ( BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ) {
        
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);         
        dz    =      sqrt(trkhit->getCovMatrix()[5]);         
        
      } 
      // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
      else if ( ( trkhit_P = dynamic_cast<TrackerHitPlane*>( hitCollection->getElementAt( ielem ) ) ) )  {
        
        // first we need to check if the measurement vectors are aligned with the global coordinates 
        gear::Vector3D U(1.0,trkhit_P->getU()[1],trkhit_P->getU()[0],gear::Vector3D::spherical);
        gear::Vector3D V(1.0,trkhit_P->getV()[1],trkhit_P->getV()[0],gear::Vector3D::spherical);
        gear::Vector3D Z(0.0,0.0,1.0);
        
        const float eps = 1.0e-07;
        // V must be the global z axis 
        if( fabs(1.0 - V.dot(Z)) > eps ) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SIT Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }
        
        // U must be normal to the global z axis
        if( fabs(U.dot(Z)) > eps ) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SIT Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }
        
        drphi = trkhit_P->getdU();
        dz    = trkhit_P->getdV();                                                 
        
      } 
      // or a simple cylindrical design, as used in the LOI      
      else if ( ( trkhit_C = dynamic_cast<TrackerHitZCylinder*>( hitCollection->getElementAt( ielem ) ) ) ) {
        
        drphi = trkhit_C->getdRPhi();
        dz    = trkhit_C->getdZ();
        
      } 
      // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
      else {
        
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);         
        dz =     sqrt(trkhit->getCovMatrix()[5]);             
        
      }
      
      // now that the hit type has been established carry on and create a 
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( trkhit );
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);
      
      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      _allSITHits.push_back( hitExt );
      mapTrackerHits[trkhit] = hitExt;
      
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = trkhit->getPosition()[i];
      }
      
      streamlog_out( DEBUG1 ) << " SIT Hit " <<  trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << std::endl;  
      
    }
    
  } catch(DataNotAvailableException &e) {
    streamlog_out( DEBUG4 ) << " collection not found : " << _SITTrackerHitCollection.c_str() << std::endl ;
  }
  
  
  
  // Reading SET Hits
  //    sprintf(strg,"SET hits ----->\n");
  //    streamlog_out(DEBUG4) << strg;
  //    sprintf(strg," id    r_hit   phi_hit    z_hit  e(r-phi)  e(z)  \n");
  //    streamlog_out(DEBUG4) << strg;
  
  //     "  0   1807.64   -1.97     26.45   0.010   0.010
  //     "  1   1812.64   -1.97     26.51   0.010   0.010
  //  try {
  //    LCCollection * col = event->getCollection(_SETTrackerHitCollection.c_str());
  //    int nelem = col->getNumberOfElements();
  //    for (int ielem=0;ielem<nelem;++ielem) {
  //      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
  //      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
  //      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
  //      hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
  //      //       float x_hit = float(hit->getPosition()[0]);
  //      //       float y_hit = float(hit->getPosition()[1]);
  //      //       float z_hit = float(hit->getPosition()[2]);
  //      //       float phi_hit = atan2( y_hit, x_hit );
  //      //       float r_hit = sqrt( x_hit*x_hit + y_hit*y_hit );
  //      //       sprintf(strg,"%3i  %8.2f  %6.3f  %8.2f  %6.3f  %6.3f   ", 
  //      //             ielem, r_hit, phi_hit, z_hit, 
  //      //             hitExt->getResolutionRPhi(), 
  //      //             hitExt->getResolutionZ() );
  //      //       streamlog_out(DEBUG4) << strg;
  //
  //      //       std::cout << hit << std::endl;
  //      // type and det are no longer used, set to INT_MAX to try and catch any missuse
  //      hitExt->setType(int(INT_MAX));
  //      hitExt->setDet(int(INT_MAX));
  //      _allSETHits.push_back( hitExt );
  //      mapTrackerHits[hit] = hitExt;
  //    }
  //  }
  //  catch( DataNotAvailableException &e ) {
  //      streamlog_out(DEBUG4) << _SETTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  //  }
  
  // Reading VTX Hits
  try {
    LCCollection * col = event->getCollection(_VTXTrackerHitCollection.c_str());
    
    int nelem = col->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of VXD hits = " << nelem << std::endl;
    
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHitPlane * trkhit = dynamic_cast<TrackerHitPlane*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(trkhit);
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(trkhit->getdU());
      hitExt->setResolutionZ(trkhit->getdV());
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));      
      hitExt->setDet(int(INT_MAX));
      _allVTXHits.push_back( hitExt );
      mapTrackerHits[trkhit] = hitExt;
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = trkhit->getPosition()[i];
      }
      
      int layer = getLayerID(trkhit);
      
      streamlog_out( DEBUG1 ) << " VXD Hit " <<  trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << std::endl;  
      
      
    }
    
    
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  // Reading TPC Tracks
  try {
    LCCollection * col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG3) << std::endl;
    streamlog_out(DEBUG3) << "Number of TPC Tracks = " << nelem << std::endl;
    streamlog_out(DEBUG3) << " Trk       p          D0         Z0       Px       Py       Pz    ntpc ndf Chi2/ndf" << std::endl;
    //           "  0  1.111   0.059      0.022    -0.54     0.61    -0.45    0.185
    
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      
      Track * tpcTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      
      TrackerHitVec hitVec = tpcTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      
      float d0TPC = tpcTrack->getD0();
      float z0TPC = tpcTrack->getZ0();
      float omegaTPC = tpcTrack->getOmega();
      float phi0TPC = tpcTrack->getPhi();
      float tanLTPC = tpcTrack->getTanLambda();
      float Chi2TPC = tpcTrack->getChi2()/float(tpcTrack->getNdf());
      const int ndfTPC = tpcTrack->getNdf();
      
      HelixClass helixTPC;
      
      helixTPC.Initialize_Canonical(phi0TPC,d0TPC,z0TPC,omegaTPC,tanLTPC,_bField);
      
      char strg[200];
      
      float pxTPC = helixTPC.getMomentum()[0];
      float pyTPC = helixTPC.getMomentum()[1];
      float pzTPC = helixTPC.getMomentum()[2];
      const float ptot = sqrt(pxTPC*pxTPC+pyTPC*pyTPC+pzTPC*pzTPC);
      sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
              ptot, d0TPC,z0TPC,pxTPC,pyTPC,pzTPC,nHits,ndfTPC,Chi2TPC);
      streamlog_out(DEBUG3) << strg << std::endl;
      
      TrackExtended * trackExt = new TrackExtended( tpcTrack );
      
      trackExt->setOmega(tpcTrack->getOmega());
      trackExt->setTanLambda(tpcTrack->getTanLambda());
      trackExt->setPhi(tpcTrack->getPhi());
      trackExt->setD0(tpcTrack->getD0());
      trackExt->setZ0(tpcTrack->getZ0());
      float cov[15];
      float param[5];
      //      float reso[4];
      param[0] = tpcTrack->getOmega();
      param[1] = tpcTrack->getTanLambda();
      param[2] = tpcTrack->getPhi();
      param[3] = tpcTrack->getD0();
      param[4] = tpcTrack->getZ0();
      
      
      const FloatVec Cov = tpcTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      }
      
      
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack->getNdf());
      trackExt->setChi2(tpcTrack->getChi2());            
      
      
      
      for (int iHit=0;iHit<nHits;++iHit) {
        TrackerHit * hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        trackExt->addTrackerHitExtended( hitExt );      
      }      
      
      
      
      _allTPCTracks.push_back( trackExt );                
    }      
  }
  catch ( DataNotAvailableException &e) {
    streamlog_out(DEBUG4) << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  // Reading Si Tracks
  try {
    LCCollection * col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG3) << std::endl;
    streamlog_out(DEBUG3) << "Number of Si Tracks = " << nelem << std::endl;
    streamlog_out(DEBUG3) << " Trk       p          D0         Z0       Px       Py       Pz   hitsSi ndf Chi2/ndf" << std::endl;
    
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * siTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( siTrack );
      TrackerHitVec hitVec = siTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(siTrack->getOmega());
      trackExt->setTanLambda(siTrack->getTanLambda());
      trackExt->setPhi(siTrack->getPhi());
      trackExt->setD0(siTrack->getD0());
      trackExt->setZ0(siTrack->getZ0());
      float cov[15];
      float param[5];
      
      param[0] = siTrack->getOmega();
      param[1] = siTrack->getTanLambda();
      param[2] = siTrack->getPhi();
      param[3] = siTrack->getD0();
      param[4] = siTrack->getZ0();      
      
      
      
      const FloatVec Cov = siTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      } 
      //      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack->getNdf());
      trackExt->setChi2(siTrack->getChi2());      
      char strg[200];
      HelixClass helixSi;
      for (int iHit=0;iHit<nHits;++iHit) {
        TrackerHit * hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        
        trackExt->addTrackerHitExtended( hitExt );      
      }
      
      
      float d0Si = trackExt->getD0();
      float z0Si = trackExt->getZ0();
      float omegaSi = trackExt->getOmega();
      float phi0Si = trackExt->getPhi();
      float tanLSi = trackExt->getTanLambda();
      helixSi.Initialize_Canonical(phi0Si,d0Si,z0Si,omegaSi,tanLSi,_bField);
      float pxSi = helixSi.getMomentum()[0];
      float pySi = helixSi.getMomentum()[1];
      float pzSi = helixSi.getMomentum()[2];
      const float pTot = sqrt(pxSi*pxSi+pySi*pySi+pzSi*pzSi);
      const int ndfSi = trackExt->getNDF();
      float Chi2Si = trackExt->getChi2()/float(trackExt->getNDF());
      sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
              pTot, d0Si,z0Si,pxSi,pySi,pzSi,nHits, ndfSi, Chi2Si);
      streamlog_out(DEBUG3) << strg << std::endl;
      
      if(nHits>0){
        _allSiTracks.push_back( trackExt );
      }else{
        delete trackExt;
      }
    }
    
    streamlog_out(DEBUG3) << std::endl;
  }
  catch ( DataNotAvailableException &e) {
    streamlog_out(DEBUG4) << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  
}

void FullLDCTracking_MarlinTrk::CleanUp(){
  
  int nNonCombTpc = int(_allNonCombinedTPCTracks.size());  
  for (int i=0;i<nNonCombTpc;++i) {
    TrackExtended * trkExt = _allNonCombinedTPCTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedTPCTracks.clear();
  
  int nNonCombSi = int(_allNonCombinedSiTracks.size());  
  for (int i=0;i<nNonCombSi;++i) {
    TrackExtended * trkExt = _allNonCombinedSiTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedSiTracks.clear();
  
  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();
  
  int nSETHits = int(_allSETHits.size());
  for (int i=0;i<nSETHits;++i) {
    TrackerHitExtended * hitExt = _allSETHits[i];
    delete hitExt;
  }
  _allSETHits.clear();
  
  int nTPCHits = int(_allTPCHits.size());
  for (int i=0;i<nTPCHits;++i) {
    TrackerHitExtended * hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();
  
  int nFTDHits = int(_allFTDHits.size());
  for (int i=0;i<nFTDHits;++i) {
    TrackerHitExtended * hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();
  
  int nETDHits = int(_allETDHits.size());
  for (int i=0;i<nETDHits;++i) {
    TrackerHitExtended * hitExt = _allETDHits[i];
    delete hitExt;
  }
  _allETDHits.clear();
  
  int nVTXHits = int(_allVTXHits.size());
  for (int i=0;i<nVTXHits;++i) {
    TrackerHitExtended * hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();
  
  int nSiTrk = int(_allSiTracks.size());
  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();
  
  int nTPCTrk = int(_allTPCTracks.size());
  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();
  
  int nCombTrk = int(_allCombinedTracks.size());
  for (int i=0;i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;    
  }
  _allCombinedTracks.clear();
  
  //   int nImplTrk = int(_trkImplVec.size());
  //   for (int i=0;i<nImplTrk;++i) {
  //     TrackExtended * trkImpl = _trkImplVec[i];
  //     delete trkImpl;
  //   }
  _trkImplVec.clear();
  
  //AS: Dont delete the individual entries, some of them are cleared elsewhere, I think
  _candidateCombinedTracks.clear();
}

void FullLDCTracking_MarlinTrk::MergeTPCandSiTracks() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      int iComp = 0;
      float angle = 0;
      float dOmega = CompareTrkII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angle);
      if ( (dOmega<_dOmegaForMerging) && (angle<_angleForMerging) && !VetoMerge(siTrackExt,tpcTrackExt)) {
        TrackExtended *combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);       
        if (combinedTrack != NULL) {
          _allCombinedTracks.push_back( combinedTrack );
          _candidateCombinedTracks.insert(tpcTrackExt);
          _candidateCombinedTracks.insert(siTrackExt);
          if (_debug >= 3 ) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }else{
          if (_debug >= 3 ) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
  
  
}



void FullLDCTracking_MarlinTrk::MergeTPCandSiTracksII() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    if(_candidateCombinedTracks.find(tpcTrackExt) != _candidateCombinedTracks.end() )continue;
    
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      if(_candidateCombinedTracks.find(siTrackExt)!= _candidateCombinedTracks.end() )continue;
      int iComp = 0;
      float angleSignificance = 0;
      float significance = CompareTrkIII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angleSignificance);
      if ( (significance<10) && (angleSignificance<5) && !VetoMerge(tpcTrackExt,siTrackExt) ) {
        TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);
        
        if (combinedTrack != NULL) {
          
          _allCombinedTracks.push_back( combinedTrack );
          if (_debug >= 3 ) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }else{
          if (_debug >= 3 ) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
}




TrackExtended * FullLDCTracking_MarlinTrk::CombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack, bool testCombinationOnly) {
  
  TrackExtended * OutputTrack = NULL;
  
  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();
  
  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());
  int nHits = nTPCHits + nSiHits;
  
  //std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nSiHits = " << nSiHits << std::endl;
  //std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nTPCHits = " << nTPCHits << std::endl;
  
  EVENT::TrackerHitVec trkHits;
  trkHits.reserve(nHits);
  
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::CombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }  
  
  for (int ih=0;ih<nTPCHits;++ih) {
    
    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::CombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }      
  
  double chi2_D;
  int ndf;
  
  if( trkHits.size() < 3 ) { 
    
    return 0 ;
    
  }
  
  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
  
  MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
  
  IMPL::TrackStateImpl pre_fit ;
  
  
  /** Provides the values of a track state from the first, middle and last hits in the hit_list. */
  int error = createPrefit( trkHits, &pre_fit, _bField, IMarlinTrack::backward);
  
  if ( error != IMarlinTrack::success ) {
    
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: creation of prefit fails with error " << error << std::endl;
    
    delete marlin_trk ;
    return 0;
    
  }
  
  // setup initial dummy covariance matrix
  EVENT::FloatVec covMatrix;
  covMatrix.resize(15);
  
  for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
    covMatrix[icov] = 0;
  }
  
  covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
  covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
  covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
  covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
  covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
  
  pre_fit.setCovMatrix(covMatrix);
  
  error = MarlinTrk::createFit( trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward , _maxChi2PerHit );
  
  if ( error != IMarlinTrack::success ) {
    
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: creation of fit fails with error " << error << std::endl;
    
    delete marlin_trk ;
    return 0;
    
  }
  
  
  const gear::Vector3D point(0.,0.,0.); // nominal IP
  int return_code = 0;
  
  TrackStateImpl trkState ;
  return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf ) ;
  
  if ( error != IMarlinTrack::success ) {
    
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: propagate to IP fails with error " << error << std::endl;
    
    delete marlin_trk ;
    return 0;
    
  }
  
  if ( ndf < 0  ) {
    
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Fit failed NDF is less that zero  " << ndf << std::endl;
    
    delete marlin_trk ;
    return 0;
    
  }
  
  
  float chi2Fit = chi2_D/float(ndf);
  
  if ( chi2Fit > _chi2FitCut ) {
    
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: track fail Chi2 cut of " << _chi2FitCut << " chi2 of track = " <<  chi2Fit << std::endl;
    
    delete marlin_trk ;
    return 0;
    
  }
  
  
  float omega = trkState.getOmega();
  float tanlambda = trkState.getTanLambda();
  float phi0 = trkState.getPhi();
  float d0 = trkState.getD0();
  float z0 = trkState.getZ0();    
  
  OutputTrack = new TrackExtended();
  GroupTracks * group = new GroupTracks();
  group->addTrackExtended(siTrack);
  group->addTrackExtended(tpcTrack);
  
  // note OutputTrack which is of type TrackExtended, only takes fits set for ref point = 0,0,0 
  OutputTrack->setGroupTracks(group);
  OutputTrack->setOmega(omega);
  OutputTrack->setTanLambda(tanlambda);
  OutputTrack->setPhi(phi0);
  OutputTrack->setZ0(z0);
  OutputTrack->setD0(d0);
  OutputTrack->setChi2(chi2_D);
  OutputTrack->setNDF(ndf);
  
  float cov[15];
  
  for (int i = 0 ; i<15 ; ++i) {
    cov[i] = trkState.getCovMatrix().operator[](i);
  }
  
  OutputTrack->setCovMatrix(cov);
  
  if ( testCombinationOnly == false ) {
    
    std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
    marlin_trk->getOutliers(outliers);
    
    
    for (int i=0;i<nSiHits;++i) {
      
      bool hit_is_outlier = false;
      
      for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
        if( outliers[ihit].first == siHitVec[i]->getTrackerHit() ){
          hit_is_outlier = true;
          break;
        }
      }
      
      if( hit_is_outlier == false ){
        TrackerHitExtended * hitExt = siHitVec[i];
        OutputTrack->addTrackerHitExtended(hitExt);
        hitExt->setUsedInFit(true);
      }
      
    }
    
    for (int i=0;i<nTPCHits;++i) {
      
      bool hit_is_outlier = false;
      
      for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
        if( outliers[ihit].first == tpcHitVec[i]->getTrackerHit() ){
          hit_is_outlier = true;
          break;
        }
      }
      
      if( hit_is_outlier == false ){
        TrackerHitExtended * hitExt = tpcHitVec[i];
        OutputTrack->addTrackerHitExtended(hitExt);
        hitExt->setUsedInFit(true);
      }
      
    }
  }
  
  delete marlin_trk ;
  
  return OutputTrack;
  
}



void FullLDCTracking_MarlinTrk::SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec) {
  
  int sizeOfVector = int(trackHitPairVec.size());
  TrackHitPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackHitPairVec[j];
      two = trackHitPairVec[j+1];
      float oneQ = one->getDistance();
      float twoQ = two->getDistance();
      if( oneQ > twoQ ) {
        Temp = trackHitPairVec[j];
        trackHitPairVec[j] = trackHitPairVec[j+1];
        trackHitPairVec[j+1] = Temp;
      }
    }  
  
  
}

void FullLDCTracking_MarlinTrk::Sorting(TrackExtendedVec & trackVec) {
  
  int sizeOfVector = int(trackVec.size());
  TrackExtended *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackVec[j];
      two = trackVec[j+1];
      float oneQ = one->getChi2()/float(one->getNDF());
      float twoQ = two->getChi2()/float(two->getNDF());
      if( oneQ > twoQ ) {
        Temp = trackVec[j];
        trackVec[j] = trackVec[j+1];
        trackVec[j+1] = Temp;
      }
    }  
}

void FullLDCTracking_MarlinTrk::SelectCombinedTracks() {
  
  int nCombTrk = int(_allCombinedTracks.size());
  
  for (int i=0; i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();
    int nTracks = int(tracks.size());
    if (nTracks == 2) {
      TrackExtended * firstTrack = tracks[0];
      TrackExtended * secondTrack = tracks[1];
      if ((firstTrack->getGroupTracks() == NULL) &&
          (secondTrack->getGroupTracks() == NULL) ) {
        firstTrack->setGroupTracks(group);
        secondTrack->setGroupTracks(group);     
        TrackerHitExtendedVec firstVec = firstTrack->getTrackerHitExtendedVec();
        TrackerHitExtendedVec secondVec = secondTrack->getTrackerHitExtendedVec();
        int nFirst = int(firstVec.size());
        int nSecond = int(secondVec.size());
        float edges[2];
        edges[0] = 1.0e+20;
        edges[1] = -1.0e+20;
        for (int iF=0;iF<nFirst;++iF) {
          TrackerHitExtended * trkHitExt = firstVec[iF];
          TrackerHit * trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;      
        }
        for (int iS=0;iS<nSecond;++iS) {
          TrackerHitExtended * trkHitExt = secondVec[iS];
          TrackerHit * trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;
        }
        group->setEdges(edges);
        _trkImplVec.push_back(trkExt);  
        
        if (_debug >= 3) {
          int iopt = 1;
          PrintOutMerging(secondTrack,firstTrack,iopt);
        }       
      }
    }else{
      if(nTracks>2) streamlog_out(DEBUG3) << " MORE THAN TWO TRACKS " << std::endl;
    }
  }
  
  
}

void FullLDCTracking_MarlinTrk::AddNotCombinedTracks() {  
  
  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());
  
  // we need some buffer vector
  TrackExtendedVec allMergedTracks;
  allMergedTracks.clear();
  
  // forcing merging of Si and TPC track segments
  if (_forceMerging==1) { 
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExtTPC = _allTPCTracks[i];
      GroupTracks * groupTPC = trkExtTPC->getGroupTracks();
      if (groupTPC == NULL) {
        float diffMin = 1.0e+20;  
        TrackExtended * siTrkToAttach = NULL;
        for (int j=0;j<nSiTrk;++j) {
          TrackExtended * trkExtSi = _allSiTracks[j];
          GroupTracks * groupSi = trkExtSi->getGroupTracks();
          if (groupSi == NULL) {
            int iComp = 0;
            //      float deltaP = CompareTrk(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp);
            float angle(0.);
            float angleSignificance(0.);
            
            float dOmega = CompareTrkII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angle);
            float significance = CompareTrkIII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angleSignificance);
            //      if (deltaP < _dPCutForForcedMerging) {
            if ( ((dOmega<_dOmegaForForcedMerging) && (angle<_angleForForcedMerging)) ||
                ((significance<5)                 && (angleSignificance<5))
                ) {
              float chi2O = dOmega/_dOmegaForForcedMerging;
              float chi2A = angle/_angleForForcedMerging;
              float deltaP = chi2O*chi2O + chi2A*chi2A; 
              if (deltaP<diffMin) {
                diffMin = deltaP;
                siTrkToAttach = trkExtSi;
              }
            }else{
              if (_debug==3) {
                int  iopt = 7;
                streamlog_out(DEBUG2) << significance << " " << angleSignificance << std::endl;
                PrintOutMerging(trkExtTPC,trkExtSi,iopt);
              }
            }
          }
        }
        
        if (siTrkToAttach!=NULL) {
          TrackExtended * trkExtSi = siTrkToAttach; 
          TrackExtended * OutputTrack = new TrackExtended();
          GroupTracks * group = new GroupTracks();
          group->addTrackExtended(trkExtSi);
          group->addTrackExtended(trkExtTPC);
          OutputTrack->setGroupTracks(group);
          //        trkExtSi->setGroupTracks(group);
          //        trkExtTPC->setGroupTracks(group);         
          OutputTrack->setOmega(trkExtTPC->getOmega());
          OutputTrack->setTanLambda(trkExtSi->getTanLambda());
          OutputTrack->setPhi(trkExtSi->getPhi());
          OutputTrack->setZ0(trkExtSi->getZ0());
          OutputTrack->setD0(trkExtSi->getD0());
          float covMatTPC[15];
          float covMatSi[15];
          float covMat[15];
          for (int iCov=0;iCov<15;++iCov) {
            covMatTPC[iCov] = trkExtTPC->getCovMatrix()[iCov];
            covMatSi[iCov] = trkExtSi->getCovMatrix()[iCov];                
            covMat[iCov] = covMatSi[iCov];
          }
          float scaling = sqrt(covMatTPC[5]/covMatSi[5]);
          covMat[5] = covMatTPC[5];
          covMat[3] = scaling*covMatSi[3];
          covMat[4] = scaling*covMatSi[4];
          covMat[8] = scaling*covMatSi[8];
          covMat[12] = scaling*covMatSi[12];          
          OutputTrack->setCovMatrix(covMat);
          TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
          TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();              
          int nTPCHits = int( tpcHitVec.size());
          int nSiHits = int( siHitVec.size());        
          float edges[2];
          edges[0] = 1.0e+20;
          edges[1] = -1.0e+20;
          for (int iH=0;iH<nSiHits;++iH) {
            TrackerHitExtended * hitExt = siHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true);
            TrackerHit * hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }       
          for (int iH=0;iH<nTPCHits;++iH) {
            TrackerHitExtended * hitExt = tpcHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true); 
            TrackerHit * hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }
          group->setEdges(edges);
          OutputTrack->setChi2(diffMin); // will be replaced if necessary
          OutputTrack->setNDF(int(1));   // will be replaced if necessary
          _allCombinedTracks.push_back( OutputTrack );
          allMergedTracks.push_back( OutputTrack );
        }           
      }
    }
    
    
    int nMerged = int(allMergedTracks.size());
    if (nMerged>0) {
      Sorting(allMergedTracks); 
      for (int iM=0;iM<nMerged;++iM) {
        TrackExtended * mergedTrack = allMergedTracks[iM];
        GroupTracks * grpTrk = mergedTrack->getGroupTracks();
        TrackExtendedVec trkVec = grpTrk->getTrackExtendedVec();
        TrackExtended * trkTPC = NULL;
        TrackExtended * trkSi = NULL;
        int nT = int(trkVec.size());
        if (nT==2) {
          trkTPC = trkVec[0];
          trkSi = trkVec[1];
          GroupTracks * groupTPC = trkTPC->getGroupTracks();
          GroupTracks * groupSi  = trkSi->getGroupTracks();
          if (groupTPC == NULL && groupSi == NULL) {
            trkTPC->setGroupTracks( grpTrk );
            trkSi->setGroupTracks( grpTrk );
            TrackerHitExtendedVec hitVec = mergedTrack->getTrackerHitExtendedVec();
            int nhits = int(hitVec.size());
            int totNdf = 2*nhits - 5;
            float totChi2 = trkTPC->getChi2() + trkSi->getChi2();
            mergedTrack->setNDF( totNdf );
            mergedTrack->setChi2( totChi2 );
            if (_debug >= 3) {
              int iopt = 2;
              PrintOutMerging(trkTPC,trkSi,iopt);
            }
            _trkImplVec.push_back( mergedTrack );
          }
        }
      }
    }
  }
  
  
  
  // clear buffer vector
  allMergedTracks.clear();
  
  // merging splitted TPC segments
  if (_mergeTPCSegments) {
    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();
    int nNonAssignedTPCSeg = 0;
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
        TrackerHitExtendedVec currentVec = trkExt->getTrackerHitExtendedVec();
        int nCur = int(currentVec.size());
        float zmin = 1e+20;
        float zmax = -1e+20;
        for (int iCur=0;iCur<nCur;++iCur) {
          TrackerHitExtended * curTrkHitExt = currentVec[iCur];
          TrackerHit * curTrkHit = curTrkHitExt->getTrackerHit();
          float zpos = float(curTrkHit->getPosition()[2]);
          if (zpos < zmin)
            zmin = zpos;
          if (zpos > zmax)
            zmax = zpos;
        }
        nNonAssignedTPCSeg++;
        int nGroups = int(TPCSegments.size());
        float dPtMin = 1.0e+10;
        GroupTracks * groupToAttach = NULL;
        TrackExtended * trkToAttach = NULL;
        for (int iG=0;iG<nGroups;++iG) {
          GroupTracks * segments = TPCSegments[iG];
          TrackExtendedVec segVec = segments->getTrackExtendedVec();
          int nTrk = int(segVec.size());
          bool consider = true;
          if (_forbidOverlapInZTPC==1) { // if overlap in Z of the two segments is forbidden
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trkInGroup = segVec[iTrk];
              TrackerHitExtendedVec hitInGroupVec = trkInGroup->getTrackerHitExtendedVec();
              int nHitsInGrp = int(hitInGroupVec.size());
              for (int iHitInGrp=0;iHitInGrp<nHitsInGrp;iHitInGrp++) {
                TrackerHitExtended * xTrkExt = hitInGroupVec[iHitInGrp];
                TrackerHit * xTrk = xTrkExt->getTrackerHit();
                float xZ = float(xTrk->getPosition()[2]);
                if (xZ>zmin&&xZ<zmax) {
                  consider = false;
                  break;
                }
              }
              if (!consider)
                break;
            }
          }
          if (consider) {
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trkInGroup = segVec[iTrk];
              int iComp = 1;
              float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);
              if (dPt < dPtMin && !VetoMerge(trkExt,trkInGroup)) {
                dPtMin = dPt;
                groupToAttach = segments;
                trkToAttach = trkInGroup;
                if (_debug==3) {
                  int iopt = 6;
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
              else {
                if (_debug==3) {
                  int iopt = 9;
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
          else {
            if (_debug >= 3) {
              int iopt = 9;
              for (int iTrk=0;iTrk<nTrk;++iTrk) {
                TrackExtended * trkInGroup = segVec[iTrk];
                int iComp = 1;
                float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);              
                if (dPt >= dPtMin) {          
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
        }
        if (dPtMin < _dPCutToMergeTPC && groupToAttach != NULL) {
          
          
          
          groupToAttach->addTrackExtended(trkExt);
          trkExt->setGroupTracks(groupToAttach);
          float zminGroup = groupToAttach->getEdges()[0]; 
          float zmaxGroup = groupToAttach->getEdges()[1];
          float edges[2];
          edges[0] = zmin;
          if (zminGroup<zmin)
            edges[0] = zminGroup;
          edges[1] = zmax;
          if (zmaxGroup>zmax)
            edges[1] = zmaxGroup;
          groupToAttach->setEdges(edges);
          if (_debug==3) {
            int iopt = 3;
            PrintOutMerging(trkExt,trkToAttach,iopt);
          }
        }
        else {
          GroupTracks * newSegment = new GroupTracks(trkExt);
          trkExt->setGroupTracks(newSegment);
          TPCSegments.push_back(newSegment);
          float edges[2];
          edges[0] = zmin;
          edges[1] = zmax;
          newSegment->setEdges(edges);
        }
      }
    }
    
    // combining splitted TPC segments with the reconstructed tracks having Si hits
    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());
    //    std::cout << "Combined tracks = " << nCombTrk << std::endl;
    //    std::cout << "nSegments = " << nSegments << std::endl;
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      TrackExtendedVec segVec = segments->getTrackExtendedVec();
      float zminTPCSeg = segments->getEdges()[0];
      float zmaxTPCSeg = segments->getEdges()[1];
      int nTrk = int(segVec.size());
      TrackExtended * CombTrkToAttach = NULL;
      TrackExtended * keyTrack = NULL;
      float deltaPtMin = _dPCutToMergeTPC;
      for (int iCTrk=0;iCTrk<nCombTrk;++iCTrk) {
        TrackExtended * combTrk = _trkImplVec[iCTrk];
        GroupTracks * groupComb = combTrk->getGroupTracks();
        bool consider = true;
        if (_forbidOverlapInZComb==1) { // if overlap in Z of the two segments is forbidden
          float zminComb = groupComb->getEdges()[0];
          float zmaxComb = groupComb->getEdges()[1];
          consider = (zminTPCSeg>zmaxComb) || (zmaxTPCSeg<zminComb);
        }
        if (consider) {
          for (int iTrk=0;iTrk<nTrk;++iTrk) {
            TrackExtended * trk = segVec[iTrk];
            int iopt = 0;
            float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
            float angleSignificance(0.);
            float significance = CompareTrkIII(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt,angleSignificance);
            if ( (dPt<deltaPtMin || significance <5 ) ) {
              if(VetoMerge(trk,combTrk)==false){
                CombTrkToAttach = combTrk;
                keyTrack = trk;
                deltaPtMin = dPt;
              }
            }
            else {
              if (_debug==3) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                PrintOutMerging(trk,dummySi,iopt_temp);
              }
            }
          }
        }
        else {
          if (_debug==3) {
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trk = segVec[iTrk];
              int iopt = 0;
              float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
              if (dPt>deltaPtMin) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                PrintOutMerging(trk,dummySi,iopt_temp);
              }
            }
          }
        }
      }
      
      if (CombTrkToAttach != NULL) { // attach TPC segment to existing Comb Track
        GroupTracks * groupToAttach = CombTrkToAttach->getGroupTracks();          
        TrackExtended * SiCombTrk = groupToAttach->getTrackExtendedVec()[0];
        TrackExtended * TpcCombTrk = groupToAttach->getTrackExtendedVec()[1];
        if (_debug==3) {
          int iopt = 4;
          PrintOutMerging(keyTrack,SiCombTrk,iopt);
          iopt = 5;
          PrintOutMerging(keyTrack,TpcCombTrk,iopt);      
        }
        for (int iTrk=0;iTrk<nTrk;iTrk++) {
          TrackExtended * trk = segVec[iTrk];
          groupToAttach->addTrackExtended( trk );
          trk->setGroupTracks( groupToAttach );
          TrackerHitExtendedVec hitVec = trk->getTrackerHitExtendedVec();
          int nHitS = int(hitVec.size());              
          for (int iHS=0;iHS<nHitS;++iHS) {
            TrackerHitExtended * hitExt = hitVec[iHS];
            hitExt->setUsedInFit(false);
            CombTrkToAttach->addTrackerHitExtended( hitExt );
          }
        }
      }
      else {
        if (nTrk==1) { //
          GroupTracks * newGrp = new GroupTracks();
          segVec[0]->setGroupTracks(newGrp);
          newGrp->addTrackExtended(segVec[0]);
          TrackerHitExtendedVec TpcHitVec = segVec[0]->getTrackerHitExtendedVec();
          int nTpcH = int(TpcHitVec.size());
          for (int iTpcH=0;iTpcH<nTpcH;++iTpcH) {
            TpcHitVec[iTpcH]->setUsedInFit( true );
          }
          _trkImplVec.push_back(segVec[0]);
          _allNonCombinedTPCTracks.push_back(segVec[0]);
        }
        else {
          float zMin = 1.0e+20;
          TrackExtended * chosenTrack = NULL;
          for (int iTrk=0;iTrk<nTrk;++iTrk) {
            TrackExtended * trk = segVec[iTrk];
            Track * track = trk->getTrack();
            TrackerHitVec hitVec = track->getTrackerHits();
            int nHits = int(hitVec.size());
            for (int iH=0;iH<nHits;++iH) {
              float zPosi = fabs(hitVec[iH]->getPosition()[2]);
              if (zPosi<zMin) {
                chosenTrack = trk;
                zMin = zPosi;
                break;
              }
            }
          }
          if (chosenTrack!=NULL) {
            GroupTracks * newGroup = new GroupTracks();
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trk = segVec[iTrk];
              trk->setGroupTracks( newGroup );
              newGroup->addTrackExtended( trk );                        
              TrackerHitExtendedVec hitVecS = trk->getTrackerHitExtendedVec();
              int nHitS = int(hitVecS.size());                  
              for (int iH=0;iH<nHitS;++iH) {
                TrackerHitExtended * trkHitExt = hitVecS[iH];
                if (trk!=chosenTrack) {
                  trkHitExt->setUsedInFit( false );
                  chosenTrack->addTrackerHitExtended( trkHitExt );                              
                }
                else {
                  trkHitExt->setUsedInFit( true );
                }
              }
            }
            _allNonCombinedTPCTracks.push_back(chosenTrack);
            _trkImplVec.push_back(chosenTrack);
          }
        }
      }
    }
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      delete segments;
    }
    TPCSegments.clear();
  }
  else { // adding all TPC segments to the list of tracks (track splitting is allowed)
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      Track * track = trkExt->getTrack();
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
        TrackerHitVec hitVec = track->getTrackerHits();       
        _trkImplVec.push_back(trkExt);
        _allNonCombinedTPCTracks.push_back( trkExt );
        GroupTracks * newGrp = new GroupTracks();
        newGrp->addTrackExtended( trkExt );
        trkExt->setGroupTracks( newGrp );
      }
    }    
  }
  
  for (int i=0;i<nSiTrk;++i) { // adding left-over Si segments to the list of tracks
    TrackExtended * trkExt = _allSiTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    if (group == NULL) {
      TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
      int nHSi = int(hitVec.size());
      for (int iHSi=0;iHSi<nHSi;++iHSi) {
        hitVec[iHSi]->setUsedInFit(true);
      }
      _trkImplVec.push_back(trkExt);
      GroupTracks * newGrp = new GroupTracks();
      newGrp->addTrackExtended( trkExt );
      trkExt->setGroupTracks( newGrp );   
      _allNonCombinedSiTracks.push_back( trkExt );
    }
  }
  
}

void FullLDCTracking_MarlinTrk::CheckTracks() {  
  
  for(unsigned int i = 0; i< _trkImplVec.size();i++){
    TrackExtended *first = _trkImplVec[i];
    if(first==NULL)continue;
    float d0First = first->getD0();
    float z0First = first->getZ0();
    float omegaFirst = first->getOmega();
    float tanLFirst = first->getTanLambda();
    float phiFirst = first->getPhi();
    HelixClass helixFirst;
    helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float momFirst[3];
    momFirst[0]= helixFirst.getMomentum()[0];
    momFirst[1]= helixFirst.getMomentum()[1];
    momFirst[2]= helixFirst.getMomentum()[2];
    float pFirst    = sqrt(momFirst[0]*momFirst[0]+momFirst[1]*momFirst[1]+momFirst[2]*momFirst[2]);
    if(std::isnan(pFirst))continue;
    TrackerHitExtendedVec firstHitVec  = first->getTrackerHitExtendedVec();
    if(firstHitVec.size()<1)continue;
    
    for(unsigned int j = i+1; j<_trkImplVec.size();j++){
      TrackExtended *second = _trkImplVec[j];
      if(second==NULL)continue;
      float d0Second = second->getD0();
      float z0Second = second->getZ0();
      float omegaSecond = second->getOmega();
      float tanLSecond = second->getTanLambda();
      float phiSecond = second->getPhi();
      HelixClass helixSecond;
      helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
      float momSecond[3];
      momSecond[0] = helixSecond.getMomentum()[0];
      momSecond[1] = helixSecond.getMomentum()[1];
      momSecond[2] = helixSecond.getMomentum()[2];
      float pSecond    = sqrt(momSecond[0]*momSecond[0]+momSecond[1]*momSecond[1]+momSecond[2]*momSecond[2]);
      if(std::isnan(pSecond))continue;
      TrackerHitExtendedVec secondHitVec  = second->getTrackerHitExtendedVec();
      if(secondHitVec.size()<1)continue;
      if(firstHitVec.size()+secondHitVec.size()<10)continue;
      
      
      float pdot = (momFirst[0]*momSecond[0]+momFirst[1]*momSecond[1]+momFirst[2]*momSecond[2])/pFirst/pSecond;
      if(pdot<0.999)continue;
      // const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      // const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      // const float deltaP = fabs(pFirst-pSecond);
      // const float sigmaDeltaP = sqrt(pFirst*sigmaPOverPFirst*pFirst*sigmaPOverPFirst+pSecond*sigmaPOverPSecond*pSecond*sigmaPOverPSecond);
      //      const float significance = deltaP/sigmaDeltaP;
      
      TrackExtended * combinedTrack = CombineTracks(first,second, true);
      if(combinedTrack != NULL){
        const int minHits = std::min(firstHitVec.size(),secondHitVec.size());
        const int maxHits = std::max(firstHitVec.size(),secondHitVec.size());
        
        if( combinedTrack->getNDF() <= 2*maxHits+minHits-5){
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
          continue;
        }
        
        float d0    = combinedTrack->getD0();
        float z0    = combinedTrack->getZ0();
        float omega = combinedTrack->getOmega();
        float tanL  = combinedTrack->getTanLambda();
        float phi   = combinedTrack->getPhi();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi,d0,z0,omega,tanL,_bField);
        float mom[3];
        mom[0]  = helix.getMomentum()[0];
        mom[1]  = helix.getMomentum()[1];
        mom[2]  = helix.getMomentum()[2];
        // float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
        // float chi2Sig =  (combinedTrack->getChi2() - combinedTrack->getNDF());
        // chi2Sig = chi2Sig/sqrt(combinedTrack->getNDF()*2);
        
        
        
        int nTpcFirst(0);
        int nUsedFirst(0);
        for(unsigned int ihit = 0;ihit<firstHitVec.size();ihit++){
          
          if( getDetectorID(firstHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) nTpcFirst++;
          
          if(firstHitVec[ihit]->getUsedInFit()==true)nUsedFirst++;
        }
        
        
        int nTpcSecond(0);
        int nUsedSecond(0);
        for(unsigned int ihit = 0;ihit<secondHitVec.size();ihit++){
          if( getDetectorID(secondHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) ++nTpcSecond;
          if( secondHitVec[ihit]->getUsedInFit()==true ) ++nUsedSecond;
        }
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }
    }
  }
  
  
}





float FullLDCTracking_MarlinTrk::CompareTrkII(TrackExtended * first, TrackExtended * second, 
                                              float d0Cut, float z0Cut,int iopt,float & Angle) {
  
  
  float result = 1.0e+20;
  Angle  = 1.0e+20; 
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();
  float deltaOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
  if(deltaOmega> 2*_dOmegaForMerging)return result;
  
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  if(!isCloseInIP)return result;
  
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  Angle = acos(Angle);
  
  result = deltaOmega;
  
  return result;
  
}




float FullLDCTracking_MarlinTrk::CompareTrkIII(TrackExtended * first, TrackExtended * second, 
                                               float d0Cut, float z0Cut,int iopt, float & AngleSignificance) {
  
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  
  //MB 2010 03
  float d0ErrFirst = sqrt(first->getCovMatrix()[0]);
  float z0ErrFirst = sqrt(first->getCovMatrix()[9]);
  //  float omegaErrFirst = sqrt(first->getCovMatrix()[5]);
  float phiErrFirst = sqrt(first->getCovMatrix()[2]);
  float qErrFirst = sqrt(cos(qFirst)*cos(qFirst)*first->getCovMatrix()[14]);
  //MB END
  //MB 2010 03
  float d0ErrSecond = sqrt(second->getCovMatrix()[0]);
  float z0ErrSecond = sqrt(second->getCovMatrix()[9]);
  //  float omegaErrSecond = sqrt(second->getCovMatrix()[5]);
  float phiErrSecond = sqrt(second->getCovMatrix()[2]);
  float qErrSecond = sqrt(cos(qSecond)*cos(qSecond)*second->getCovMatrix()[14]);
  //MB END
  
  
  //  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  //isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  //MB 2010 03
  bool isCloseInIP = (fabs(d0First-d0Second)/sqrt(d0ErrFirst*d0ErrFirst+d0ErrSecond*d0ErrSecond)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)/sqrt(z0ErrFirst*z0ErrFirst+z0ErrSecond*z0ErrSecond)<z0Cut);
  
  if (!isCloseInIP)return result;
  
  float Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float momFirst = 0;
  float momSecond = 0;
  
  for (int iC=0;iC<3;++iC) {
    pFirst[iC] = helixFirst.getMomentum()[iC];
    pSecond[iC] = helixSecond.getMomentum()[iC];
    momFirst += pFirst[iC]* pFirst[iC];
    momSecond += pSecond[iC]*pSecond[iC];
  }
  momFirst = sqrt(momFirst);
  momSecond = sqrt(momSecond);
  
  
  float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
  
  
  const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
  const float deltaP = fabs(momFirst-momSecond);
  const float sigmaPFirst = momFirst*sigmaPOverPFirst;
  const float sigmaPSecond = momSecond*sigmaPOverPSecond;
  const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
  const float significance = deltaP/sigmaDeltaP;
  
  //MB 2010 03
  float errorAngle =sin(phiFirst)*sin(phiFirst)*phiErrFirst*phiErrFirst*cos(phiSecond)*cos(phiSecond)+
  sin(phiSecond)*sin(phiSecond)*phiErrSecond*phiErrSecond*cos(phiFirst)*cos(phiFirst)+
  sin(qFirst)*sin(qFirst)*qErrFirst*qErrFirst*cos(qSecond)*cos(qSecond)+
  sin(qSecond)*sin(qSecond)*qErrSecond*qErrSecond*cos(qFirst)*cos(qFirst)+
  cos(phiFirst)*cos(phiFirst)*phiErrFirst*phiErrFirst*(sin(phiSecond)*sin(qFirst)*sin(qSecond))*(sin(phiSecond)*sin(qFirst)*sin(qSecond))+
  cos(phiSecond)*cos(phiSecond)*phiErrSecond*phiErrSecond*(sin(phiFirst)*sin(qFirst)*sin(qSecond))*(sin(phiFirst)*sin(qFirst)*sin(qSecond))+
  cos(qFirst)*cos(qFirst)*qErrFirst*qErrFirst*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))+
  cos(qSecond)*cos(qSecond)*qErrSecond*qErrSecond*(sin(phiFirst)*sin(phiSecond)*sin(qFirst))*(sin(phiFirst)*sin(phiSecond)*sin(qFirst));
  
  if(Angle<1.){
    errorAngle = sqrt(1./(1.-Angle*Angle)*errorAngle);
  }else{
    errorAngle = sqrt(errorAngle);
  }
  
  if(errorAngle<1.e-6)errorAngle=1.e-6;
  
  AngleSignificance = fabs(acos(Angle)/errorAngle);
  if(pdot<0.999)return result;
  
  return significance;
  
}


float FullLDCTracking_MarlinTrk::CompareTrk(TrackExtended * first, TrackExtended * second, 
                                            float d0Cut, float z0Cut,int iopt) {
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  
  if (iopt>0) isCloseInIP = isCloseInIP || (fabs(d0First+d0Second)<d0Cut);
  
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float dPminus[3];
  float dPplus[3];
  float momFirst = 0;
  float momSecond = 0;
  float momMinus = 0;
  float momPlus = 0;
  
  if ( isCloseInIP ) {
    
    for (int iC=0;iC<3;++iC) {
      pFirst[iC] = helixFirst.getMomentum()[iC];
      pSecond[iC] = helixSecond.getMomentum()[iC];
      momFirst += pFirst[iC]* pFirst[iC];
      momSecond += pSecond[iC]*pSecond[iC];
      dPminus[iC] = pFirst[iC] - pSecond[iC];
      dPplus[iC] = pFirst[iC] + pSecond[iC];
      momMinus += dPminus[iC]*dPminus[iC];
      momPlus += dPplus[iC]*dPplus[iC];
    }
    momFirst = sqrt(momFirst);
    momSecond = sqrt(momSecond);
    
    float ptFirst = sqrt(pFirst[0]*pFirst[0]+pFirst[1]*pFirst[1]);
    float ptSecond = sqrt(pSecond[0]*pSecond[0]+pSecond[1]*pSecond[1]);
    
    
    if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {
      
      momMinus = sqrt(momMinus);
      momPlus = sqrt(momPlus);
      float nom = momMinus;
      if (momPlus<nom && iopt>0)
        nom = momPlus;
      float den = momFirst;
      if (momSecond<momFirst)
        den = momSecond;
      
      result = nom/den;     
      
    }
    
    else {
      
      
      // check for cases where PatRec splits non-looping TPC tracks 
      // look for two tracks where total tpc hits are not more than total number
      // of pad rows and that the hits on one track are close to the helix of the
      // other track
      
      float dpOverP = 2.0*fabs(momFirst-momSecond)/(momFirst+momSecond);
      const float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
      const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      const float deltaP = fabs(momFirst-momSecond);
      const float sigmaPFirst = momFirst*sigmaPOverPFirst;
      const float sigmaPSecond = momSecond*sigmaPOverPSecond;
      const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
      const float significance = deltaP/sigmaDeltaP;
      
      
      //compare angle between the two vectors (cos theta) and their momentum
      //    if( ( pdot>0.99 && dpOverP<0.01 ) || ( pdot>0.998 && dpOverP<0.25 ) ){
      if( (pdot>_cosThetaCutHighPtMerge && dpOverP<_momDiffCutHighPtMerge) 
         || 
         (pdot>_cosThetaCutSoftHighPtMerge && dpOverP<_momDiffCutSoftHighPtMerge) 
         || (pdot > 0.9999 && significance <10) 
         ){
        
        
        int nTrkGrpFirst = 0;
        int nTrkGrpSecond = 0;
        TrackerHitVec hitvecFirst;
        TrackerHitVec hitvecSecond;
        GroupTracks * groupFirst = first->getGroupTracks();
        GroupTracks * groupSecond = second->getGroupTracks();
        
        if(groupFirst!=NULL){
          
          TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
          nTrkGrpFirst = int(tracksInGroupFirst.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
            
            TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
            TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i =0; i<hitVec.size(); ++i){
              hitvecFirst.push_back(hitVec[i]->getTrackerHit());          
            }
          }
        }
        
        if(groupSecond!=NULL){
          
          TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
          nTrkGrpSecond = int(tracksInGroupSecond.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
            TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
            TrackerHitExtendedVec hitVec = 
            trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i=0;i<hitVec.size();++i){
              hitvecSecond.push_back(hitVec[i]->getTrackerHit());
            }
          }
        }
        
        
        
        // for non-looping tracks 
        int nhitsFirst  = (int)hitvecFirst.size();
        int nhitsSecond = (int)hitvecSecond.size();
        int ntpcFirst   = 0;
        int ntpcSecond  = 0;
        float hitxyz[3];
        float dist[3];
        float maxdistFirst=0;
        float maxdistSecond=0;
        int ncloseFirst = 0;
        int ncloseSecond = 0;
        float zminFirst = 99999;
        float zminSecond = 99999;
        float zmaxFirst = -99999;
        float zmaxSecond = -99999;
        
        for(int ih =0;ih<nhitsFirst;++ih){
          
          float x = (float) hitvecFirst[ih]->getPosition()[0];
          float y = (float) hitvecFirst[ih]->getPosition()[1];
          float z = (float) hitvecFirst[ih]->getPosition()[2];
          
          if(fabs(z)<zminFirst) zminFirst=fabs(z);
          if(fabs(z)>zmaxFirst) zmaxFirst=fabs(z);
          
          float r = sqrt(x*x+y*y);
          
          if(r>_tpc_inner_r) ntpcFirst++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixSecond.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistFirst) maxdistFirst=dist[2];
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseFirst++;
        }
        
        for(int ih =0;ih<nhitsSecond;++ih){
          
          float x = (float)hitvecSecond[ih]->getPosition()[0];
          float y = (float)hitvecSecond[ih]->getPosition()[1];
          float z = (float)hitvecSecond[ih]->getPosition()[2];
          
          if(fabs(z)<zminSecond) zminSecond=fabs(z);
          if(fabs(z)>zmaxSecond) zmaxSecond=fabs(z);
          
          float r = sqrt(x*x+y*y);
          
          if(r>_tpc_inner_r) ntpcSecond++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixFirst.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistSecond) maxdistSecond=dist[2];
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseSecond++;
        }
        
        float fcloseFirst  = (float)ncloseFirst/(float)nhitsFirst;
        float fcloseSecond = (float)ncloseSecond/(float)nhitsSecond;
        
        
        
        bool split = false;
        //std::cout << "Momenta = " << momFirst << " " << momSecond << std::endl;
        //std::cout << "MaxDist = " << maxdistSecond << " " << maxdistFirst << " " << _maxHitDistanceCutHighPtMerge << std::endl;
        //std::cout << "close   = " << fcloseSecond << " " << fcloseFirst << " " << _maxFractionOfOutliersCutHighPtMerge << std::endl;
        //std::cout << "ntpc    = " << ntpcFirst << " " << ntpcSecond << " " << _tpc_pad_height+10 << std::endl;
        //std::cout << "pdot    = " << pdot << " significance " << significance << std::endl;
        
        TrackExtended * combinedTrack = CombineTracks(first,second, true);
        
        if(combinedTrack != NULL){
          //std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << first->getNDF()+second->getNDF()+5 << std::endl;
          if(combinedTrack->getNDF()+10>first->getNDF()+second->getNDF()+5){
            split = true;
            dpOverP = 0;
            //std::cout << " Forcing MERGE " << std::endl;
          }
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
        } 
        else {
          //std::cout << "Could not combine track " << std::endl;
          if(significance<5 && fcloseFirst>_maxFractionOfOutliersCutHighPtMerge){
            split = true;
            dpOverP = 0;
            //      int overlap = SegmentRadialOverlap(first,second);
            //std::cout << " Forcing MERGE " << overlap << std::endl;
          }
        }
        
        // criteria for split track
        // old criterion
        if( maxdistSecond < _maxHitDistanceCutHighPtMerge && maxdistFirst < _maxHitDistanceCutHighPtMerge 
           && 
           (fcloseSecond > _maxFractionOfOutliersCutHighPtMerge || fcloseFirst > _maxFractionOfOutliersCutHighPtMerge) 
           && 
           ntpcFirst+ntpcSecond < _tpc_pad_height+10.) {
          
          split = true;        
          
        }
        
        if(split){
          result = dpOverP;
        }
        
      }
    }
  }
  
  return result;
  
}

void FullLDCTracking_MarlinTrk::AddNotAssignedHits() {
  
  
  // Creating helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CreateExtrapolations();
  
  //  if (_assignSETHits>0) { // Assignment of SET Hits
  //    
  //    const gear::GearParameters& pSETDet = Global::GEAR->getGearParameters("SET");  
  //    int nLayersSET = int(pSETDet.getDoubleVals("SETLayerRadius").size());
  //    
  //    int nSETHits = _allSETHits.size();
  //    std::vector<TrackerHitExtendedVec> SETHits;
  //    SETHits.resize(nLayersSET);
  //    
  //    for (int iSET=0;iSET<nSETHits;++iSET) {
  //      TrackerHitExtended * trkHit = _allSETHits[iSET];
  //      TrackerHit * hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0&&layer<nLayersSET) 
  //        SETHits[layer].push_back(trkHit);   
  //    }
  //    for (int iL=0; iL<nLayersSET; ++iL) { // loop over SET layers
  //      TrackerHitExtendedVec hitVec = SETHits[iL];
  //      int refit = 1;
  //      AssignOuterHitsToTracks(hitVec,_distCutForSETHits,refit);
  //    }
  //  }
  
  //  if (_assignETDHits>0) { // Assignment of ETD Hits
  //    
  //    const gear::GearParameters& pETDDet = Global::GEAR->getGearParameters("ETD");  
  //    int nLayersETD = int(pETDDet.getDoubleVals("ETDLayerZ").size());
  //    
  //    int nETDHits = _allETDHits.size();
  //    std::vector<TrackerHitExtendedVec> ETDHits;
  //    ETDHits.resize(nLayersETD);
  //    
  //    for (int iETD=0;iETD<nETDHits;++iETD) {
  //      TrackerHitExtended * trkHit = _allETDHits[iETD];
  //      TrackerHit * hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0 && layer < nLayersETD) 
  //        ETDHits[layer].push_back(trkHit);
  //    }
  //    for (int iL=0; iL<nLayersETD; ++iL) { // loop over ETD layers
  //      TrackerHitExtendedVec hitVec = ETDHits[iL];
  //      int refit = 0;
  //      AssignOuterHitsToTracks( hitVec, _distCutForETDHits, refit );
  //    }
  //    
  //  }
  
  //  // Cleaning up helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CleanUpExtrapolations();
  
  if (_assignSITHits>0) { // Treatment of left-over SIT hits 
    
    streamlog_out(DEBUG4) << "Assign SIT hits *********************************" << std::endl;
    
    std::vector<TrackerHitExtendedVec> nonAssignedSITHits;    
    nonAssignedSITHits.resize(_nLayersSIT);
    
    int nSITHits = int(_allSITHits.size());    
    
    // loop over all SIT hits ...
    for (int iH=0;iH<nSITHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allSITHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
        
        int layer = getLayerID(trkHit);
        
        if (layer >=0 && (unsigned)layer < _nLayersSIT) {
          nonAssignedSITHits[layer].push_back(trkHitExt);
        }
      }
    }       
    
    for (int iL=_nLayersSIT-1;iL>=0;--iL) { // reverse loop over layers in Si
      
      TrackerHitExtendedVec hitVec = nonAssignedSITHits[iL];
      
      if ( hitVec.empty() == false ) {
        streamlog_out(DEBUG4) << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec,
                             _distCutForSITHits);
        
      }
      
    }
  }
  
  if (_assignFTDHits>0) { // Treatment of left-over FTD hits
    
    streamlog_out(DEBUG4) << "Assign FTD hits *********************************" << std::endl;
    
    std::vector<TrackerHitExtendedVec> nonAssignedFTDHits;
    nonAssignedFTDHits.resize(_nLayersFTD);
    
    int nFTDHits = int(_allFTDHits.size());
    
    // loop over all FTD hits ...
    for (int iH=0;iH<nFTDHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allFTDHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
        
        
        // get the layer number
        int layer = getLayerID(trkHit);
        int petalIndex = getModuleID(trkHit);
        
        if ( _petalBasedFTDWithOverlaps == true ) {
          
          // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
          // the layers will follow the even odd numbering of the petals 
          if ( petalIndex % 2 == 0 ) {
            layer = 2*layer;
          }
          else {
            layer = 2*layer + 1;
          }
          
        }
        
        if (layer >=0 && layer < (int)_nLayersFTD)
          nonAssignedFTDHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersFTD-1;iL>=0;--iL) {
      if ( nonAssignedFTDHits[iL].empty() == false ) {
        
        TrackerHitExtendedVec hitVec = nonAssignedFTDHits[iL];
        streamlog_out(DEBUG4) << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec,
                             _distCutForFTDHits);     
        
      }
    }
  }
  
  
  
  if (_assignVTXHits>0) { // Treatment of left-over VTX hits
    
    streamlog_out(DEBUG4) << "Assign VXD hits *********************************" << std::endl;
    
    std::vector<TrackerHitExtendedVec> nonAssignedVTXHits;
    nonAssignedVTXHits.resize(_nLayersVTX);
    
    int nVTXHits = int(_allVTXHits.size());
    
    // loop over all VXD hits ...
    for (int iH=0;iH<nVTXHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allVTXHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
        
        int layer = getLayerID(trkHit);
        
        if (layer >=0 && layer < (int)_nLayersVTX)
          nonAssignedVTXHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersVTX-1;iL>=0;--iL) {
      if ( nonAssignedVTXHits[iL].empty() == false ) {
        
        TrackerHitExtendedVec hitVec = nonAssignedVTXHits[iL];
        streamlog_out(DEBUG4) << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec,
                             _distCutForVTXHits);     
      }
    }
  }
  
  streamlog_out(DEBUG4) << "Assign TPC hits *********************************" << std::endl;
  
  if (_assignTPCHits) {// Treatment of left-over TPC hits
    TrackerHitExtendedVec nonAssignedTPCHits;
    int nTPCHits = int(_allTPCHits.size());
    for (int iH=0;iH<nTPCHits;++iH) {
      TrackerHitExtended * trkHitExt = _allTPCHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      if (trkExt == NULL) {
        nonAssignedTPCHits.push_back(trkHitExt);
      }
    }
    streamlog_out(DEBUG4) << "AddNotAssignedHits : Number of Non Assigned TPC hits = " <<  nonAssignedTPCHits.size() << std::endl;
    AssignTPCHitsToTracks(nonAssignedTPCHits,
                          _distCutForTPCHits);
  }
  
  
}


void FullLDCTracking_MarlinTrk::CreateExtrapolations() {
  
  _trackExtrapolatedHelix.clear();
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix = GetExtrapolationHelix( trk );
    _trackExtrapolatedHelix[trk] = helix;
  }
  
}

void FullLDCTracking_MarlinTrk::CleanUpExtrapolations() {
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix =  _trackExtrapolatedHelix[trk];
    delete helix;
  }  
  _trackExtrapolatedHelix.clear();
  
}


//void FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit) {
//  
//  int nHits = int(hitVec.size());
//  int nTrk = int(_trkImplVec.size());
//  
//  std::map <TrackExtended*,bool> flagTrack;
//  std::map <TrackerHitExtended*,bool> flagHit;
//  TrackHitPairVec pairs;
//  flagTrack.clear();
//  flagHit.clear();
//  pairs.clear();
//  
//  for (int iH=0;iH<nHits;++iH) {
//    float pos[3];
//    TrackerHitExtended * trkHitExt = hitVec[iH];
//    TrackerHit * hit = trkHitExt->getTrackerHit();
//    for (int ip=0;ip<3;++ip)
//      pos[ip] = float(hit->getPosition()[ip]);
//    //      float r_hit = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
//    //      float phi_hit = atan2(pos[1],pos[0]);
//    //      float z_hit = pos[2];
//    //      std::cout << r_hit << " " << phi_hit << " "
//    //                << z_hit << " " 
//    //                << trkHitExt->getResolutionRPhi() << " "
//    //                << trkHitExt->getResolutionZ() << std::endl;
//    
//    for (int iT=0;iT<nTrk;++iT) {
//      TrackExtended * trkExt = _trkImplVec[iT];
//      float tanLambda = trkExt->getTanLambda();           
//      float product = pos[2]*tanLambda;
//      if (product>0) {
//        HelixClass * helix = _trackExtrapolatedHelix[trkExt];
//        float distance = helix->getDistanceToPoint(pos,dcut);
//        //      std::cout << "Dist = " << dist[2] << std::endl;
//        if (distance<dcut) {
//          TrackHitPair * trkHitPair = 
//          new TrackHitPair(trkExt,trkHitExt,distance);
//          pairs.push_back(trkHitPair);
//          flagTrack[trkExt] = true;
//          flagHit[trkHitExt] = true;
//        }
//      }
//    }
//  }
//  
//  int nPairs = int(pairs.size());
//  if (nPairs>0) {
//    SortingTrackHitPairs(pairs);
//    for (int iP=0;iP<nPairs;++iP) {
//      TrackHitPair * trkHitPair = pairs[iP];
//      TrackExtended * trkExt = trkHitPair->getTrackExtended();
//      TrackerHitExtended * trkHitExt = 
//      trkHitPair->getTrackerHitExtended();       
//      if (flagTrack[trkExt] && flagHit[trkHitExt]) {
//        if (refit==0) {
//          trkExt->addTrackerHitExtended( trkHitExt );
//          trkHitExt->setUsedInFit( false );
//          trkHitExt->setTrackExtended( trkExt );
//        }
//        else {
//          TrackerHitExtendedVec hitsInTrack = 
//          trkExt->getTrackerHitExtendedVec();
//          int nTotH = int(hitsInTrack.size());
//          int nHitsInFit = 0;
//          for (int iTH=0;iTH<nTotH;++iTH) {
//            TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
//            if (hitInTrack->getUsedInFit())
//              nHitsInFit++;
//          }
//          float * x_h = new float[nHitsInFit+1];
//          float * y_h = new float[nHitsInFit+1];
//          float * z_h = new float[nHitsInFit+1];
//          int * idet_h = new int[nHitsInFit+1];
//          int * ityp_h = new int[nHitsInFit+1];
//          int * lhits = new int[nHitsInFit+1];
//          float * rR_h = new float[nHitsInFit+1];
//          float * rZ_h = new float[nHitsInFit+1];
//          int iHitInFit = 0;
//          for (int iHit=0;iHit<nTotH;++iHit) {
//            TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
//            if (hitInTrack->getUsedInFit()) {
//              TrackerHit * hit = hitInTrack->getTrackerHit();
//              x_h[iHitInFit] = float(hit->getPosition()[0]);
//              y_h[iHitInFit] = float(hit->getPosition()[1]);
//              z_h[iHitInFit] = float(hit->getPosition()[2]);
//              idet_h[iHitInFit] = hitInTrack->getDet();
//              ityp_h[iHitInFit] = hitInTrack->getType();
//              rR_h[iHitInFit] = hitInTrack->getResolutionRPhi();
//              rZ_h[iHitInFit] = hitInTrack->getResolutionZ();
//              iHitInFit++;
//            }
//          }
//          TrackerHit * remainHit = trkHitExt->getTrackerHit();
//          x_h[iHitInFit] = float(remainHit->getPosition()[0]);
//          y_h[iHitInFit] = float(remainHit->getPosition()[1]);
//          z_h[iHitInFit] = float(remainHit->getPosition()[2]);
//          idet_h[iHitInFit] = trkHitExt->getDet();
//          ityp_h[iHitInFit] = trkHitExt->getType();
//          rR_h[iHitInFit] = trkHitExt->getResolutionRPhi();
//          rZ_h[iHitInFit] = trkHitExt->getResolutionZ();          
//          iHitInFit++;
//          
//          int NPT = iHitInFit;
//          float chi2_D;
//          int ndf_D;
//          float chi2rphi,chi2z;
//          float par[5];
//          float epar[15];
//          float refPoint[3];
//          
//          int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,
//                                         _bField,idet_h,ityp_h,
//                                         _chi2PrefitCut,
//                                         x_h,y_h,z_h,rR_h,rZ_h,
//                                         par,epar,refPoint,chi2_D,ndf_D,
//                                         chi2rphi,chi2z,lhits);              
//          
//          float chi2ndf = chi2_D/float(ndf_D);
//          if (ierr==0 && chi2ndf<_chi2FitCut) {
//            trkExt->addTrackerHitExtended( trkHitExt );
//            trkExt->setCovMatrix(epar);
//            trkExt->setOmega(par[0]);   
//            trkExt->setTanLambda(par[1]);
//            trkExt->setPhi(par[2]);
//            trkExt->setD0(par[3]);
//            trkExt->setZ0(par[4]);
//            trkExt->setChi2(chi2_D);
//            trkExt->setNDF(ndf_D);
//            trkHitExt->setTrackExtended( trkExt );
//            trkHitExt->setUsedInFit( true );
//            flagTrack[trkExt] = false;
//            flagHit[trkHitExt] = false;
//            if (_optFit == 4) {
//              GroupTracks * groupTracks = trkExt->getGroupTracks();
//              if (groupTracks!=NULL) {
//                TrackExtendedVec group = groupTracks->getTrackExtendedVec();
//                int nTracks = int(group.size());
//                if (nTracks>0) {
//                  TrackExtended * siTrack = group[0];
//                  // preserving positiveness of matrix
//                  float eD0 = siTrack->getCovMatrix()[0];
//                  float eZ0 = siTrack->getCovMatrix()[9];
//                  float sD0 = sqrt(eD0/epar[0]);
//                  float sZ0 = sqrt(eZ0/epar[9]);
//                  epar[0]  = eD0;
//                  epar[9]  = eZ0;
//                  epar[1]  = sD0*epar[1];
//                  epar[3]  = sD0*epar[3];
//                  epar[10] = sD0*epar[10];
//                  epar[7]  = sZ0*epar[7];
//                  epar[8]  = sZ0*epar[8];
//                  epar[13] = sZ0*epar[13];
//                  epar[6]  = sD0*sZ0*epar[6];
//                  trkExt->setZ0(siTrack->getZ0());
//                  trkExt->setD0(siTrack->getD0());
//                  trkExt->setCovMatrix(epar);
//                }
//              }
//            }
//          }
//          delete[] x_h;
//          delete[] y_h;
//          delete[] z_h;
//          delete[] rR_h;
//          delete[] rZ_h;
//          delete[] idet_h;
//          delete[] ityp_h;
//          delete[] lhits;
//        }
//      }
//    }
//    
//    for (int iP=0;iP<nPairs;++iP) {
//      TrackHitPair * trkHitPair = pairs[iP];
//      delete trkHitPair;
//    }
//    
//    pairs.clear();
//    
//  }
//  
//  
//}

HelixClass * FullLDCTracking_MarlinTrk::GetExtrapolationHelix( TrackExtended * track) {
  
  TrackerHitExtendedVec hitVec = track->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());
  
  HelixClass * helix = new HelixClass();
  
  int nHitsFit = nHits;
  if (nHits>_nHitsExtrapolation)
    nHitsFit = _nHitsExtrapolation;
  
  float * ampl = new float[nHitsFit];
  float * xhit = new float[nHitsFit];
  float * yhit = new float[nHitsFit];
  float * zhit = new float[nHitsFit];
  
  
  if (nHits<=_nHitsExtrapolation) {
    for (int iH=0;iH<nHitsFit;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }    
    
  }
  else {
    int * index = new int[nHits];
    float * zcoor = new float[nHits];        
    for (int iH=0;iH<nHits;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      zcoor[iH] = fabs(float(trkHit->getPosition()[2]));
    }
    int order = 1;
    GeneralSorting(index,zcoor,order,nHits);
    for (int iH=0;iH<nHitsFit;++iH) {
      int idx = index[iH];
      TrackerHit * trkHit = hitVec[idx]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }
    delete[] index;
    delete[] zcoor;
  }
  
  ClusterShapes * shapes = new ClusterShapes(nHitsFit,ampl,xhit,yhit,zhit);
  
  float parSh[5];
  float dparSh[5];
  float distmax = 0;
  float chi2Sh = 0;
  // do fitting
  shapes->FitHelix(500, 0, 1, parSh, dparSh, chi2Sh, distmax);  
  
  float x0Sh = parSh[0];
  float y0Sh = parSh[1];
  float r0Sh = parSh[2];
  float bzSh = parSh[3];
  float phi0Sh = parSh[4];
  float signPz = 1;
  float zBegin = zhit[0];
  float zEnd   = zhit[0];
  float zMax   = fabs(zBegin);
  float zMin   = fabs(zEnd);
  for (int iH=1;iH<nHitsFit;++iH) {
    float zCurrent = fabs(zhit[iH]);
    if (zCurrent>zMax) {
      zMax = zCurrent;
      zEnd = zhit[iH];
    }
    if (zCurrent<zMin) {
      zMin = zCurrent;
      zBegin = zhit[iH];
    }      
  }
  
  delete shapes;
  
  if (zEnd<zBegin)
    signPz = -1;
  helix->Initialize_BZ(x0Sh, y0Sh, r0Sh, 
                       bzSh, phi0Sh, _bField,signPz,
                       zBegin);
  
  delete[] ampl;
  delete[] xhit;
  delete[] yhit;
  delete[] zhit;
  
  return helix;
  
}


void FullLDCTracking_MarlinTrk::AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec,
                                                      float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      TrackerHitExtendedVec hitVecGrp = trkGrp->getTrackerHitExtendedVec();
      int nHits_Grp = int(hitVecGrp.size());
      float zMin = 1.0e+20;
      float zMax = -1.0e+20;
      float startPoint[3] = {0.,0.,0.};
      float endPoint[3]   = {0.,0.,0.};
      for (int iH=0;iH<nHits_Grp;++iH) {
        TrackerHitExtended * trkHitExt = hitVecGrp[iH];
        float pos[3] = {0.,0.,0.};
        for (int iC=0;iC<3;++iC) 
          pos[iC] = float(trkHitExt->getTrackerHit()->getPosition()[iC]);         
        if (pos[2]>zMax) {
          zMax = pos[2];
          for (int iC=0;iC<3;++iC)
            endPoint[iC] = pos[iC];           
        }
        if (pos[2]<zMin) {
          zMin = pos[2];
          for (int iC=0;iC<3;++iC)
            startPoint[iC] = pos[iC];
        }
      }
      trkGrp->setStart(startPoint);
      trkGrp->setEnd(endPoint);
    }
  }
  
  
  // replace previous version with faster loop ordering
  
  std::vector<float> minDistances(nHits, dcut);
  std::vector<TrackExtended*> tracksToAttach(nHits);
  std::vector< std::vector<float> > HitPositions(nHits);
  std::vector<int> HitSign(nHits);//Positive or Negative side
  for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
    tracksToAttach[iH]=NULL;
    //Get all TrackerHit positions, so we only have to get them once
    TrackerHit* temphit = hitVec[iH]->getTrackerHit();
    const double *temppos = temphit->getPosition();
    HitPositions[iH].push_back(float(temppos[0]));
    HitPositions[iH].push_back(float(temppos[1]));
    HitPositions[iH].push_back(float(temppos[2]));
    HitSign[iH]=std::signbit(temppos[2]);
  }    
  
  streamlog_out(DEBUG4) << " Starting loop " << nTrk << " tracks   and  " << nHits << " hits" << std::endl;
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    int tanlambdaSign = std::signbit(foundTrack->getTanLambda());//we only care about positive or negative
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      float tanLambda = trkGrp->getTanLambda();
      float omega = trkGrp->getOmega();
      float d0 = trkGrp->getD0();
      float z0 = trkGrp->getZ0();
      float phi0 = trkGrp->getPhi();
      float startPointZ = trkGrp->getStart()[2];
      float endPointZ   = trkGrp->getEnd()[2];
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      float OnePFivehalfPeriodZ = 1.5*fabs(acos(-1.)*tanLambda/omega);
      
      for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
        
        //check if the hit and the track or on the same side
        //xor return 1, if hits are different
        if ( tanlambdaSign^HitSign[iH] ) continue;
        
        float DeltaStart = fabs(HitPositions[iH][2]-startPointZ);
        float DeltaEnd = fabs(HitPositions[iH][2]-endPointZ);
        bool consider = DeltaStart <= OnePFivehalfPeriodZ;
        consider = consider || (DeltaEnd <= OnePFivehalfPeriodZ);
        consider = consider || ( (HitPositions[iH][2]>=startPointZ) && (HitPositions[iH][2]<=endPointZ) );
        
        if(consider){
          float distance = helix.getDistanceToPoint(HitPositions[iH], minDistances[iH]);
          if (distance < minDistances[iH]) {
            minDistances[iH] = distance;
            tracksToAttach[iH] = foundTrack;
          }
        }
      } // loop over leftover TPC hits
    } //groups in tracks
  } // loop over all tracks
  
  for (int iH=0;iH<nHits;++iH) {
    TrackerHitExtended * trkHitExt = hitVec[iH];
    if (tracksToAttach[iH]!=NULL) {
      tracksToAttach[iH]->addTrackerHitExtended(trkHitExt);
      trkHitExt->setTrackExtended( tracksToAttach[iH] );
      trkHitExt->setUsedInFit( false );
    }
  }
  
  streamlog_out(DEBUG4) << " Fast loop done " << std::endl;
  
  
  //     for (int iH=0;iH<nHits;iH++) { // loop over leftover TPC hits
  //    TrackerHitExtended * hitExt = hitVec[iH];
  //    float pos[3];
  //    for (int ip=0;ip<3;++ip) 
  //        pos[ip] = float(hitExt->getTrackerHit()->getPosition()[ip]);
  //    float minDist = 1.0e+20;
  //    TrackExtended * trackToAttach = NULL;
  //    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
  //        TrackExtended * foundTrack = _trkImplVec[iT];
  //        float tanLambdaFound = foundTrack->getTanLambda();
  //        float product = tanLambdaFound*pos[2];
  //        if (product>0) {
  //          GroupTracks * group = foundTrack->getGroupTracks();
  //          TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
  //          int nTrkGrp = int(tracksInGroup.size());
  //          for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
  //            TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
  //            float tanLambda = trkGrp->getTanLambda();
  //            float omega = trkGrp->getOmega();
  //            float d0 = trkGrp->getD0();
  //            float z0 = trkGrp->getZ0();
  //            float phi0 = trkGrp->getPhi();
  //            float dist[3];
  //            float startPoint[3];
  //            float endPoint[3];
  //            for (int iC=0;iC<3;++iC) {
  //              startPoint[iC] = trkGrp->getStart()[iC];
  //              endPoint[iC] = trkGrp->getEnd()[iC];
  //            }
  //            HelixClass helix;
  //            helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
  //            float halfPeriodZ = fabs(acos(-1.)*tanLambda/omega);
  //            helix.getDistanceToPoint(pos,dist);
  //            float DeltaStart = fabs(pos[2]-startPoint[2]);
  //            float DeltaEnd = fabs(pos[2]-endPoint[2]);
  //            bool consider = DeltaStart <= 1.5*halfPeriodZ;
  //            consider = consider || (DeltaEnd <= 1.5*halfPeriodZ);
  //            consider = consider || ( (pos[2]>=startPoint[2]) && (pos[2]<=endPoint[2]) );
  // //                 float ZMin = DeltaStart;
  // //                 if (DeltaEnd<ZMin)
  // //                   ZMin = DeltaEnd;
  //            if (dist[2]<dcut && consider && dist[2]<minDist) {
  //              minDist = dist[2];
  //              trackToAttach = foundTrack;
  //            }
  //          }
  //        }
  //    }
  
  //    if (trackToAttach!=NULL) {
  //   trackToAttach->addTrackerHitExtended(hitExt);
  //  hitExt->setTrackExtended( trackToAttach );
  //  hitExt->setUsedInFit( false );
  //  if(trackToAttach!=tracksToAttach[iH])std::cout << " Check Failed" << trackToAttach << "  " << tracksToAttach[iH] << std::endl;
  //
  //}
  //else {
  ///     std::cout << iH << " hit is not assigned : distance to closest track = " << minDist << std::endl;
  ///}
  //}
  //std::cout << " Slow loop done " << std::endl;
  
}

void FullLDCTracking_MarlinTrk::AssignSiHitsToTracks(TrackerHitExtendedVec hitVec,
                                                     float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_allNonCombinedTPCTracks.size());
  
  streamlog_out(DEBUG4) << "AssignSiHitsToTracks : Number of hits to assign " <<  hitVec.size() << " : Number of available tracks = " << nTrk << std::endl;
  
  std::map <TrackExtended*,bool> flagTrack;
  std::map <TrackerHitExtended*,bool> flagHit;
  TrackHitPairVec pairs;
  flagTrack.clear();
  flagHit.clear();
  pairs.clear();
  
  for (int iH=0;iH<nHits;++iH) {
    
    float pos[3];
    TrackerHitExtended * trkHitExt = hitVec[iH];
    TrackerHit * hit = trkHitExt->getTrackerHit();
    
    for (int ip=0;ip<3;++ip) {
      pos[ip] = float(hit->getPosition()[ip]);
    }
    
    for (int iT=0;iT<nTrk;++iT) {
      
      TrackExtended * trkExt = _allNonCombinedTPCTracks[iT];
      
      float tanLambda = trkExt->getTanLambda();       
      float product = pos[2]*tanLambda;
      
      streamlog_out(DEBUG2) << "AssignSiHitsToTracks : product =  " << product << " z hit = " << pos[2] <<  std::endl;
      
      
      if (product>0) {
        
        float d0 = trkExt->getD0();
        float z0 = trkExt->getZ0();
        float phi0 = trkExt->getPhi();
        float omega = trkExt->getOmega();
        tanLambda = trkExt->getTanLambda();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
        float distance = helix.getDistanceToPoint(pos,dcut);
        
        streamlog_out(DEBUG2) << "AssignSiHitsToTracks : distance =  " << distance << " cut = " << dcut << std::endl;
        
        if (distance<dcut) {
          TrackHitPair * trkHitPair = 
          new TrackHitPair(trkExt,trkHitExt,distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;
        }
      }
    }
  }
  
  int nPairs = int(pairs.size());
  
  streamlog_out(DEBUG4) << "AssignSiHitsToTracks : Number of track hit pairs to try =  " << nPairs << std::endl;
  
  if (nPairs>0) {
    
    SortingTrackHitPairs(pairs);
    
    for (int iP=0;iP<nPairs;++iP) {
      
      TrackHitPair * trkHitPair = pairs[iP];
      TrackExtended * trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended * trkHitExt = 
      
      trkHitPair->getTrackerHitExtended();
      
      if (flagTrack[trkExt] && flagHit[trkHitExt]) {              
        
        // get the hits already assigned to the track
        TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();
        
        int nTotH = int(hitsInTrack.size());
        int nHitsInFit = 0;
        
        for (int iTH=0;iTH<nTotH;++iTH) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
          
          // count the number of hits used in the fit
          if (hitInTrack->getUsedInFit()) {
            nHitsInFit++; 
          }
        }
        
        int iHitInFit = 0;
        
        // add the previously used hits from the track to the vectors 
        
        EVENT::TrackerHitVec trkHits;
        
        for (int iHit=0;iHit<nTotH;++iHit) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
          if (hitInTrack->getUsedInFit()) {
            TrackerHit * hit = hitInTrack->getTrackerHit();
            iHitInFit++;
            if(hit) { 
              trkHits.push_back(hit);   
            }
            else{
              throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: TrackerHit pointer == NULL ")  ) ;
            }
          }
        }
        
        // add the hit to be attached to the vectors 
        TrackerHit * remainHit = trkHitExt->getTrackerHit();
        iHitInFit++;
        trkHits.push_back(remainHit);
        
        
        double chi2_D;
        int ndf;
        
        if( trkHits.size() < 3 ) return ;
        
        // sort the hits in R
        std::vector< std::pair<float, EVENT::TrackerHit*> > r2_values;
        r2_values.reserve(trkHits.size());
        
        for (TrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
          EVENT::TrackerHit* h = *it;
          float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
          r2_values.push_back(std::make_pair(r2, *it));
        }
        
        sort(r2_values.begin(),r2_values.end());
        
        trkHits.clear();
        trkHits.reserve(r2_values.size());
        
        for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
          trkHits.push_back(it->second);
        }
        
        
        EVENT::TrackerHitVec::iterator it = trkHits.begin();
        
        streamlog_out(DEBUG2) << "AssignSiHitsToTracks: Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
                
        
        MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
        
        IMPL::TrackStateImpl pre_fit ;
        
        
        /** Provides the values of a track state from the first, middle and last hits in the hit_list. */
        int error = createPrefit( trkHits, &pre_fit, _bField, IMarlinTrack::backward);
        
        if ( error != IMarlinTrack::success ) {
          
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: creation of prefit fails with error " << error << std::endl;
          
          delete marlin_trk ;
          return;
          
        }
        
        // setup initial dummy covariance matrix
        EVENT::FloatVec covMatrix;
        covMatrix.resize(15);
        
        for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
          covMatrix[icov] = 0;
        }
        
        covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
        covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
        covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
        covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
        covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
        
        pre_fit.setCovMatrix(covMatrix);
        
        error = MarlinTrk::createFit( trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward , _maxChi2PerHit );
        
        if ( error != IMarlinTrack::success ) {
          
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: creation of fit fails with error " << error << std::endl;
          
          delete marlin_trk ;
          return ;
          
        }
        
        
        const gear::Vector3D point(0.,0.,0.); // nominal IP
        int return_code = 0;
        
        TrackStateImpl trkState ;
        return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf ) ;
        
        delete marlin_trk ;
        
        if ( error != IMarlinTrack::success ) {
          
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: propagate to IP fails with error " << error << std::endl;
          
          return ;
          
        }
        
        if ( ndf < 0  ) {
          
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: Fit failed : NDF is less that zero  " << ndf << std::endl;
          
          return ;
          
        }
        
        
        float chi2Fit = chi2_D/float(ndf);
        
        if ( chi2Fit > _chi2FitCut ) {
          
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: track fail Chi2 cut of " << _chi2FitCut << " chi2 of track = " <<  chi2Fit << std::endl;
          
          return ;
          
        }
              
        
        
        // note trackAR which is of type TrackExtended, only takes fits with ref point = 0,0,0 
        trkExt->setOmega(trkState.getOmega());
        trkExt->setTanLambda(trkState.getTanLambda());
        trkExt->setPhi(trkState.getPhi());
        trkExt->setD0(trkState.getD0());
        trkExt->setZ0(trkState.getZ0());
        
        
        float cov[15];
        
        for (int i = 0 ; i<15 ; ++i) {
          cov[i] = trkState.getCovMatrix().operator[](i);
        }
        
        trkExt->setCovMatrix(cov);
        trkExt->setChi2(chi2_D);
        trkExt->setNDF(ndf);
        
        trkExt->addTrackerHitExtended( trkHitExt );
        trkHitExt->setTrackExtended( trkExt );
        trkHitExt->setUsedInFit( true );
        flagTrack[trkExt] = false;
        flagHit[trkHitExt] = false;
        
        
        streamlog_out(DEBUG2) << "AssignSiHitsToTracks: Start Fitting: Hit " << trkHitExt << " successfully assigned to track " << trkExt << std::endl;
        
      }
    }
    
    for (int iP=0;iP<nPairs;++iP) {
      TrackHitPair * trkHitPair = pairs[iP];
      delete trkHitPair;
    }
    
    pairs.clear();
    
  }
}

void FullLDCTracking_MarlinTrk::PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, int iopt) {
  // iopt = 1 false Si and TPC merging
  // iopt = 2 false Si and TPC forced merging
  // iopt = 3 false TPC segments merging
  // iopt = 4 false Comb Si and TPC merging
  // iopt = 5 false Comb TPC and TPC merging
  // iopt = 6 unmerged TPC and Si segments (soft merging)
  // iopt = 7 unmerged TPC and Si segments (forced merging)
  // iopt = 8 unmerged Comb and TPC
  // iopt = 9 unmerged TPC segments
  
  char strg[200];
  
  try {
    
    Track * firstTrack = firstTrackExt->getTrack();
    Track * secondTrack = secondTrackExt->getTrack();
    
    std::string firstColName = _TPCTrackMCPCollName;
    std::string secondColName = _TPCTrackMCPCollName;
    
    if (iopt==1) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==2) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==3) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==4) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==5) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==6) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==7) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==8) {
      secondColName = _SiTrackMCPCollName;
    }    
    else {
      secondColName = _TPCTrackMCPCollName;
    }
    
    
    LCCollection * firstCol = _evt->getCollection(firstColName.c_str());
    LCCollection * secondCol = _evt->getCollection(secondColName.c_str());
    
    
    LCRelationNavigator firstNav(firstCol);
    LCRelationNavigator secondNav(secondCol);
    LCObjectVec firstVec = firstNav.getRelatedToObjects(firstTrack);
    FloatVec firstWeights = firstNav.getRelatedToWeights(firstTrack);
    LCObject * firstMCP = NULL;
    float firstWght = 0;
    int nObj = firstVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (firstWeights[iObj]>firstWght) {
        firstWght = firstWeights[iObj];
        firstMCP = firstVec[iObj];
      }
    }
    
    LCObjectVec secondVec = secondNav.getRelatedToObjects(secondTrack);
    FloatVec secondWeights = secondNav.getRelatedToWeights(secondTrack);
    LCObject * secondMCP = NULL;
    float secondWght = 0;
    nObj = secondVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (secondWeights[iObj]>secondWght) {
        secondWght = secondWeights[iObj];
        secondMCP = secondVec[iObj];
      }
    }
    
    
    float d0First = firstTrackExt->getD0();
    float z0First = firstTrackExt->getZ0();
    float omegaFirst = firstTrackExt->getOmega();
    float tanLFirst = firstTrackExt->getTanLambda();
    float phi0First = firstTrackExt->getPhi();
    
    float d0Second = secondTrackExt->getD0();
    float z0Second = secondTrackExt->getZ0();
    float omegaSecond = secondTrackExt->getOmega();
    float tanLSecond = secondTrackExt->getTanLambda();
    float phi0Second = secondTrackExt->getPhi();            
    
    HelixClass firstHelix;
    firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float pxFirst = firstHelix.getMomentum()[0];
    float pyFirst = firstHelix.getMomentum()[1];
    float pzFirst = firstHelix.getMomentum()[2];            
    
    HelixClass secondHelix;
    secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
    float pxSecond = secondHelix.getMomentum()[0];
    float pySecond = secondHelix.getMomentum()[1];
    float pzSecond = secondHelix.getMomentum()[2];          
    
    float dPx = pxFirst - pxSecond;
    float dPy = pyFirst - pySecond;
    float dPz = pzFirst - pzSecond;
    
    float dPplus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    dPx = pxFirst + pxSecond;
    dPy = pyFirst + pySecond;
    dPz = pzFirst + pzSecond;
    
    float dPminus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
    float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
    
    //SJA:FIXME Hardcoded cut here should be removed 
    if(pFirst < 1.5 && pSecond < 1.5 )return; 
    
    const float sigmaPOverPFirst  = sqrt(firstTrackExt->getCovMatrix()[5])/fabs(omegaFirst);
    const float sigmaPOverPSecond = sqrt(secondTrackExt->getCovMatrix()[5])/fabs(omegaSecond);
    //    const float deltaP = fabs(pFirst-pSecond);
    const float sigmaPFirst = pFirst*sigmaPOverPFirst;
    const float sigmaPSecond = pSecond*sigmaPOverPSecond;
    //const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
    //    const float significance = deltaP/sigmaDeltaP;
    
    
    float den = pFirst;
    if (pSecond<pFirst)
      den = pSecond;
    
    dPplus  = dPplus/den;
    dPminus = dPminus/den; 
    
    if (firstMCP!=secondMCP && iopt < 6) {
      
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Erroneous combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Erroneous merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Erroneous merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Erroneous merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      
      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }
    else if (firstMCP==secondMCP && ( (iopt==8) || (iopt==9) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      if (iopt==8) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Comb segments (iopt=8) --->" << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC segments (iopt=9) --->" << std::endl;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << std::endl;
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, true);
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }    
    else if (firstMCP==secondMCP && ( (iopt == 6) || (iopt == 7) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      
      if (iopt ==6) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=6) --->" << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=7) --->" << std::endl;
        deltaOmega = _dOmegaForForcedMerging;
        deltaAngle = _angleForForcedMerging;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << std::endl;
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      streamlog_out(DEBUG4) << std::endl;      
      
    }else if (firstMCP==secondMCP && iopt < 6) {
      //      return;
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Correct combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Correct merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Correct merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Correct merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Correct merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      streamlog_out(DEBUG3) << "    p         error       D0      Z0     Px      Py      Pz      wieght" << std::endl;
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt,true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      delete combinedTrack->getGroupTracks();
      delete combinedTrack;
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }
  }
  
  catch(DataNotAvailableException &e){};
  
}    

void FullLDCTracking_MarlinTrk::GeneralSorting(int * index, float * val, int direct, int nVal) {
  /**
   Sorting of index vector in ascending (0) /descending (!=0) order of val
   */
  
  float valOne, valTwo, valTemp;
  int   indTemp;
  for (int i=0; i<nVal; ++i) {
    index[i] = i;
  }
  
  for (int i = 0 ; i < nVal-1; i++) {
    for (int j = 0; j < nVal-i-1; j++) {      
      valOne = val[j];
      valTwo = val[j+1];
      bool order = valOne > valTwo;
      if (direct>0) 
        order = valOne <= valTwo;
      if( order )
      {
        valTemp = val[j];
        val[j] = val[j+1];
        val[j+1] = valTemp;
        indTemp = index[j];
        index[j] = index[j+1];
        index[j+1] = indTemp;
      }
    }  
  }
  
  
}

int FullLDCTracking_MarlinTrk::SegmentRadialOverlap(TrackExtended* first, TrackExtended* second){
  
  
  int nTrkGrpFirst = 0;
  int nTrkGrpSecond = 0;
  TrackerHitVec hitvecFirst;
  TrackerHitVec hitvecSecond;
  GroupTracks * groupFirst = first->getGroupTracks();
  GroupTracks * groupSecond = second->getGroupTracks();
  
  if(groupFirst!=NULL){
    
    TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
    nTrkGrpFirst = int(tracksInGroupFirst.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
      
      TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i =0; i<hitVec.size(); ++i){
        hitvecFirst.push_back(hitVec[i]->getTrackerHit());        
      }
    }
  }
  if(groupSecond!=NULL){
    
    TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
    nTrkGrpSecond = int(tracksInGroupSecond.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
      TrackerHitExtendedVec hitVec = 
      trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i=0;i<hitVec.size();++i){
        hitvecSecond.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }
  
  
  int nhitsFirst = (int)hitvecFirst.size();
  int nhitsSecond = (int)hitvecSecond.size();
  int count = 0;
  for(int i =0;i<nhitsFirst;++i){
    float xi = (float)hitvecFirst[i]->getPosition()[0];
    float yi = (float)hitvecFirst[i]->getPosition()[1];
    float ri = sqrt(xi*xi+yi*yi);
    if(ri < _tpc_inner_r || ri > _tpc_pad_height)continue;
    for(int j =0;j<nhitsSecond;++j){
      float xj = (float)hitvecSecond[j]->getPosition()[0];
      float yj = (float)hitvecSecond[j]->getPosition()[1];
      float rj = sqrt(xj*xj+yj*yj);
      if(fabs(ri-rj)<_tpc_pad_height/2.0)count++;
    }
  }  
  return count;
}

bool FullLDCTracking_MarlinTrk::VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt){
  
  
  const float d0First = firstTrackExt->getD0();
  const float z0First = firstTrackExt->getZ0();
  const float omegaFirst = firstTrackExt->getOmega();
  const float tanLFirst = firstTrackExt->getTanLambda();
  const float phi0First = firstTrackExt->getPhi();
  
  const float d0Second = secondTrackExt->getD0();
  const float z0Second = secondTrackExt->getZ0();
  const float omegaSecond = secondTrackExt->getOmega();
  const float tanLSecond = secondTrackExt->getTanLambda();
  const float phi0Second = secondTrackExt->getPhi();        
  
  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
  const float pxFirst = firstHelix.getMomentum()[0];
  const float pyFirst = firstHelix.getMomentum()[1];
  const float pzFirst = firstHelix.getMomentum()[2];        
  
  HelixClass secondHelix;
  secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  const float pxSecond = secondHelix.getMomentum()[0];
  const float pySecond = secondHelix.getMomentum()[1];
  const float pzSecond = secondHelix.getMomentum()[2];      
  const float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
  const float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
  
  //SJA:FIXME hardcoded cut 
  if(pFirst<2.5 || pSecond<2.5)return false;
  
  TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt,true);
  bool veto = false;
  if(combinedTrack!=NULL){
    //SJA:FIXME hardcoded cut 
    if(combinedTrack->getNDF()+15<firstTrackExt->getNDF()+secondTrackExt->getNDF()+5)veto=true;
    delete combinedTrack->getGroupTracks();
    delete combinedTrack;
  }else{
    veto = true;
  }
  if(SegmentRadialOverlap(firstTrackExt,secondTrackExt)>10)veto=true;
  return veto;
  
}


void FullLDCTracking_MarlinTrk::check(LCEvent * evt) { }

void FullLDCTracking_MarlinTrk::end() { 
  
  delete _encoder ;
  
}

void FullLDCTracking_MarlinTrk::setupGearGeom( const gear::GearMgr* gearMgr ){
  
  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  //-- VXD Parameters--
  _nLayersVTX = 0 ;
  const gear::VXDParameters* pVXDDetMain = 0;
  const gear::VXDLayerLayout* pVXDLayerLayout = 0;
  
  try{
    
    streamlog_out( DEBUG9 ) << " filling VXD parameters from gear::SITParameters " << std::endl ;
    
    pVXDDetMain = &Global::GEAR->getVXDParameters();
    pVXDLayerLayout = &(pVXDDetMain->getVXDLayerLayout());
    _nLayersVTX = pVXDLayerLayout->getNLayers();
  }
  catch( gear::UnknownParameterException& e){
    
    streamlog_out( DEBUG9 ) << " ### gear::VXDParameters Not Present in GEAR FILE" << std::endl ;
    
  }
  
  
  
  //-- SIT Parameters--
  _nLayersSIT = 0 ;
  const gear::ZPlanarParameters* pSITDetMain = 0;
  const gear::ZPlanarLayerLayout* pSITLayerLayout = 0;
  
  try{
    
    streamlog_out( DEBUG9 ) << " filling SIT parameters from gear::SITParameters " << std::endl ;
    
    pSITDetMain = &Global::GEAR->getSITParameters();
    pSITLayerLayout = &(pSITDetMain->getZPlanarLayerLayout());
    _nLayersSIT = pSITLayerLayout->getNLayers();
    
  }
  catch( gear::UnknownParameterException& e){
    
    streamlog_out( DEBUG9 ) << " ### gear::SITParameters Not Present in GEAR FILE" << std::endl ;
    
  }
  
  if( _nLayersSIT == 0 ){
    // try the old LOI style key value pairs as defined in the SSit03 Mokka drive
    try{
      
      streamlog_out( MESSAGE ) << "  FullLDCTracking_MarlinTrk - Simple Cylinder Based SIT using parameters defined by SSit03 Mokka driver " << std::endl ;
      
      // SIT
      
      const gear::GearParameters& pSIT = gearMgr->getGearParameters("SIT");
      
      const EVENT::DoubleVec& SIT_r   =  pSIT.getDoubleVals("SITLayerRadius" )  ;
      const EVENT::DoubleVec& SIT_hl  =  pSIT.getDoubleVals("SITSupportLayerHalfLength" )  ;
      
      _nLayersSIT = SIT_r.size() ; 
      
      if (_nLayersSIT != SIT_r.size() || _nLayersSIT != SIT_hl.size()) {
        
        streamlog_out( ERROR ) << "ILDSITCylinderKalDetector miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << std::endl ;
        exit(1);
        
      }
    }
    catch( gear::UnknownParameterException& e){
      
      streamlog_out( DEBUG9 ) << " ### gear::SIT Parameters from as defined in SSit03 Not Present in GEAR FILE" << std::endl ;
      
    } 
    
  }
  
  
  
  //-- FTD Parameters--
  _petalBasedFTDWithOverlaps = false;  
  _nLayersFTD = 0;
  
  try{
    
    streamlog_out( DEBUG9 ) << " filling FTD parameters from gear::FTDParameters " << std::endl ;
    
    const gear::FTDParameters&   pFTD      = Global::GEAR->getFTDParameters();
    const gear::FTDLayerLayout&  ftdlayers = pFTD.getFTDLayerLayout() ;
    
    _nLayersFTD = ftdlayers.getNLayers() ;
    
    for (unsigned int disk=0; disk < _nLayersFTD; ++disk) {
      
      _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 0, 1) ); // front petal even numbered
      
      if ( ftdlayers.getNPetals(disk) > 0) {
        _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 1, 1) );  // front petal odd numbered
        _petalBasedFTDWithOverlaps = true;
      }
      
    }
    
    // SJA: Here we increase the size of _nlayersFTD as we are treating the 
    _nLayersFTD =_zLayerFTD.size() ;     
    
  }
  catch( gear::UnknownParameterException& e){
    
    streamlog_out( DEBUG9 ) << " ### gear::FTDParameters Not Present in GEAR FILE" << std::endl ;
    
  } 
  
  if( _nLayersFTD == 0 ){
    
    // FTD
    try{
      
      streamlog_out( MESSAGE ) << "  FullLDCTracking_MarlinTrk - Simple Disc Based FTD using parameters defined by SFtd05 Mokka driver " << std::endl ;
      
      const gear::GearParameters& pFTD = gearMgr->getGearParameters("FTD");
      
      const EVENT::DoubleVec* pFTD_z   = NULL;
      
      streamlog_out( MESSAGE ) << " For FTD using parameters defined by SFtd05 Mokka driver " << std::endl ;
      
      pFTD_z = &pFTD.getDoubleVals("FTDZCoordinate" )  ;
      
      _nLayersFTD = pFTD_z->size();
      
      for (unsigned int i = 0; i<_nLayersFTD; ++i) {
        _zLayerFTD.push_back((*pFTD_z)[i]);
      }
    }
    catch( gear::UnknownParameterException& e){
      
      streamlog_out( DEBUG9 ) << " ### gear::FTD Parameters as defined in SFtd05 Not Present in GEAR FILE" << std::endl ;
      
    } 
  }
  
}
