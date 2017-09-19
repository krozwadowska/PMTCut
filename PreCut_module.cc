////////////////////////////////////////////////////////////////////////
// Class:       PreCut
// Plugin Type: analyzer (art v2_06_03)
// File:        PreCut_module.cc
//
// Generated at Fri Sep  8 10:27:23 2017 by Karolina M. Rozwadowska using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

// Art include
  #include "art/Framework/Core/EDAnalyzer.h"
  #include "art/Framework/Core/ModuleMacros.h"
  #include "art/Framework/Principal/Event.h"
  #include "art/Framework/Principal/Handle.h"
  #include "art/Framework/Principal/Run.h"
  #include "art/Framework/Principal/SubRun.h"
  #include "canvas/Utilities/InputTag.h"
  #include "fhiclcpp/ParameterSet.h"
  #include "messagefacility/MessageLogger/MessageLogger.h"
  #include "art/Framework/Services/Optional/TFileService.h"
  #include "art/Framework/Services/Optional/TFileDirectory.h"
  #include "canvas/Persistency/Common/FindManyP.h"
  

  // Data products include
  #include "lardataobj/MCBase/MCTrack.h"
  #include "lardataobj/RecoBase/PFParticle.h"
  #include "lardataobj/RecoBase/OpHit.h"
  #include "lardataobj/RecoBase/OpFlash.h"
  #include "lardataobj/RecoBase/Track.h"
  #include "lardataobj/RecoBase/Cluster.h"
  #include "lardataobj/RecoBase/Hit.h"
  #include "lardataobj/RecoBase/SpacePoint.h"
  #include "lardataobj/AnalysisBase/T0.h"
  

  // LArSoft include
  #include "uboone/UBFlashFinder/PECalib.h"
  #include "larsim/MCCheater/BackTracker.h"
  #include "lardata/Utilities/AssociationUtil.h"
  #include "larcore/Geometry/Geometry.h"
  #include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
  #include "lardataobj/AnalysisBase/T0.h"
  #include "lardataobj/MCBase/MCDataHolder.h"
  #include "lardataobj/MCBase/MCHitCollection.h"
  
  // Root include
  #include "TString.h"
  #include "TTree.h"
  #include "TH1D.h"

class PreCut;


class PreCut : public art::EDAnalyzer {
public:
  explicit PreCut(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PreCut(PreCut const &) = delete;
  PreCut(PreCut &&) = delete;
  PreCut & operator = (PreCut const &) = delete;
  PreCut & operator = (PreCut &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

  TTree* _tree1;
  int _run, _subrun, _event;
  
  int _fv, _ccnc, _nupdg;
  
  std::vector<double> _ophit_beam_pe, _ophit_beam_time, _ophit_beam_timeabs, _ophit_beam_tick;
  int _n_ophit_beam;

  
  int _nbeamfls;
  std::vector<double> _beamfls_time, _beamfls_pe, _beamfls_z;
  
  bool _oh9_passed40pe;
  bool _oh9_passed50pe;
  bool _fls_passed40pe;
  bool _fls_passed50pe;
  bool _numu_cc_fv;
};


PreCut::PreCut(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
{
  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree", "");
  _tree1->Branch("run",					&_run,				"run/I");
  _tree1->Branch("subrun",				&_subrun,			"subrun/I");
  _tree1->Branch("event",				&_event,			"event/I");
  
  _tree1->Branch("fv",                   		&_fv,                    	"fv/I");
  _tree1->Branch("ccnc",                 		&_ccnc,                  	"ccnc/I");
  _tree1->Branch("nupdg",                		&_nupdg,                 	"nupdg/I");  
  
  _tree1->Branch("n_ophit_beam", 			&_n_ophit_beam, 		"n_ophit_beam/I");
  _tree1->Branch("ophit_beam_pe",			"std::vector<double>", 		&_ophit_beam_pe);
  _tree1->Branch("ophit_beam_time",			"std::vector<double>", 		&_ophit_beam_time);
  _tree1->Branch("ophit_beam_timeabs",			"std::vector<double>", 		&_ophit_beam_timeabs);
  _tree1->Branch("ophit_beam_tick",			"std::vector<double>", 		&_ophit_beam_tick);
  
  _tree1->Branch("nbeamfls",                   		&_nbeamfls,          		"nbeamfls/I");
  _tree1->Branch("beamfls_time",               		"std::vector<double>", 		&_beamfls_time);
  _tree1->Branch("beamfls_pe",                 		"std::vector<double>",          &_beamfls_pe);
  _tree1->Branch("beamfls_z",                  		"std::vector<double>",          &_beamfls_z);

  
  _tree1->Branch("oh9_passed40pe", 			&_oh9_passed40pe, 		"oh9_passed40pe/O");
  _tree1->Branch("oh9_passed50pe", 			&_oh9_passed50pe, 		"oh9_passed50pe/O");
  _tree1->Branch("fls_passed40pe", 			&_fls_passed40pe, 		"fls_passed40pe/O");
  _tree1->Branch("fls_passed50pe", 			&_fls_passed50pe, 		"fls_passed50pe/O");

  _tree1->Branch("numu_cc_fv", 				&_numu_cc_fv, 			"numu_cc_fv/O");
  
}

void PreCut::analyze(art::Event const & evt)
{
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();
  
    double this_bin9_pe = 0;
    
    _oh9_passed40pe = false;
    _oh9_passed50pe = false;
    _fls_passed40pe = false;
    _fls_passed50pe = false;
    _numu_cc_fv = false;
    
//-------------------------------------
//     check if in FV, CC, nu mu
//-------------------------------------
      art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      std::vector<art::Ptr<simb::MCTruth> > mclist;
      if (evt.getByLabel("generator",mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
    
      int iList = 0; // 1 nu int per spill
      double nu_vertex_xyz[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                                mclist[iList]->GetNeutrino().Nu().Vy(),
                                mclist[iList]->GetNeutrino().Nu().Vz()};
      _ccnc    = mclist[iList]->GetNeutrino().CCNC();
      _nupdg   = mclist[iList]->GetNeutrino().Nu().PdgCode();
    
    
    double x = nu_vertex_xyz[0];
    double y = nu_vertex_xyz[1];
    double z = nu_vertex_xyz[2];
  
    //This defines our current settings for the fiducial volume
    double FVx = 256.35;
    double FVy = 233;
    double FVz = 1036.8;
    double borderx = 10.;
    double bordery = 20.;
    double borderz = 10.;
  
    _fv=0;
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)){
      _fv=1;
    }
    
//--------------------------------------------------
// 		beam OpHits
//--------------------------------------------------
  art::Handle<std::vector<recob::OpHit>> ophit_beam_h;
  evt.getByLabel("ophitBeam", ophit_beam_h);
  if(!ophit_beam_h.isValid()){
    std::cerr<<"Cannot locate OpHit"<<std::endl;
  }
 _n_ophit_beam = ophit_beam_h->size();
 _ophit_beam_pe.resize(_n_ophit_beam);
 _ophit_beam_time.resize(_n_ophit_beam);
 _ophit_beam_timeabs.resize(_n_ophit_beam);
 _ophit_beam_tick.resize(_n_ophit_beam);
 
   for (size_t oh = 0; oh < ophit_beam_h->size(); oh++){ 
    
    auto const & ophit_beam = (*ophit_beam_h)[oh];
    
   _ophit_beam_pe[oh] 		= ophit_beam.PE();
   _ophit_beam_time[oh] 	= ophit_beam.PeakTime();
   _ophit_beam_timeabs[oh] 	= ophit_beam.PeakTimeAbs();
   _ophit_beam_tick[oh] 	= ophit_beam.PeakTime()/0.015625;
 }	
 
 
//----------------------------------------------------
// 		beam Flashes
//----------------------------------------------------
    ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
    evt.getByLabel("opflashBeam",beamflash_h);
    if( !beamflash_h.isValid() || beamflash_h->empty() ) {
      std::cerr << "Don't have good flashes." << std::endl;
    }
    
    _nbeamfls = beamflash_h->size();
    _beamfls_pe.resize(_nbeamfls);
    _beamfls_time.resize(_nbeamfls);
    _beamfls_z.resize(_nbeamfls);

  
    for (size_t n = 0; n < beamflash_h->size(); n++) {
      auto const& flash_beam = (*beamflash_h)[n];
      _beamfls_pe[n]   = flash_beam.TotalPE();
      _beamfls_time[n] = flash_beam.Time();
      _beamfls_z[n]    = flash_beam.ZCenter();

    }
    
//-------------------------------------------
// 		SELECTION
//-------------------------------------------
 if (_fv==1 && _ccnc==0 && _nupdg==14){
   _numu_cc_fv = true;
 }
   // OpHits 9 tick binning
   for (int tick=205; tick<307; tick=tick+9){
    for (int x=0; x<9; x++){
    for (size_t oh = 0; oh < ophit_beam_h->size(); oh++){ 
      if (_ophit_beam_tick[oh]==(tick+x)){
	this_bin9_pe = this_bin9_pe + _ophit_beam_pe[oh];
      }
    }
    if (this_bin9_pe>=40. && _oh9_passed40pe == false) {
      _oh9_passed40pe = true;
    }
    if (this_bin9_pe>=50. && _oh9_passed50pe == false) {
      _oh9_passed50pe = true;
    }
  }
    this_bin9_pe=0;
 } 

// Flashes
      for (size_t n = 0; n < beamflash_h->size(); n++) {
	if (_beamfls_time[n]>=3.2 && _beamfls_time[n]<=4.8){
	    if (_beamfls_pe[n] >= 40. && _fls_passed40pe == false) {
	      _fls_passed40pe = true;
	    }
	    if (_beamfls_pe[n] >= 50. && _fls_passed50pe == false) {
	      _fls_passed50pe = true;
	    }
	}
    }
  
  _tree1->Fill();
  
} 

DEFINE_ART_MODULE(PreCut)
