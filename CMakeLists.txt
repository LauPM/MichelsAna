art_make( MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${ROOT_BASIC_LIB_LIST}
  
  cetlib cetlib_except
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas

  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardata_RecoObjects
  lardata_ArtDataHelper
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_MCBase
  lardataobj_AnalysisBase
  lardataobj_Simulation
  nusimdata_SimulationBase
  nug4_ParticleNavigation

  larsim_Simulation nug4_ParticleNavigation
  larsim_Utils
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service

  larreco_RecoAlg
  larpandora_LArPandoraInterface
	dunereco_AnaUtils
)

install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)