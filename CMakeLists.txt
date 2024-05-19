set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-4")
add_subdirectory(fcl)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
art_make(BASENAME_ONLY
	LIB_LIBRARIES
	${ROOT_CORE}
	${ROOT_EVE}
	cetlib cetlib_except
	art_Persistency_Common canvas
	art_Persistency_Provenance canvas
	art_Utilities canvas
	ART_PERSISTENCY_COMMON
	ART_PERSISTENCY_PROVENANCE
	ART_UTILITIES

	larcore::Geometry_Geometry_service
	larcore_Geometry_Geometry_service
	larcorealg::Geometry
	larcorealg_Geometry
	lardata_RecoObjects
	lardata_ArtDataHelper
	lardata_Utilities
	lardataobj_RawData
	lardataobj_RecoBase
	lardataobj_MCBase
	lardataobj_AnalysisBase
	lardataobj_Simulation
	nug4_ParticleNavigation
	larsim_Simulation nug4_ParticleNavigation
	larsim_Utils
	larsim_MCCheater_BackTrackerService_service
	larsim_MCCheater_PhotonBackTrackerService_service
	larsim_MCCheater_ParticleInventoryService_service
	larreco_RecoAlg
	larreco_RecoAlg_PMAlg

	messagefacility::MF_MessageLogger

	SERVICE_LIBRARIES
	duneana_MichelsAna
	${ART_ROOT_IO_ROOTINPUT_SOURCE}

	MODULE_LIBRARIES
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
	messagefacility::MF_MessageLogger
	ROOT::Core
	art_Persistency_Common canvas
	art_Persistency_Provenance canvas
	art_Utilities canvas
	ART_PERSISTENCY_COMMON
	ART_PERSISTENCY_PROVENANCE
	ART_UTILITIES

	larcore_Geometry_Geometry_service
	larcorealg_Geometry
	lardata_RecoObjects
	lardata_ArtDataHelper
	lardata_Utilities
	lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
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
	larreco_RecoAlg_PMAlg
	larpandora_LArPandoraInterface
	dunereco_AnaUtils
	dunecore_DuneObj
	duneprototypes::ProtoDUNEDataUtils
	ProtoDUNEUtilities
	duneana_MichelsAna
	ROOT_BASIC_LIB_LIST
)

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)