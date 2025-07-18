#----------------------------------------------------------------------------
# Set up colors
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(GRSISort)

#----------------------------------------------------------------------------
# Check that we are building out of source
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(FATAL_ERROR "In-source builds are not recommended, please re-run from a separate build directory, or specify the build directory via the -B flag.")
endif()

#----------------------------------------------------------------------------
# find git (to set version correctly)
find_package(Git REQUIRED)

#----------------------------------------------------------------------------
# finde required packages
find_package(X11 REQUIRED)

#----------------------------------------------------------------------------
# find the root package (and add COMPONENTS)
find_package(ROOT REQUIRED COMPONENTS Gui GuiHtml Minuit Spectrum OPTIONAL_COMPONENTS MathMore Proof)

# TODO check if using add_compile_definitions would be better here
if(${ROOT_mathmore_FOUND})
	add_compile_options(-DHAS_MATHMORE)
else()
	message("${Red}Mathmore feature for ROOT not found (ROOT_mathmore_FOUND = ${ROOT_mathmore_FOUND})!${ColourReset}")
endif()
if(${ROOT_xml_FOUND})
	add_compile_options(-DHAS_XML)
else()
	message("${Red}XML feature for ROOT not found (ROOT_xml_FOUND = ${ROOT_xml_FOUND})!${ColourReset}")
endif()

find_package(Python)
if(${Python_VERSION} STREQUAL "2.7")
	add_compile_options(-DHAS_CORRECT_PYTHON_VERSION)
else()
	message(CONFIGURE_LOG "${Yellow}Python version found is ${Python_VERSION} and not 2.7, GUI won't be active${ColourReset}")
endif()

# check the root version (ROOT_VERSION=6.32.02, ROOT_VERSION_COUNT=3, ROOT_VERSION_MAJOR=6, ROOT_VERSION_MINOR=32, ROOT_VERSION_PATCH=2, ROOT_VERSION_TWEAK=0)
if(NOT ${ROOT_VERSION_MAJOR} EQUAL "6")
	message(FATAL_ERROR "The ROOT version installed (${ROOT_VERSION}) is wront, the major version needs to be 6 for GRSISort to work (ROOT_VERSION_MAJOR=${ROOT_VERSION_MAJOR})!")
endif()

#----------------------------------------------------------------------------
# set c++ standard to the one used for root (adopted from agataselector's CMakeLists.txt)
if(NOT ROOT_CXX_VERSION)
	if(ROOT_CXX_STANDARD)
      set(ROOT_CXX_VERSION ${ROOT_CXX_STANDARD})
	else()
		string(REGEX REPLACE ".*-std=.+\\+\\+([0-9][0-9]).*" "\\1" ROOT_CXX_VERSION ${ROOT_CXX_FLAGS})
		if(${CMAKE_MATCH_COUNT} EQUAL "1")
			message("Got ROOT_CXX_VERSION=${ROOT_CXX_VERSION} from \"${ROOT_CXX_FLAGS}\"")
		else()
			message(FATAL_ERROR "The ROOT cxx flags seem to be misformed (${ROOT_CXX_FLAGS}) (CMAKE_MATCH_COUNT=${CMAKE_MATCH_COUNT}), consider providing a cxx version by hand using \"-DROOT_CXX_VERSION=??\"")
		endif()
	endif()
else()
	message("Got user provided ROOT_CXX_VERSION=${ROOT_CXX_VERSION}")
endif()

if(${ROOT_CXX_VERSION} LESS 14)
	message(FATAL_ERROR "${CMAKE_PROJECT_NAME} requires at least c++14, please consider installing a newer ROOT version that was compiled with at least c++14")
endif()

set(CMAKE_CXX_STANDARD ${ROOT_CXX_VERSION})
message("${BoldBlue}C++ version set to: ${CMAKE_CXX_STANDARD}\n${ColourReset}")
set(CMAKE_CXX_STANDARD_REQUIRED True)

#----------------------------------------------------------------------------
# configure the include directories and set other variables (could maybe have just one include_directories for both linux and macOS?)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${X11_INCLUDE_DIR})
	add_compile_options(-DOS_DARWIN)
else()
	include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_BINDIR ${CMAKE_BINARY_DIR}/bin)

#----------------------------------------------------------------------------
# copy some necessary text files
# this will always overwrite the files, including any changes in them
# tried using file(COPY, file(GENERATE, and add_custom_command to get around this for grsirc, no success
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/thisgrsi.sh  ${CMAKE_CURRENT_BINARY_DIR}/thisgrsi.sh  COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/thisgrsi.csh ${CMAKE_CURRENT_BINARY_DIR}/thisgrsi.csh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/.grsi.mimes  ${CMAKE_CURRENT_BINARY_DIR}/.grsi.mimes  COPYONLY)
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/.grsirc")
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/.grsirc   ${CMAKE_CURRENT_BINARY_DIR}/.grsirc      COPYONLY)
endif()

#----------------------------------------------------------------------------
# copy includes, examples, and SourceData
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/libraries/TAnalysis/SourceData DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/libraries/TAnalysis/SRIMData DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# define all libraries and root dictionaries
# for some we also create dependencies on other libraries to remove linking errors later on

add_library(MakeFragmentHistograms SHARED
	${PROJECT_SOURCE_DIR}/histos/MakeFragmentHistograms.cxx
	)
target_link_libraries(MakeFragmentHistograms TFormat ${ROOT_LIBRARIES})

add_library(TBetaDecay SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TBetaDecay/TBetaDecay.cxx
	)
root_generate_dictionary(G__TBetaDecay TBetaDecay.h MODULE TBetaDecay LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TBetaDecay/LinkDef.h)
target_link_libraries(TBetaDecay TNucleus ${ROOT_LIBRARIES})

add_library(TBgo SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TBgo/TBgo.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TBgo/TBgoHit.cxx
	)
root_generate_dictionary(G__TBgo TBgo.h TBgoHit.h MODULE TBgo LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TBgo/LinkDef.h)
target_link_libraries(TBgo TFormat ${ROOT_LIBRARIES})

add_library(TCal SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCFDCal.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCalGraph.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TTimeCal.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TEfficiencyGraph.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCalList.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TEfficiencyCalibration.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCalPoint.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TEfficiencyCalibrator.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCal.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TEnergyCal.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TSourceList.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TGainMatch.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TEfficiencyCal.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/TCalManager.cxx
	)
root_generate_dictionary(G__TCal TCal.h TCalManager.h TCFDCal.h TEfficiencyCal.h TEnergyCal.h TGainMatch.h TTimeCal.h TCalPoint.h TCalList.h TSourceList.h TCalGraph.h TEfficiencyGraph.h TEfficiencyCalibration.h TEfficiencyCalibrator.h MODULE TCal LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TCal/LinkDef.h)
target_link_libraries(TCal TFormat TGRSIFit TGUI TNucleus TPeakFitting ${ROOT_LIBRARIES})

add_library(TDataParser SHARED
	${PROJECT_SOURCE_DIR}/libraries/TDataParser/TFragmentMap.cxx
	${PROJECT_SOURCE_DIR}/libraries/TDataParser/TDataParser.cxx
	)
target_link_libraries(TDataParser TFormat ${ROOT_LIBRARIES})

add_library(TFormat SHARED
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TMnemonic.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/Globals.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TGRSIUtilities.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TGRSISortInfo.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TSortingDiagnostics.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TDataFrameLibrary.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TPPG.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TScalerQueue.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TScaler.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TChannel.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TParsingDiagnostics.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TRunInfo.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TParserLibrary.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TUserSettings.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TFragment.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TEpicsFrag.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TBadFragment.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TDetector.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TDetectorHit.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/GValue.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TGRSIFrame.cxx
	${PROJECT_SOURCE_DIR}/libraries/TFormat/TGRSIHelper.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/TGRSIint.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/grsixx.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/TAnalysisOptions.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/ArgParser.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/TGRSIOptions.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGRSIint/FullPath.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TTreeFillMutex.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TUnpackedEvent.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TUnpackingLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TDataLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/StoppableThread.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TEventBuildingLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TAnalysisWriteLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TFragWriteLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TDetBuildingLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TFragmentChainLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TRuntimeObjects.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/DynamicLibrary.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TCompiledHistograms.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TFragHistLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/TLoops/TAnalysisHistLoop.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GSnapshot.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GHSym.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GRootBrowser.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GPopup.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GRootCommands.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GNotifier.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GH2I.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GCube.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GRootCanvas.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GH1D.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GGaus.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GH2D.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GCanvas.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/TCalibrator.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GPeak.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GRootFunctions.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GH2Base.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/GRootGuiFactory.cxx
	${PROJECT_SOURCE_DIR}/libraries/GROOT/TLevelScheme.cxx
	)
root_generate_dictionary(G__TFormat TSingleton.h TFragment.h TBadFragment.h TChannel.h TRunInfo.h TGRSISortInfo.h TPPG.h TEpicsFrag.h TScaler.h TScalerQueue.h TParsingDiagnostics.h TGRSIUtilities.h TMnemonic.h TSortingDiagnostics.h TTransientBits.h TPriorityValue.h TSingleton.h TDetectorInformation.h TParserLibrary.h TDataFrameLibrary.h TUserSettings.h TDetector.h TDetectorHit.h GValue.h TGRSIFrame.h TGRSIHelper.h TGRSIOptions.h TGRSIint.h TAnalysisOptions.h TDataLoop.h StoppableThread.h TFragWriteLoop.h TTerminalLoop.h TEventBuildingLoop.h TDetBuildingLoop.h TAnalysisWriteLoop.h TFragHistLoop.h TCompiledHistograms.h TRuntimeObjects.h TAnalysisHistLoop.h GRootGuiFactory.h GRootFunctions.h GRootCommands.h GRootCanvas.h GRootBrowser.h GCanvas.h GH2Base.h GH2I.h GH2D.h GPeak.h GGaus.h GH1D.h GNotifier.h GPopup.h GSnapshot.h TCalibrator.h GHSym.h GCube.h TLevelScheme.h MODULE TFormat LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TFormat/CMakeLinkDef.h)
target_link_libraries(TFormat TGRSIFit TNucleus ${ROOT_LIBRARIES} ${X11_LIBRARIES} ${X11_Xpm_LIB})
add_dependencies(TFormat GVersion)

if(${ROOT_proof_FOUND})
	add_library(TGRSIProof SHARED
		${PROJECT_SOURCE_DIR}/libraries/TGRSIProof/TGRSISelector.cxx
		)
	root_generate_dictionary(G__TGRSIProof TGRSIProof.h TObjectWrapper.h TGRSISelector.h MODULE TGRSIProof LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TGRSIProof/LinkDef.h)
	target_link_libraries(TGRSIProof TFormat ${ROOT_LIBRARIES})
endif()

add_library(TGUI SHARED
	${PROJECT_SOURCE_DIR}/libraries/TGUI/TCalibrationGraph.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGUI/TBGSubtraction.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGUI/TCalibrateDescant.cxx
	${PROJECT_SOURCE_DIR}/libraries/TGUI/TSourceCalibration.cxx
	)
root_generate_dictionary(G__TGUI TBGSubtraction.h TCalibrateDescant.h TSourceCalibration.h TCalibrationGraph.h MODULE TGUI LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TGUI/LinkDef.h)
target_link_libraries(TGUI TFormat TNucleus TPeakFitting ${ROOT_LIBRARIES})

add_library(TKinematics SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TKinematics/TKinematics.cxx
	)
root_generate_dictionary(G__TKinematics TKinematics.h MODULE TKinematics LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TKinematics/LinkDef.h)
target_link_libraries(TKinematics TNucleus ${ROOT_LIBRARIES})

add_library(TPeakFitting SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TRWPeak.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TPeakFitter.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TGauss.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TAB3Peak.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TABPeak.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/TSinglePeak.cxx
	)
root_generate_dictionary(G__TPeakFitting TSinglePeak.h TRWPeak.h TPeakFitter.h TABPeak.h TAB3Peak.h TGauss.h MODULE TPeakFitting LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPeakFitting/LinkDef.h)
target_link_libraries(TPeakFitting TGRSIFit ${ROOT_LIBRARIES})

add_library(TPulseAnalyzer SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPulseAnalyzer/TPulseAnalyzer.cxx
	)
root_generate_dictionary(G__TPulseAnalyzer TPulseAnalyzer.h MODULE TPulseAnalyzer LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TPulseAnalyzer/LinkDef.h)
target_link_libraries(TPulseAnalyzer TGRSIFit ${ROOT_LIBRARIES})

add_library(TReaction SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TReaction/TReaction.cxx
	)
root_generate_dictionary(G__TReaction TReaction.h MODULE TReaction LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TReaction/LinkDef.h)
target_link_libraries(TReaction TNucleus ${ROOT_LIBRARIES})

add_library(TSuppressed SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TSuppressed/TSuppressed.cxx
	)
root_generate_dictionary(G__TSuppressed TSuppressed.h MODULE TSuppressed LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TSuppressed/LinkDef.h)
target_link_libraries(TSuppressed TFormat ${ROOT_LIBRARIES})

add_library(TGRSIFit SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TDecay.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TMultiPeak.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TGRSIFit.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TLMFitter.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TPeak.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/TGRSIFunctions.cxx
	)
root_generate_dictionary(G__TGRSIFit TGRSIFit.h TGRSIFunctions.h TMultiPeak.h TPeak.h TDecay.h TLMFitter.h MODULE TGRSIFit LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TGRSIFit/LinkDef.h)
target_link_libraries(TGRSIFit ${ROOT_LIBRARIES})

add_library(TNucleus SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TNucleus/TNucleus.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TNucleus/TTransition.cxx
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TNucleus/TGRSITransition.cxx
	)
root_generate_dictionary(G__TNucleus TNucleus.h TTransition.h TGRSITransition.h MODULE TNucleus LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TNucleus/LinkDef.h)
target_link_libraries(TNucleus ${ROOT_LIBRARIES})

add_library(TSRIM SHARED
	${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TSRIM/TSRIM.cxx
	)
root_generate_dictionary(G__TSRIM TSRIM.h MODULE TSRIM LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TAnalysis/TSRIM/LinkDef.h)
target_link_libraries(TSRIM TFormat ${ROOT_LIBRARIES})

#add_library(TRawFile INTERFACE)
add_library(TRawFile SHARED
	)
root_generate_dictionary(G__TRawFile TRawFile.h TRawEvent.h MODULE TRawFile LINKDEF ${PROJECT_SOURCE_DIR}/libraries/TRawFile/LinkDef.h)
target_link_libraries(TRawFile ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the generated libraries and the ROOT libraries
add_executable(grsisort ${PROJECT_SOURCE_DIR}/src/grsisort.cxx)

if(${ROOT_proof_FOUND})
	set(GRSI_LIBRARIES TFormat TPulseAnalyzer TSuppressed TReaction TPeakFitting TGRSIFit TKinematics TNucleus TBgo TBetaDecay TSRIM TCal TGUI TGRSIProof TDataParser TRawFile MakeFragmentHistograms)
else()
	set(GRSI_LIBRARIES TFormat TPulseAnalyzer TSuppressed TReaction TPeakFitting TGRSIFit TKinematics TNucleus TBgo TBetaDecay TSRIM TCal TGUI TDataParser TRawFile MakeFragmentHistograms)
endif()
target_link_libraries(grsisort ${GRSI_LIBRARIES} ${ROOT_LIBRARIES} ${X11_LIBRARIES} ${X11_Xpm_LIB})

#----------------------------------------------------------------------------
# add all executable in util
set(UTIL_NAMES AddToChannelNumber AnalyzeDataLoss checkSort Compton_Scale_Plot converTH1 EfficCal EnableCrosstalkCorrection EnergyCal gadd grsiframe LeanAnalyzeComptonMatrices PlotVsRun ReadCalFromRoot Root2Ascii Root2Rad ScaleCalibration SetOdbFromData WriteCalToRoot WriteNonlinearitiesCal)
foreach(UTIL IN LISTS UTIL_NAMES)
	add_executable(${UTIL} ${PROJECT_SOURCE_DIR}/util/${UTIL}.cxx)
	target_link_libraries(${UTIL} ${GRSI_LIBRARIES} ${ROOT_LIBRARIES} ${X11_LIBRARIES} ${X11_Xpm_LIB})
endforeach()

#----------------------------------------------------------------------------
# copy scripts
configure_file(${PROJECT_SOURCE_DIR}/util/grsi-config    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grsi-config COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/util/ErrorReport.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ErrorReport.sh COPYONLY)

#----------------------------------------------------------------------------
# clean up all copied files and directories
# we're using grsisort as target here, because most (all?) of these do not belong to a specific target
set_property(TARGET grsisort APPEND PROPERTY ADDITIONAL_CLEAN_FILES
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grsi-config
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ErrorReport.sh
	${CMAKE_CURRENT_BINARY_DIR}/thisgrsi.sh
	${CMAKE_CURRENT_BINARY_DIR}/thisgrsi.csh
	${CMAKE_CURRENT_BINARY_DIR}/.grsirc
	${CMAKE_CURRENT_BINARY_DIR}/.grsi.mimes
	${CMAKE_CURRENT_BINARY_DIR}/examples
	${CMAKE_CURRENT_BINARY_DIR}/module.modulemap
	${CMAKE_INSTALL_BINDIR}
	${CMAKE_INSTALL_LIBDIR}
	)

#----------------------------------------------------------------------------
# custom target to generate the version header
add_custom_target(GVersion
	${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/include/GVersion.h.in
	                 -D DST=${CMAKE_SOURCE_DIR}/include/GVersion.h
                    -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                    -P ${CMAKE_SOURCE_DIR}/GenerateVersionHeader.cmake
						  BYPRODUCTS ${CMAKE_SOURCE_DIR}/include/GVersion.h
						  COMMENT "Generating ${CMAKE_SOURCE_DIR}/include/GVersion.h from ${CMAKE_SOURCE_DIR}/include/GVersion.h.in"
						  #WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )

add_custom_command(TARGET GVersion POST_BUILD
	  COMMAND ${CMAKE_COMMAND} -E copy
	  ${CMAKE_SOURCE_DIR}/include/GVersion.h
	  ${CMAKE_CURRENT_BINARY_DIR}/include/GVersion.h
)

#----------------------------------------------------------------------------
# add grsiproof (but exclude it from target all)
if(${ROOT_proof_FOUND})
	add_executable(grsiproof EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/GRSIProof/grsiproof.cxx)
	target_link_libraries(grsiproof ${GRSI_LIBRARIES} ${ROOT_LIBRARIES} ${X11_LIBRARIES} ${X11_Xpm_LIB})
endif()

#----------------------------------------------------------------------------
# all parser libraries

#----------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
	GRSIData
	EXCLUDE_FROM_ALL
	GIT_REPOSITORY https://github.com/GRIFFINCollaboration/GRSIData.git
	GIT_TAG main
	GIT_REMOTE_NAME upstream
	)
FetchContent_MakeAvailable(
	GRSIData
	)
#set_property(TARGET GRSIData APPEND PROPERTY ADDITIONAL_CLEAN_FILES
#	${CMAKE_CURRENT_BINARY_DIR}/_deps
#	)
add_dependencies(GRSIData grsisort ${UTIL_NAMES})
if(${ROOT_proof_FOUND})
	add_dependencies(GRSIData grsiproof)
endif()

#----------------------------------------------------------------------------
FetchContent_Declare(
	ILLData
	EXCLUDE_FROM_ALL
	GIT_REPOSITORY https://github.com/UoG-Nuclear-Physics-Group/ILLData.git
	GIT_TAG main
	GIT_REMOTE_NAME upstream
	)
FetchContent_MakeAvailable(
	ILLData
	)
add_dependencies(ILLData grsisort ${UTIL_NAMES})
if(${ROOT_proof_FOUND})
	add_dependencies(ILLData grsiproof)
endif()

#----------------------------------------------------------------------------
FetchContent_Declare(
	HILData
	EXCLUDE_FROM_ALL
	GIT_REPOSITORY https://github.com/UoG-Nuclear-Physics-Group/HILData.git
	GIT_TAG main
	GIT_REMOTE_NAME upstream
	)
FetchContent_MakeAvailable(
	HILData
	)
add_dependencies(HILData grsisort ${UTIL_NAMES})
if(${ROOT_proof_FOUND})
	add_dependencies(HILData grsiproof)
endif()

#----------------------------------------------------------------------------
FetchContent_Declare(
	iThembaData
	EXCLUDE_FROM_ALL
	GIT_REPOSITORY https://github.com/UoG-Nuclear-Physics-Group/iThembaData.git
	GIT_TAG main
	GIT_REMOTE_NAME upstream
	)
FetchContent_MakeAvailable(
	iThembaData
	)
add_dependencies(iThembaData grsisort ${UTIL_NAMES})
if(${ROOT_proof_FOUND})
	add_dependencies(iThembaData grsiproof)
endif()

#----------------------------------------------------------------------------
# create a compile_commands.json file which can be used for clang-tidy
if(NOT DEFINED ENV{CMAKE_EXPORT_COMPILE_COMMANDS})
	message("ENV(CMAKE_EXPORT_COMPILE_COMMANDS) NOT DEFINED, setting it to TRUE")
	set($ENV{CMAKE_EXPORT_COMPILE_COMMANDS} TRUE)
endif()
