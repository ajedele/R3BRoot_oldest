##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019 Members of R3B Collaboration                          #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################

# Create a library called "libR3B" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(SYSTEM_INCLUDE_DIRECTORIES
${SYSTEM_INCLUDE_DIRECTORIES}
${BASE_INCLUDE_DIRECTORIES}
${ucesb_INCLUDE_DIR}
)
Set(R3BSOF_SOURCE_DIR ${R3BROOT_SOURCE_DIR}/sofia)

set(INCLUDE_DIRECTORIES
#put here all directories where header files are located
${R3BROOT_SOURCE_DIR}/r3bsource
${R3BROOT_SOURCE_DIR}/r3bsource/ext
${R3BROOT_SOURCE_DIR}/r3bbase
${R3BROOT_SOURCE_DIR}/r3bdata/beammonitorData
${R3BROOT_SOURCE_DIR}/r3bdata/califaData
${R3BROOT_SOURCE_DIR}/r3bdata/neulandData
${R3BROOT_SOURCE_DIR}/r3bdata/pspData
${R3BROOT_SOURCE_DIR}/r3bdata/fibData
${R3BROOT_SOURCE_DIR}/r3bdata/sfibData
${R3BROOT_SOURCE_DIR}/r3bdata/losData
${R3BROOT_SOURCE_DIR}/r3bdata/roluData
${R3BROOT_SOURCE_DIR}/r3bdata/sci2Data
${R3BROOT_SOURCE_DIR}/r3bdata/sci8Data
${R3BROOT_SOURCE_DIR}/r3bdata/tofData
${R3BROOT_SOURCE_DIR}/r3bdata/pdcData
${R3BROOT_SOURCE_DIR}/r3bdata/strawtubesData
#${R3BROOT_SOURCE_DIR}/psp/par
${R3BROOT_SOURCE_DIR}/r3bdata/startrackData
${R3BROOT_SOURCE_DIR}/r3bdata/amsData
${R3BROOT_SOURCE_DIR}/r3bdata/wrData
${R3BROOT_SOURCE_DIR}/r3bdata/musicData
${R3BROOT_SOURCE_DIR}/r3bdata/sampData
)

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${ucesb_LIBRARY_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})
link_directories( ${LINK_DIRECTORIES})

set(SRCS
R3BUcesbSource.cxx
R3BReader.cxx
R3BUnpackReader.cxx
#R3BWhiterabbitReader.cxx
R3BWhiterabbitMasterReader.cxx
R3BWhiterabbitCalifaReader.cxx
R3BWhiterabbitNeulandReader.cxx
R3BWhiterabbitAmsReader.cxx
R3BWhiterabbitS2Reader.cxx
R3BWhiterabbitS8Reader.cxx
R3BTrloiiTpatReader.cxx
R3BTrloiiSampReader.cxx
R3BTimestampMasterReader.cxx
R3BPspxReader.cxx
R3BPspxReader_s515.cxx
R3BBunchedFiberReader.cxx
R3BBunchedFiberSPMTTrigReader.cxx
R3BBeamMonitorReader.cxx
R3BFi0Reader.cxx
R3BFi1aReader.cxx
R3BFi1bReader.cxx
#R3BFi2aReader.cxx
#R3BFi2bReader.cxx
R3BFi3aReader.cxx
R3BFi3bReader.cxx
R3BFi4Reader.cxx
R3BFi5Reader.cxx
R3BFi6Reader.cxx
R3BFi7Reader.cxx
R3BFi8Reader.cxx
#R3BFi9Reader.cxx
R3BFi10Reader.cxx
R3BFi11Reader.cxx
R3BFi12Reader.cxx
R3BFi13Reader.cxx
R3BLosReader.cxx
R3BRoluReader.cxx
R3BSci2Reader.cxx
R3BSci8Reader.cxx
R3BSfibReader.cxx
R3BTofdReader.cxx
R3BPdcReader.cxx
R3BCalifaFebexReader.cxx
R3BPtofReader.cxx
R3BNeulandTamexReader.cxx
R3BNeulandTacquilaReader.cxx
R3BStrawtubesReader.cxx
R3BStartrackReader.cxx
R3BAmsReader.cxx
R3BMusicReader.cxx
)

Set(STRUCT_HEADERS
ext/ext_h101_unpack.h
ext/ext_h101_whiterabbit.h
ext/ext_h101_wrmaster.h
ext/ext_h101_wrcalifa.h
ext/ext_h101_wrneuland.h
ext/ext_h101_wrams.h
ext/ext_h101_wrs2.h
ext/ext_h101_wrs8.h
ext/ext_h101_tpat.h
ext/ext_h101_timestamp_master.h
ext/ext_h101_bmon.h
ext/ext_h101_fibzero.h
ext/ext_h101_fibonea.h
ext/ext_h101_fiboneb.h
#ext_h101_fibtwoa.h
#ext_h101_fibtwob.h
ext/ext_h101_fibthreea.h
ext/ext_h101_fibthreeb.h
ext/ext_h101_fibfour.h
ext/ext_h101_fibfive.h
ext/ext_h101_fibsix.h
ext/ext_h101_fibseven.h
ext/ext_h101_fibeight.h
#ext_h101_fibnine.h
ext/ext_h101_fibten.h
ext/ext_h101_fibeleven.h
ext/ext_h101_fibtwelve.h
ext/ext_h101_fibthirteen.h
ext/ext_h101_los.h
ext/ext_h101_rolu.h
ext/ext_h101_sci2.h
ext/ext_h101_sci8.h
ext/ext_h101_sfib.h
ext/ext_h101_tofd.h
ext/ext_h101_ptof.h
ext/ext_h101_pspx.h
ext/ext_h101_pspx_s515.h
ext/ext_h101_raw_califa_febex.h
ext/ext_h101_raw_nnp_tacquila.h
ext/ext_h101_raw_nnp_tamex.h
ext/ext_h101_strawtubes.h
ext/ext_h101_star.h
ext/ext_h101_fibfour.h
ext/ext_h101_ams.h
ext/ext_h101_music.h
ext/ext_h101_samp.h
ext/ext_h101_pdc.h
ext/ext_h101_samplos.h
ext/ext_h101_sampms.h
ext/ext_h101_samplosms.h
)

# fill list of header files from list of source files
# by exchanging the file extension
CHANGE_FILE_EXTENSION(*.cxx *.h HEADERS "${SRCS}")
Set(HEADERS ${STRUCT_HEADERS} ${HEADERS})

set(LINKDEF_FILE R3BLinkDef.h)
set(DEPENDENCIES
    GeoBase ParBase MbsAPI Base FairTools
    R3Bbase R3BData Core Geom GenVector
    Physics Matrix MathCore
    #R3BPsp
    ${ucesb_LIBRARY_SHARED})
set(LIBRARY_NAME R3Bsource)

If(ROOT_FOUND_VERSION LESS 59999)
  Set(DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/G__${LIBRARY_NAME}Dict.cxx)
  ROOT_GENERATE_DICTIONARY_OLD_EXTRA("${HEADERS}" "${LINKDEF_FILE}" "${DICTIONARY}" "${INCLUDE_DIRECTORIES}${SYSTEM_INCLUDE_DIRECTORIES}")
  Set(SRCS ${SRCS} ${DICTIONARY})
  CHANGE_FILE_EXTENSION(*.cxx *.h HEADERS "${SRCS}")
  Set(HEADERS ${STRUCT_HEADERS} ${HEADERS})
  Set(LINKDEF ${LINKDEF_FILE})
  ROOT_GENERATE_ROOTMAP()
  Set(LINKDEF )
Else(ROOT_FOUND_VERSION LESS 59999)
  Set(LINKDEF ${LINKDEF_FILE})
EndIf(ROOT_FOUND_VERSION LESS 59999)

GENERATE_LIBRARY()
