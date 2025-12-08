# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2_add_library(FD3Workflow
               SOURCES src/DigitReaderSpec.cxx
                       src/FD3DigitWriterSpec.cxx
               PUBLIC_LINK_LIBRARIES O2::DataFormatsFD3
                                     O2::DetectorsCommonDataFormats
                                     O2::Framework
                                     O2::DPLUtils
                                     O2::CommonUtils)


o2_add_executable(digit-reader-workflow
                  SOURCES src/digits-reader-workflow.cxx
                  COMPONENT_NAME fd3
                  PUBLIC_LINK_LIBRARIES O2::FD3Workflow)

o2_add_executable(digits-writer-workflow
                  SOURCES src/digits-writer-workflow.cxx
                  COMPONENT_NAME fd3
                  PUBLIC_LINK_LIBRARIES O2::FD3Workflow)


if(NOT APPLE)

 set_property(TARGET ${fitrecoexe} PROPERTY LINK_WHAT_YOU_USE ON)

endif()

