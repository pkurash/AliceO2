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

o2_add_library(DataFormatsFD3
  SOURCES src/CalibParam.cxx
          src/ChannelData.cxx
	  src/CTF.cxx
	  src/Digit.cxx
	  src/DigitFilterParam.cxx
          src/Hit.cxx
          src/Triggers.cxx
	  src/RecPoints.cxx
          
  PUBLIC_LINK_LIBRARIES O2::FD3Base
                        O2::SimulationDataFormat
                        O2::CommonDataFormat
                        Microsoft.GSL::GSL
                        O2::DetectorsCommonDataFormats
)

o2_target_root_dictionary(DataFormatsFD3
  HEADERS include/DataFormatsFD3/CalibParam.h
	  include/DataFormatsFD3/ChannelData.h
          include/DataFormatsFD3/CTF.h
	  include/DataFormatsFD3/Digit.h
	  include/DataFormatsFD3/DigitFilterParam.h
          include/DataFormatsFD3/Hit.h
	  include/DataFormatsFD3/MCLabel.h
          include/DataFormatsFD3/Triggers.h
	  include/DataFormatsFD3/RecPoints.h
)
