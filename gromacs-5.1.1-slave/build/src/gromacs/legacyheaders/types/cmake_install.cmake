# Install script for directory: /home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/export/base/shisuan/swyaow/online/luhao/zx/APP/gmx511")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/legacyheaders/types" TYPE FILE FILES
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/commrec_fwd.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/constr.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/energy.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/enums.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/fcdata.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/force_flags.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/forcerec.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/genborn.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/group.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/hw_info.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/ifunc.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/inputrec.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/interaction_const.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/mdatom.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/membedt.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/nblist.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/nrnb.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/ns.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/nsgrid.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/oenv.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/qmmmrec.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/rgb.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/shellfc.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/simple.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/types/state.h"
    )
endif()

