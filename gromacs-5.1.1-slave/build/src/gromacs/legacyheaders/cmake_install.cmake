# Install script for directory: /home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/legacyheaders" TYPE FILE FILES
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/calcgrid.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/calcmu.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/chargegroup.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/checkpoint.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/constr.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/copyrite.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/disre.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/ebin.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/force.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/genborn.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/gmx_cpuid.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/gmx_detect_hardware.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/gmx_omp_nthreads.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/gmx_thread_affinity.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/inputrec.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/macros.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/main.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/md_logging.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/md_support.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/mdatoms.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/mdebin.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/mdrun.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/names.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/network.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/nonbonded.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/nrnb.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/ns.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/nsgrid.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/oenv.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/orires.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/perf_est.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/qmmm.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/rbin.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/readinp.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/shellfc.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/sighandler.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/sim_util.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/splitter.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/tables.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/tgroup.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/txtdump.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/typedefs.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/update.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/vcm.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/viewit.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/vsite.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/legacyheaders/warninp.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/src/gromacs/legacyheaders/types/cmake_install.cmake")

endif()

