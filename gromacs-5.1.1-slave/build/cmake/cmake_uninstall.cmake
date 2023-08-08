IF(NOT EXISTS "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/install_manifest.txt")
  MESSAGE(FATAL_ERROR "Cannot find install manifest: \"/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/install_manifest.txt\"")
ENDIF(NOT EXISTS "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/install_manifest.txt")

FILE(READ "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/install_manifest.txt" files)
STRING(REGEX REPLACE "\n" ";" files "${files}")
FOREACH(file ${files})
  MESSAGE(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
  IF(EXISTS "$ENV{DESTDIR}${file}" OR IS_SYMLINK "$ENV{DESTDIR}${file}")
    EXEC_PROGRAM(
      "/home/export/online1/mdt00/shisuan/swyaow/luhao/zx/APP/cmake-3.21.0-rc2-linux-x86_64/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF(NOT "${rm_retval}" STREQUAL 0)
      MESSAGE(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ELSE(EXISTS "$ENV{DESTDIR}${file}" OR IS_SYMLINK "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File \"$ENV{DESTDIR}${file}\" does not exist.")
  ENDIF(EXISTS "$ENV{DESTDIR}${file}" OR IS_SYMLINK "$ENV{DESTDIR}${file}")
ENDFOREACH(file)

