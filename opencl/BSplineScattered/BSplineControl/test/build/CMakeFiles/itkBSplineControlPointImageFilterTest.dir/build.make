# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yao/works/itk-work/BSplineScattered/BSplineControl/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build

# Include any dependencies generated for this target.
include CMakeFiles/itkBSplineControlPointImageFilterTest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/itkBSplineControlPointImageFilterTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/itkBSplineControlPointImageFilterTest.dir/flags.make

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/flags.make
CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o: ../itkBSplineControlPointImageFilterTest.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o -c /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/itkBSplineControlPointImageFilterTest.cxx

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/itkBSplineControlPointImageFilterTest.cxx > CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.i

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/itkBSplineControlPointImageFilterTest.cxx -o CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.s

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.requires:
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.requires

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.provides: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.requires
	$(MAKE) -f CMakeFiles/itkBSplineControlPointImageFilterTest.dir/build.make CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.provides.build
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.provides

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.provides.build: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o

# Object files for target itkBSplineControlPointImageFilterTest
itkBSplineControlPointImageFilterTest_OBJECTS = \
"CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o"

# External object files for target itkBSplineControlPointImageFilterTest
itkBSplineControlPointImageFilterTest_EXTERNAL_OBJECTS =

itkBSplineControlPointImageFilterTest: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o
itkBSplineControlPointImageFilterTest: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/build.make
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkdouble-conversion-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitksys-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkvnl_algo-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkvnl-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkv3p_netlib-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKCommon-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkNetlibSlatec-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKStatistics-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOImageBase-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOBMP-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOBioRad-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKEXPAT-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkopenjpeg-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkzlib-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmDICT-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmMSFF-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGDCM-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGIPL-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkjpeg-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOJPEG-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitktiff-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTIFF-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOLSM-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKMetaIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOMeta-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKznz-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKniftiio-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIONIFTI-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKNrrdIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIONRRD-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkpng-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOPNG-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOStimulate-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOVTK-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKMesh-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKSpatialObjects-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKPath-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKLabelMap-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKQuadEdgeMesh-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKOptimizers-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKPolynomials-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKBiasCorrection-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKBioCell-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKDICOMParser-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOXML-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOSpatialObjects-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKFEM-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKgiftiio-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOMesh-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkhdf5_cpp-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkhdf5-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOCSV-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOIPL-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGE-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOSiemens-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOHDF5-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformBase-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformHDF5-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformInsightLegacy-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformMatlab-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKKLMRegionGrowing-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVTK-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKWatersheds-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKOptimizersv4-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKReview-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVideoCore-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVideoIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVideoBridgeOpenCV-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKgiftiio-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOBMP-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOBioRad-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGDCM-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmMSFF-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkopenjpeg-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmDICT-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmIOD-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmDSED-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkgdcmCommon-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGIPL-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOJPEG-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTIFF-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitktiff-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkjpeg-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOMeta-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIONIFTI-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKniftiio-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKznz-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIONRRD-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKNrrdIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOPNG-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkpng-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOStimulate-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOVTK-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKLabelMap-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKQuadEdgeMesh-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKBiasCorrection-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKPolynomials-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKBioCell-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOSpatialObjects-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOXML-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKEXPAT-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKFEM-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKMetaIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKOptimizers-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOSiemens-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOGE-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOIPL-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformHDF5-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkhdf5_cpp-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkhdf5-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkzlib-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformInsightLegacy-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformMatlab-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOTransformBase-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKKLMRegionGrowing-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVTK-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKWatersheds-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKSpatialObjects-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKMesh-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKPath-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKStatistics-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkNetlibSlatec-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVideoIO-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKIOImageBase-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVideoCore-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKCommon-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkdouble-conversion-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitksys-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libITKVNLInstantiation-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkvnl_algo-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkv3p_lsqr-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkvnl-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkvcl-4.4.so.1
itkBSplineControlPointImageFilterTest: /tmp/itk-build/lib/libitkv3p_netlib-4.4.so.1
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_calib3d.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_contrib.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_core.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_features2d.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_flann.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_gpu.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_highgui.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_imgproc.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_legacy.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_ml.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_nonfree.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_objdetect.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_ocl.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_photo.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_stitching.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_superres.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_ts.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_video.so
itkBSplineControlPointImageFilterTest: /usr/local/lib/libopencv_videostab.so
itkBSplineControlPointImageFilterTest: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable itkBSplineControlPointImageFilterTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/itkBSplineControlPointImageFilterTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/itkBSplineControlPointImageFilterTest.dir/build: itkBSplineControlPointImageFilterTest
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/build

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/requires: CMakeFiles/itkBSplineControlPointImageFilterTest.dir/itkBSplineControlPointImageFilterTest.cxx.o.requires
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/requires

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/itkBSplineControlPointImageFilterTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/clean

CMakeFiles/itkBSplineControlPointImageFilterTest.dir/depend:
	cd /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yao/works/itk-work/BSplineScattered/BSplineControl/test /home/yao/works/itk-work/BSplineScattered/BSplineControl/test /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build /home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build/CMakeFiles/itkBSplineControlPointImageFilterTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/itkBSplineControlPointImageFilterTest.dir/depend

