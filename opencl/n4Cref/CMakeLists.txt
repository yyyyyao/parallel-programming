cmake_minimum_required(VERSION 2.8)
project (itkN4BiasFieldCorrectionImageFilterTest)

set(ITK_DIR /home/yao/works/itk-build)
find_package(ITK 4 REQUIRED)
include (${ITK_USE_FILE})

add_definitions(-DDEBUG)

add_executable(itkN4BiasFieldCorrectionImageFilterTest
  itkN4BiasFieldCorrectionImageFilterTest.cxx)
target_link_libraries(itkN4BiasFieldCorrectionImageFilterTest ${ITK_LIBRARIES})
