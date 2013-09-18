cmake_minimum_required(VERSION 2.8)
project (itkBSplineControlPointImageFilterTest)

set(ITK_DIR /home/yao/works/itk-build)
find_package(ITK 4 REQUIRED)
include (${ITK_USE_FILE})

add_definitions(-DDEBUG)

add_executable(itkBSplineControlPointImageFilterTest
  itkBSplineControlPointImageFilterTest.cxx)
target_link_libraries(itkBSplineControlPointImageFilterTest ${ITK_LIBRARIES})
