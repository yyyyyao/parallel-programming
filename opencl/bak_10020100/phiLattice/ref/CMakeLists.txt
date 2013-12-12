cmake_minimum_required(VERSION 2.8)
project (itkBSplineScatteredDataPointSetToImageFilterTest)

set(ITK_DIR /home/yao/works/itk-build)
find_package(ITK 4 REQUIRED)
find_package(OpenCV REQUIRED)
include (${ITK_USE_FILE} ${OPENCV_INCLUDE_DIR})

add_definitions(-DDEBUG)

add_executable(itkBSplineScatteredDataPointSetToImageFilterTest
  itkBSplineScatteredDataPointSetToImageFilterTest.cxx)
target_link_libraries(itkBSplineScatteredDataPointSetToImageFilterTest ${ITK_LIBRARIES} ${OpenCV_LIBS})
