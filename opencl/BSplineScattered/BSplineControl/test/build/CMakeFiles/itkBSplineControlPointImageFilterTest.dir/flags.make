# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# compile CXX with /usr/bin/c++
CXX_FLAGS =   -msse2 -g -I/home/yao/works/itk-work/BSplineScattered/BSplineControl/test/build/ITKIOFactoryRegistration -I/home/yao/Downloads/ITK/Modules/Video/Filtering/include -I/usr/local/include -I/usr/local/include/opencv -I/home/yao/Downloads/ITK/Modules/Video/BridgeOpenCV/include -I/home/yao/Downloads/ITK/Modules/Video/IO/include -I/home/yao/Downloads/ITK/Modules/Video/Core/include -I/home/yao/Downloads/ITK/Modules/Nonunit/Review/include -I/home/yao/Downloads/ITK/Modules/Registration/RegistrationMethodsv4/include -I/home/yao/Downloads/ITK/Modules/Registration/Metricsv4/include -I/home/yao/Downloads/ITK/Modules/Numerics/Optimizersv4/include -I/home/yao/Downloads/ITK/Modules/Segmentation/LevelSetsv4/include -I/home/yao/Downloads/ITK/Modules/Segmentation/Watersheds/include -I/home/yao/Downloads/ITK/Modules/Segmentation/Voronoi/include -I/home/yao/Downloads/ITK/Modules/Bridge/VTK/include -I/home/yao/Downloads/ITK/Modules/Filtering/SpatialFunction/include -I/home/yao/Downloads/ITK/Modules/Segmentation/RegionGrowing/include -I/home/yao/Downloads/ITK/Modules/Filtering/QuadEdgeMeshFiltering/include -I/home/yao/Downloads/ITK/Modules/Numerics/NeuralNetworks/include -I/home/yao/Downloads/ITK/Modules/Segmentation/MarkovRandomFieldsClassifiers/include -I/home/yao/Downloads/ITK/Modules/Segmentation/LabelVoting/include -I/home/yao/Downloads/ITK/Modules/Segmentation/KLMRegionGrowing/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageFusion/include -I/home/yao/Downloads/ITK/Modules/IO/TransformMatlab/include -I/home/yao/Downloads/ITK/Modules/IO/TransformInsightLegacy/include -I/home/yao/Downloads/ITK/Modules/IO/TransformHDF5/include -I/home/yao/Downloads/ITK/Modules/IO/TransformBase/include -I/home/yao/Downloads/ITK/Modules/IO/RAW/include -I/home/yao/Downloads/ITK/Modules/IO/HDF5/include -I/home/yao/Downloads/ITK/Modules/IO/Siemens/include -I/home/yao/Downloads/ITK/Modules/IO/GE/include -I/home/yao/Downloads/ITK/Modules/IO/IPL/include -I/home/yao/Downloads/ITK/Modules/IO/CSV/include -I/tmp/itk-build/Modules/ThirdParty/HDF5/src -I/home/yao/Downloads/ITK/Modules/ThirdParty/HDF5/src -I/home/yao/Downloads/ITK/Modules/Filtering/GPUThresholding/include -I/home/yao/Downloads/ITK/Modules/Filtering/GPUSmoothing/include -I/home/yao/Downloads/ITK/Modules/Registration/GPUPDEDeformable/include -I/home/yao/Downloads/ITK/Modules/Registration/GPUCommon/include -I/home/yao/Downloads/ITK/Modules/Filtering/GPUImageFilterBase/include -I/home/yao/Downloads/ITK/Modules/Filtering/GPUAnisotropicSmoothing/include -I/home/yao/Downloads/ITK/Modules/Core/GPUFiniteDifference/include -I/home/yao/Downloads/ITK/Modules/Core/GPUCommon/include -I/home/yao/Downloads/ITK/Modules/IO/Mesh/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/GIFTI/src/gifticlib -I/home/yao/Downloads/ITK/Modules/Registration/FEM/include -I/home/yao/Downloads/ITK/Modules/Registration/PDEDeformable/include -I/home/yao/Downloads/ITK/Modules/Numerics/FEM/include -I/home/yao/Downloads/ITK/Modules/Registration/Common/include -I/home/yao/Downloads/ITK/Modules/IO/SpatialObjects/include -I/home/yao/Downloads/ITK/Modules/IO/XML/include -I/home/yao/Downloads/ITK/Modules/Numerics/Eigen/include -I/home/yao/Downloads/ITK/Modules/Filtering/DisplacementField/include -I/home/yao/Downloads/ITK/Modules/Filtering/DiffusionTensorImage/include -I/home/yao/Downloads/ITK/Modules/Filtering/Denoising/include -I/home/yao/Downloads/ITK/Modules/Segmentation/DeformableMesh/include -I/home/yao/Downloads/ITK/Modules/Filtering/Deconvolution/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/DICOMParser/src/DICOMParser -I/tmp/itk-build/Modules/ThirdParty/DICOMParser/src/DICOMParser -I/home/yao/Downloads/ITK/Modules/Filtering/Convolution/include -I/home/yao/Downloads/ITK/Modules/Filtering/FFT/include -I/home/yao/Downloads/ITK/Modules/Filtering/Colormap/include -I/home/yao/Downloads/ITK/Modules/Segmentation/Classifiers/include -I/home/yao/Downloads/ITK/Modules/Segmentation/BioCell/include -I/home/yao/Downloads/ITK/Modules/Filtering/BiasCorrection/include -I/home/yao/Downloads/ITK/Modules/Numerics/Polynomials/include -I/home/yao/Downloads/ITK/Modules/Filtering/AntiAlias/include -I/home/yao/Downloads/ITK/Modules/Segmentation/LevelSets/include -I/home/yao/Downloads/ITK/Modules/Segmentation/SignedDistanceFunction/include -I/home/yao/Downloads/ITK/Modules/Numerics/Optimizers/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageFeature/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageSources/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageGradient/include -I/home/yao/Downloads/ITK/Modules/Filtering/Smoothing/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageCompare/include -I/home/yao/Downloads/ITK/Modules/Filtering/FastMarching/include -I/home/yao/Downloads/ITK/Modules/Core/QuadEdgeMesh/include -I/home/yao/Downloads/ITK/Modules/Filtering/DistanceMap/include -I/home/yao/Downloads/ITK/Modules/Numerics/NarrowBand/include -I/home/yao/Downloads/ITK/Modules/Filtering/BinaryMathematicalMorphology/include -I/home/yao/Downloads/ITK/Modules/Filtering/LabelMap/include -I/home/yao/Downloads/ITK/Modules/Filtering/MathematicalMorphology/include -I/home/yao/Downloads/ITK/Modules/Segmentation/ConnectedComponents/include -I/home/yao/Downloads/ITK/Modules/Filtering/Thresholding/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageLabel/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageIntensity/include -I/home/yao/Downloads/ITK/Modules/Filtering/Path/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageStatistics/include -I/home/yao/Downloads/ITK/Modules/Core/SpatialObjects/include -I/home/yao/Downloads/ITK/Modules/Core/Mesh/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageCompose/include -I/home/yao/Downloads/ITK/Modules/Core/TestKernel/include -I/home/yao/Downloads/ITK/Modules/IO/VTK/include -I/home/yao/Downloads/ITK/Modules/IO/Stimulate/include -I/home/yao/Downloads/ITK/Modules/IO/PNG/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/PNG/src -I/tmp/itk-build/Modules/ThirdParty/PNG/src -I/home/yao/Downloads/ITK/Modules/IO/NRRD/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/NrrdIO/src/NrrdIO -I/tmp/itk-build/Modules/ThirdParty/NrrdIO/src/NrrdIO -I/home/yao/Downloads/ITK/Modules/IO/NIFTI/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/NIFTI/src/nifti/znzlib -I/home/yao/Downloads/ITK/Modules/ThirdParty/NIFTI/src/nifti/niftilib -I/home/yao/Downloads/ITK/Modules/IO/Meta/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/MetaIO/src/MetaIO -I/tmp/itk-build/Modules/ThirdParty/MetaIO/src/MetaIO -I/home/yao/Downloads/ITK/Modules/IO/LSM/include -I/home/yao/Downloads/ITK/Modules/IO/TIFF/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/TIFF/src -I/tmp/itk-build/Modules/ThirdParty/TIFF/src/itktiff -I/tmp/itk-build/Modules/ThirdParty/TIFF/src -I/home/yao/Downloads/ITK/Modules/IO/JPEG/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/JPEG/src -I/tmp/itk-build/Modules/ThirdParty/JPEG/src -I/home/yao/Downloads/ITK/Modules/IO/GIPL/include -I/home/yao/Downloads/ITK/Modules/IO/GDCM/include -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/DataStructureAndEncodingDefinition -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/MessageExchangeDefinition -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/InformationObjectDefinition -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/Common -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/DataDictionary -I/home/yao/Downloads/ITK/Modules/ThirdParty/GDCM/src/gdcm/Source/MediaStorageAndFileFormat -I/tmp/itk-build/Modules/ThirdParty/GDCM/src/gdcm/Source/Common -I/tmp/itk-build/Modules/ThirdParty/GDCM -I/home/yao/Downloads/ITK/Modules/ThirdParty/ZLIB/src -I/tmp/itk-build/Modules/ThirdParty/ZLIB/src -I/home/yao/Downloads/ITK/Modules/ThirdParty/OpenJPEG/src/openjpeg -I/tmp/itk-build/Modules/ThirdParty/OpenJPEG/src/openjpeg -I/home/yao/Downloads/ITK/Modules/ThirdParty/Expat/src/expat -I/tmp/itk-build/Modules/ThirdParty/Expat/src/expat -I/home/yao/Downloads/ITK/Modules/IO/BioRad/include -I/home/yao/Downloads/ITK/Modules/IO/BMP/include -I/home/yao/Downloads/ITK/Modules/IO/ImageBase/include -I/tmp/itk-build/Modules/IO/ImageBase -I/home/yao/Downloads/ITK/Modules/Filtering/AnisotropicSmoothing/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageGrid/include -I/home/yao/Downloads/ITK/Modules/Core/ImageFunction/include -I/home/yao/Downloads/ITK/Modules/Core/Transform/include -I/home/yao/Downloads/ITK/Modules/Numerics/Statistics/include -I/tmp/itk-build/Modules/ThirdParty/Netlib -I/home/yao/Downloads/ITK/Modules/Core/ImageAdaptors/include -I/home/yao/Downloads/ITK/Modules/Filtering/CurvatureFlow/include -I/home/yao/Downloads/ITK/Modules/Filtering/ImageFilterBase/include -I/home/yao/Downloads/ITK/Modules/Core/FiniteDifference/include -I/home/yao/Downloads/ITK/Modules/Core/Common/include -I/tmp/itk-build/Modules/Core/Common -I/home/yao/Downloads/ITK/Modules/ThirdParty/VNLInstantiation/include -I/tmp/itk-build/Modules/ThirdParty/VNL/src/vxl/core -I/tmp/itk-build/Modules/ThirdParty/VNL/src/vxl/vcl -I/tmp/itk-build/Modules/ThirdParty/VNL/src/vxl/v3p/netlib -I/home/yao/Downloads/ITK/Modules/ThirdParty/VNL/src/vxl/core -I/home/yao/Downloads/ITK/Modules/ThirdParty/VNL/src/vxl/vcl -I/home/yao/Downloads/ITK/Modules/ThirdParty/VNL/src/vxl/v3p/netlib -I/tmp/itk-build/Modules/ThirdParty/KWSys/src -I/home/yao/Downloads/ITK/Modules/ThirdParty/DoubleConversion/src/double-conversion   

CXX_DEFINES = -DDEBUG -DITK_IO_FACTORY_REGISTER_MANAGER

