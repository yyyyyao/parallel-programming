/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkConstantPadImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkN4BiasFieldCorrectionImageFilterOrig.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkOpenCVImageBridge.h"

#include <cv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <itkMedianImageFilter.h>

#include <iostream>
#include <string>
#include "itkImageRegionIterator.h"
#include "itkImageScanlineConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {}

public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    const TFilter * filter =
      dynamic_cast< const TFilter * >( object );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
                                        { return; }
    if( filter->GetElapsedIterations() == 1 )
      {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1
                << std::endl;
      }
    std::cout << "  Iteration " << filter->GetElapsedIterations()
              << " (of "
              //<< filter->GetMaximumNumberOfIterations()[
      //filter->GetCurrentLevel()]
              << ").  ";
    std::cout << " Current convergence value = "
              << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
    }

};


template<class TValue>
TValue Convert( std::string optionString )
{
  TValue             value;
  std::istringstream iss( optionString );

  iss >> value;
  return value;
}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string        element = optionString.substr( 0, crosspos );
    TValue             value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos );
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

template<unsigned int ImageDimension>
int N4( int argc, char *argv[] )
{
  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension>  ImageType;
  typedef typename ImageType::Pointer           ImagePointer;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  ImagePointer inputImage = reader->GetOutput();
  inputImage->DisconnectPipeline();

  printf("----------------------------------\n");
  itk::ImageScanlineConstIterator<ImageType> it(inputImage, inputImage->GetRequestedRegion());
  it.GoToBegin();
  while(!it.IsAtEnd()) {
    while(!it.IsAtEndOfLine()) {
      std::cout << it.Get() << " ";
      ++it;
    }
    it.NextLine();
    std::cout << std::endl;
  }

#if 0
  //testing itk to opencv 
  typedef itk::Image<unsigned char, 2>  CvImageType;
  typedef itk::OpenCVImageBridge  BridgeType;  
  typedef itk::MedianImageFilter <CvImageType, CvImageType> CvFilterType;

  cv::Mat cvInput = BridgeType::ITKImageToCVMat<ImageType>(inputImage);
  cv::Mat cvFloatInput;
  cvFloatInput.create(cvInput.cols, cvInput.rows, CV_32FC1);
  cv::Mat floatImg;
  cvInput.convertTo(floatImg, CV_32F);

  IplImage* cvIpl = BridgeType::ITKImageToIplImage<ImageType>(inputImage);
  printf("row:%d col:%d\n", cvInput.rows, cvInput.cols);
  for(int i = 0; i < cvInput.rows; i++) {
    for(int j = 0; j < cvInput.cols; j++) {
      printf("%d ",cvInput.data[i * cvInput.cols+ j]);
#if 0
      unsigned char temp = cvInput.data[i * cvInput.cols + j];
      printf("te:%d ", temp);
      cvFloatInput.data[i * cvInput.cols + j] = 11.11; //static_cast<float>(temp);
        //(float)cvInput.data[i * cvInput.cols + j];
      printf("%f ", cvFloatInput.data[i * cvInput.cols+ j]);
      std::cout << floatImg.data[i * cvInput.cols+ j] << " ";
#endif
    }
    printf("\n");
  }
  printf("depth:%d\n", cvInput.depth());
  for(int i = 0; i < cvIpl->height; i++) {
    for(int j = 0; j < cvIpl->width; j++) {
      printf("%d ",cvIpl->imageData[i * cvIpl->widthStep]);
    }
    printf("\n");
  }
#if 0
  cv::Mat input = cv::imread("./BiasNonCorrect32x28.png", 0);
  cv::imshow("test", input);
  cv::waitKey(0);
  CvImageType::Pointer itkImage = BridgeType::CVMatToITKImage<CvImageType>(input);
  CvFilterType::Pointer cvFilter = CvFilterType::New();
  CvImageType::SizeType radius;
  radius[0] = 9;
  radius[1] = 9;
  cvFilter->SetRadius(radius);

  cvFilter->SetInput(itkImage);
  cvFilter->Update();
  cv::Mat resultImage = BridgeType::ITKImageToCVMat<CvImageType>(cvFilter->GetOutput());
  printf("row:%d col:%d\n", resultImage.rows, resultImage.cols);
  for(int i = 0; i < resultImage.rows; i++) {
    for(int j = 0; j < resultImage.cols; j++) {
      printf("%d ",resultImage.data[i * resultImage.cols+ j]);
    }
    printf("\n");
  }

  cv::imshow("test", input);
  cv::waitKey(0);
  return 0;
#endif
#endif
  // handle the mask image
  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage = NULL;

  if( argc > 6 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( argv[6] );
    try
      {
      maskreader->Update();
      maskImage = maskreader->GetOutput();
      maskImage->DisconnectPipeline();
      }
    catch( ... )
      {
      maskImage = NULL;
      }
    }

  if( !maskImage )
    {
    std::cout << "Mask not read.  Creating Otsu mask." << std::endl;
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType>
    ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput( inputImage );
    // otsu->SetNumberOfHistogramBins( 200 );
    otsu->SetInsideValue( 0 );
    //otsu->SetOutsideValue( 255 );
    otsu->SetOutsideValue( 1);

    otsu->Update();
    maskImage = otsu->GetOutput();
    maskImage->DisconnectPipeline();
    }

#define WRITE_MASK_IMAGE 1
#ifdef WRITE_MASK_IMAGE
  typedef itk::ImageFileWriter<MaskImageType> MaskImgWriterType;
  typename MaskImgWriterType::Pointer maskWriter = MaskImgWriterType::New();
  std::string maskFileName = (std::string)argv[2] + ".mask.png";
  maskWriter->SetFileName(maskFileName);
  maskWriter->SetInput(maskImage);
  maskWriter->Update();
#endif



  // instantiate N4 and assign variables not exposed to the user in this test.
  typedef itk::N4BiasFieldCorrectionImageFilterOrig<ImageType, MaskImageType,
                                                ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetMaskLabel( 1 );
  correcter->SetSplineOrder( 3 );
  correcter->SetWienerFilterNoise( 0.01 );
  correcter->SetBiasFieldFullWidthAtHalfMaximum( 0.15 );
  correcter->SetConvergenceThreshold( 0.0000001 );

  // handle the number of iterations
  std::vector<unsigned int> numIters = ConvertVector<unsigned int>(
      std::string( "100x50x50" ) );
  if( argc > 5 )
    {
    numIters = ConvertVector<unsigned int>( argv[5] );
    }
#define MAX_ITE_SIZE 10
  unsigned int maximumNumberOfIterations[MAX_ITE_SIZE];
  unsigned int iteSize = numIters.size();
  for( unsigned int d = 0; d < numIters.size(); d++ )
    {
    maximumNumberOfIterations[d] = numIters[d];
    }
  correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations, iteSize);

  //typename CorrecterType::ArrayType numberOfFittingLevels;
  ////numberOfFittingLevels.Fill( numIters.size() );
  correcter->SetNumberOfFittingLevels(iteSize);

  /* B-spline options -- we place this here to take care of the case where
   * the user wants to specify things in terms of the spline distance.
   *  1. need to pad the images to get as close to possible to the
   *     requested domain size.
   */
  typename ImageType::PointType newOrigin = inputImage->GetOrigin();

  typename CorrecterType::ArrayType numberOfControlPoints;
#define DIMENSION 2
  unsigned int numberOfControlPoints_[DIMENSION];

  float splineDistance = 200;
  if( argc > 7 )
    {
    splineDistance = atof( argv[7] );
    }

  itk::SizeValueType lowerBound[ImageDimension];
  itk::SizeValueType upperBound[ImageDimension];

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    float domain = static_cast<RealType>( inputImage->
      GetLargestPossibleRegion().GetSize()[d] - 1 ) *
      inputImage->GetSpacing()[d];
    unsigned int numberOfSpans = static_cast<unsigned int>(
      vcl_ceil( domain / splineDistance ) );
    unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans *
      splineDistance - domain ) / inputImage->GetSpacing()[d] + 0.5 );
    lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
    upperBound[d] = extraPadding - lowerBound[d];
    newOrigin[d] -= ( static_cast<RealType>( lowerBound[d] ) *
                      inputImage->GetSpacing()[d] );
    numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
    numberOfControlPoints_[d] = numberOfSpans + correcter->GetSplineOrder();
    }

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  typename PadderType::Pointer padder = PadderType::New();
  padder->SetInput( inputImage );
  padder->SetPadLowerBound( lowerBound );
  padder->SetPadUpperBound( upperBound );
  padder->SetConstant( 0 );
  padder->Update();

  inputImage = padder->GetOutput();
  inputImage->DisconnectPipeline();

  typedef itk::ConstantPadImageFilter<MaskImageType, MaskImageType>
  MaskPadderType;
  typename MaskPadderType::Pointer maskPadder = MaskPadderType::New();
  maskPadder->SetInput( maskImage );
  maskPadder->SetPadLowerBound( lowerBound );
  maskPadder->SetPadUpperBound( upperBound );
  maskPadder->SetConstant( 0 );
  maskPadder->Update();

  maskImage = maskPadder->GetOutput();
  maskImage->DisconnectPipeline();

  correcter->SetNumberOfControlPoints( numberOfControlPoints_ );

  // handle the shrink factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( inputImage );
  shrinker->SetShrinkFactors( 1 );

  typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType>
  MaskShrinkerType;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput( maskImage );
  maskshrinker->SetShrinkFactors( 1 );

  if( argc > 4 )
    {
    shrinker->SetShrinkFactors( atoi( argv[4] ) );
    maskshrinker->SetShrinkFactors( atoi( argv[4] ) );
    }
  shrinker->Update();
  inputImage = shrinker->GetOutput();
  inputImage->DisconnectPipeline();

  maskshrinker->Update();
  maskImage = maskshrinker->GetOutput();
  maskImage->DisconnectPipeline();

  // set the input image and mask image
  correcter->SetInput( inputImage );
  correcter->SetMaskImage( maskImage );

  typedef CommandIterationUpdate<CorrecterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  correcter->AddObserver( itk::IterationEvent(), observer );

#if 0
  //check input image by yao
  cv::Mat checkInput = BridgeType::ITKImageToCVMat<ImageType>(inputImage);
  int rows = checkInput.rows;
  int cols = checkInput.cols;


  printf("----------------------------------\n");
  printf("row:%d col:%d\n", rows, cols);
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      if(checkInput.data[i * cols + j] != 0.0f) {
        printf("row:%d col:%d %f ",i, j, checkInput.data[i * cols+ j]);
      }
    }
  }

  printf("----------------------------------\n");
  itk::ImageRegionIterator<ImageType> it2(inputImage, inputImage->GetRequestedRegion());
  it2.GoToBegin();
  while(!it2.IsAtEnd()) {
    std::cout << it2.Get() << std::endl;
    ++it2;
  }
  //cv::imshow("test", checkInput);
  //cv::waitKey(0);
#endif
#if 0
  printf("----------------------------------\n");
  itk::ImageScanlineConstIterator<ImageType> it2(inputImage, inputImage->GetRequestedRegion());
  it2.GoToBegin();
  while(!it2.IsAtEnd()) {
    while(!it2.IsAtEndOfLine()) {
      std::cout << it2.Get() << " ";
      ++it2;
    }
    it2.NextLine();
    std::cout << std::endl;
  }
#endif


  try
    {
    correcter->Update();
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  //for Debug.
  //return 0;

  correcter->Print( std::cout, 3 );

#define WRITE_IMAGE 1
#ifdef WRITE_IMAGE
  typedef itk::ImageFileWriter<ImageType> ImgWriterType;
  typename ImgWriterType::Pointer writer = ImgWriterType::New();
  std::string outputFileName = (std::string)argv[2] + ".out.nii.gz";
  writer->SetFileName(outputFileName);
  writer->SetInput( correcter->GetOutput() );
  writer->Update();
#endif

  // Test the reconstruction of the log bias field
  ImagePointer originalInputImage = reader->GetOutput();
  reader->UpdateOutputInformation();
  typedef itk::BSplineControlPointImageFilter
  <typename CorrecterType::BiasFieldControlPointLatticeType, typename
   CorrecterType::ScalarImageType> BSplinerType;
  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
  bspliner->SetSplineOrder( correcter->GetSplineOrder() );
  bspliner->SetSize(
    originalInputImage->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( originalInputImage->GetOrigin() );
  bspliner->SetDirection( originalInputImage->GetDirection() );
  bspliner->SetSpacing( originalInputImage->GetSpacing() );
  bspliner->Update();


#ifdef WRITE_CONTROL_POINTS
  // output the log bias field control point lattice
  typedef itk::ImageFileWriter<
  typename CorrecterType::BiasFieldControlPointLatticeType> WriterType;
  typename WriterType::Pointer imgWriter = WriterType::New();
  imgWriter->SetFileName( argv[3] );
  imgWriter->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
  imgWriter->Update();
#endif

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
//int itkN4BiasFieldCorrectionImageFilterTest( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage "
              << "outputLogControlPointLattice [shrinkFactor,default=1] "
              << "[numberOfIterations,default=100x50x50] "
              << " [maskImageWithLabelEqualTo1] [splineDistance,default=200]"
              << std::endl;
    exit( EXIT_FAILURE );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      return N4<2>( argc, argv );
      break;
    case 3:
      return N4<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
  return EXIT_SUCCESS;
}
