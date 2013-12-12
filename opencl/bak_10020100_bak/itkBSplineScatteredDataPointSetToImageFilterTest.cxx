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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilterOrig.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"

#define MASK_VALUE 128
/**
 * In this test, we approximate a 2-D scalar field.
 * The scattered data is derived from a segmented
 * image.  We write the output to an image for
 * comparison.
 */
//int itkBSplineScatteredDataPointSetToImageFilterTest( int argc, char * argv [] )
int main(int argc, char const* argv[])
{
  if ( argc != 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int ParametricDimension = 2;
  const unsigned int DataDimension = 1;

  typedef int                                           PixelType;
  typedef itk::Image<PixelType, ParametricDimension>    InputImageType;
  typedef float                                         RealType;
  typedef itk::Vector<RealType, DataDimension>          VectorType;
  typedef itk::Image<VectorType, ParametricDimension>   VectorImageType;
  typedef itk::PointSet
    <VectorImageType::PixelType, ParametricDimension>   PointSetType;

  PointSetType::Pointer pointSet = PointSetType::New();

  typedef itk::ImageFileReader<InputImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::Image<unsigned char, ParametricDimension> MaskImageType;
  MaskImageType::Pointer maskImage = NULL;
  std::cout << "Creating Otsu mask." << std::endl;
  typedef itk::OtsuThresholdImageFilter<InputImageType, MaskImageType> ThresholderType;
  ThresholderType::Pointer otsu = ThresholderType::New();
  otsu->SetInput( reader->GetOutput());
  // otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetInsideValue( 0 );
  //otsu->SetOutsideValue( 255 );
  otsu->SetOutsideValue( MASK_VALUE);

  otsu->Update();
  maskImage = otsu->GetOutput();
  maskImage->DisconnectPipeline();

#define WRITE_MASK_IMAGE 1
#ifdef WRITE_MASK_IMAGE
  typedef itk::ImageFileWriter<MaskImageType> MaskImgWriterType;
  MaskImgWriterType::Pointer maskWriter = MaskImgWriterType::New();
  std::string maskFileName = (std::string)argv[2] + ".mask.png";
  maskWriter->SetFileName(maskFileName);
  maskWriter->SetInput(maskImage);
  maskWriter->Update();
#endif



  itk::ImageRegionIteratorWithIndex<InputImageType>
    It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

  // Iterate through the input image which consists of multivalued
  // foreground pixels (=nonzero) and background values (=zero).
  // The foreground pixels comprise the input point set.

  unsigned int pointCount = 0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if (maskImage->GetPixel( It.GetIndex()) == MASK_VALUE)
    //if ( It.Get() != itk::NumericTraits<PixelType>::Zero )
      {
        pointCount++;
      }
    }

  DataSet_* inputData = new DataSet_[pointCount];
  int j = 0;

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if (maskImage->GetPixel( It.GetIndex()) == MASK_VALUE)
    //if ( It.Get() != itk::NumericTraits<PixelType>::Zero )
      {
   
      // We extract both the 2-D location of the point
      // and the pixel value of that point.

      PointSetType::PointType point;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

      unsigned long i = pointSet->GetNumberOfPoints();
      pointSet->SetPoint( i, point );

      PointSetType::PixelType V( DataDimension );
      V[0] = static_cast<RealType>( It.Get() );
      pointSet->SetPointData( i, V );

      inputData[j].index[0] = point[0];
      inputData[j].index[1] = point[1];
      inputData[j].data[0] = static_cast<RealType>(It.Get());
      j++;
      }
    }


  // Instantiate the B-spline filter and set the desired parameters.
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <PointSetType, VectorImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetSplineOrder( 3 );
  unsigned int* ncps = new unsigned int[DIM];
  for(int i = 0; i < DIM; i++) {
    ncps[i] = 8;
  }
  //ncps.Fill( 8 );
  filter->SetNumberOfControlPoints( ncps );

  //TODO:at first, I ignore muliti level supports.
  //filter->SetNumberOfLevels( 3 );

  unsigned int* close = new unsigned int[DIM];
  for(int i = 0; i < DIM; i++) {
    close[i] = 0;
  }
  //close.Fill( 0 );
  filter->SetCloseDimension( close );

  // Define the parametric domain.
  filter->SetOrigin( reader->GetOutput()->GetOrigin() );
  filter->SetSpacing( reader->GetOutput()->GetSpacing() );
  filter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  filter->SetDirection( reader->GetOutput()->GetDirection() );

  filter->SetInput( pointSet );

  /* original codes by yao. */
  filter->AllocDataSet(pointCount);
  filter->SetDataSet(inputData);

  /* execute */
  try
    {
    filter->Update();
    }
  catch (...)
    {
    std::cerr << "Test 1: itkBSplineScatteredDataImageFilter exception thrown"
              << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Origin: " << filter->GetOrigin() << std::endl;
  std::cout << "Spacing: " << filter->GetSpacing() << std::endl;
  std::cout << "Size: " << filter->GetSize() << std::endl;
  std::cout << "Direction: " << filter->GetDirection() << std::endl;

#if 0
  std::cout << "Number of control points: " <<
    filter->GetNumberOfControlPoints() << std::endl;
  std::cout << "Current number of control points: " <<
    filter->GetCurrentNumberOfControlPoints() << std::endl;
  std::cout << "Number of levels: " <<
    filter->GetNumberOfLevels() << std::endl;
  std::cout << "Close dimension: " <<
    filter->GetCloseDimension() << std::endl;
  std::cout << "Spline order: " << filter->GetSplineOrder() << std::endl;
#endif

  typedef itk::Image<RealType, ParametricDimension> RealImageType;
#if 0
  VectorImageType *outputImage = filter->GetOutput();
  // Write the output to an image.
  RealImageType::Pointer image = RealImageType::New();
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image->Allocate();
  itk::ImageRegionIteratorWithIndex<RealImageType>
    Itt( image, image->GetLargestPossibleRegion() );

  for ( Itt.GoToBegin(); !Itt.IsAtEnd(); ++Itt )
    {
    Itt.Set( outputImage->GetPixel( Itt.GetIndex() )[0] );
    //std::cout << "output" <<outputImage->GetPixel(Itt.GetIndex()) << std::endl;
    //std::cout << "phiLattice" <<phiLattice->GetPixel(Itt.GetIndex()) << std::endl;
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( argv[2] );
  writer->Update();
#endif

#if 1
  VectorImageType *phiLattice = filter->GetPhiLattice();

  typedef itk::BSplineControlPointImageFilter
    <FilterType::PointDataImageType, VectorImageType> BSplineReconstructerType;
  BSplineReconstructerType::Pointer reconstructer =
    BSplineReconstructerType::New();
  reconstructer->SetInput(phiLattice);
  reconstructer->SetOrigin( reader->GetOutput()->GetOrigin() );
  reconstructer->SetSpacing( reader->GetOutput()->GetSpacing() );
  reconstructer->SetDirection( reader->GetOutput()->GetDirection() );
  reconstructer->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  reconstructer->Update();
  VectorImageType* out = reconstructer->GetOutput();

  // Write the output to an image.
  RealImageType::Pointer image_ = RealImageType::New();
  image_->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image_->Allocate();
  itk::ImageRegionIteratorWithIndex<RealImageType>
    Itt_( image_, image_->GetLargestPossibleRegion() );

  for ( Itt_.GoToBegin(); !Itt_.IsAtEnd(); ++Itt_ )
    {
    Itt_.Set( out->GetPixel( Itt_.GetIndex() )[0] );
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer_ = WriterType::New();
  writer_->SetInput( image_ );
  writer_->SetFileName(argv[2]);
  writer_->Update();


#endif

  return EXIT_SUCCESS;
}
