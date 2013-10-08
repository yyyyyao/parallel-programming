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
#ifndef __itkBSplineScatteredDataPointSetToImageFilter_hxx
#define __itkBSplineScatteredDataPointSetToImageFilter_hxx

#include "itkBSplineScatteredDataPointSetToImageFilterOrig.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"
#include "itkCastImageFilter.h"
#include "itkNumericTraits.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_vector.h"
#include "vcl_limits.h"

#include <stdio.h>
namespace itk
{
/**
 * \author Nicholas J. Tustison
 *
 * This code was contributed in the Insight Journal paper:
 * "N-D C^k B-Spline Scattered Data Approximation"
 * by Nicholas J. Tustison, James C. Gee
 * http://hdl.handle.net/1926/140
 * http://www.insight-journal.org/browse/publication/57
 *
 */

template<class TInputPointSet, class TOutputImage>
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::BSplineScatteredDataPointSetToImageFilter() //constructor
{

  //malloc.
  this->m_NumberOfControlPoints = new unsigned int[ImageDimension];
  this->m_CurrentNumberOfControlPoints = new unsigned int[ImageDimension];
  this->m_CloseDimension = new unsigned int[ImageDimension];
  this->m_SplineOrder = new unsigned int[ImageDimension];
  this->m_NumberOfLevels = new unsigned int[ImageDimension];

  for(int i = 0; i < DIM; i++) {
    this->m_SplineOrder[i] = 3; //default is 3.
  }

  for( unsigned int i = 0; i < ImageDimension; i++ ) //the following value is depends on dim.
    {
    this->m_NumberOfControlPoints[i] = this->m_SplineOrder[i] + 1;
#ifdef MULTI_LEVEL
    this->m_Kernel[i] = KernelType::New();
    this->m_Kernel[i]->SetSplineOrder( this->m_SplineOrder[i] );
#endif
    }
#if 0
  this->m_KernelOrder0 = KernelOrder0Type::New();
  this->m_KernelOrder1 = KernelOrder1Type::New();
  this->m_KernelOrder2 = KernelOrder2Type::New();
  this->m_KernelOrder3 = KernelOrder3Type::New();
#endif

  for(int i = 0; i < DIM; i++) {
    this->m_CloseDimension[i] = 0; //what is this value?
  }
  this->m_DoMultilevel = false;
  //this->m_GenerateOutputImage = true;
  this->m_GenerateOutputImage = false;
  for(int i = 0; i < DIM; i++) {
    this->m_NumberOfLevels[i] = 1; //what is this value?
  }
  this->m_MaximumNumberOfLevels = 1; //what is this?

  this->m_PhiLattice = NULL;

#ifdef MULTI_LEVEL
  this->m_PsiLattice = PointDataImageType::New();
#endif
  //this->m_InputPointData = PointDataContainerType::New();
  //this->m_OutputPointData = PointDataContainerType::New();

  //this->m_PointWeights = WeightsContainerType::New(); //all of weight value is same.
  this->m_UsePointWeights = false;

  this->m_BSplineEpsilon = vcl_numeric_limits<RealType>::epsilon(); //this is the constant value.

  this->m_IsFittingComplete = false; //what is this?
}

template<class TInputPointSet, class TOutputImage>
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::~BSplineScatteredDataPointSetToImageFilter()
{}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetSplineOrder( unsigned int order )
{
  for(int i = 0; i < DIM; i++) { 
    this->m_SplineOrder[i] = order;
  }
  this->SetSplineOrder( this->m_SplineOrder );
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetSplineOrder(unsigned int* order )
{
  itkDebugMacro("Setting m_SplineOrder to " << order);

  this->m_SplineOrder = order;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_SplineOrder[i] == 0 )
      {
      itkExceptionMacro(
        "The spline order in each dimension must be greater than 0");
      }

#ifdef MULTI_LEVEL
    this->m_Kernel[i] = KernelType::New();
    this->m_Kernel[i]->SetSplineOrder( this->m_SplineOrder[i] );

    if( this->m_DoMultilevel )
      {
      typename KernelType::MatrixType C;
      C = this->m_Kernel[i]->GetShapeFunctionsInZeroToOneInterval();

      vnl_matrix<RealType> R;
      vnl_matrix<RealType> S;
      R.set_size( C.rows(), C.cols() );
      S.set_size( C.rows(), C.cols() );
      for( unsigned int j = 0; j < C.rows(); j++ )
        {
        for( unsigned int k = 0; k < C.cols(); k++ )
          {
          R(j, k) = S(j, k) = static_cast<RealType>( C(j, k) );
          }
        }
      for( unsigned int j = 0; j < C.cols(); j++ )
        {
        RealType c = vcl_pow( static_cast<RealType>( 2.0 ),
          static_cast<RealType>( C.cols() ) - j - 1 );

        for( unsigned int k = 0; k < C.rows(); k++ )
          {
          R(k, j) *= c;
          }
        }
      R = R.transpose();
      R.flipud();
      S = S.transpose();
      S.flipud();

      this->m_RefinedLatticeCoefficients[i] =
        ( vnl_svd<RealType>( R ).solve( S ) ).extract( 2, S.cols() );
      }
#endif
    }
  this->Modified();
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::AllocDataSet(unsigned int pointCount)
{
  this->m_InputPointDataSet_ = new DataSet_[pointCount];
  this->m_OutputPointDataSet_ = new DataSet_[pointCount];
  this->m_DataSet_Weight = new float[pointCount];
  this->m_PointCount = pointCount;
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetDataSet(DataSet_* input)
{
  for(int i = 0; i < m_PointCount; i++) {
    this->m_InputPointDataSet_[i].index[0] = input[i].index[0];
    this->m_InputPointDataSet_[i].index[1] = input[i].index[1];
    this->m_InputPointDataSet_[i].data[0] = input[i].data[0];

    this->m_OutputPointDataSet_[i].index[0] = input[i].index[0];
    this->m_OutputPointDataSet_[i].index[1] = input[i].index[1];
    this->m_OutputPointDataSet_[i].data[0] = input[i].data[0];

    this->m_DataSet_Weight[i] = 1;
  }
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetNumberOfLevels( unsigned int levels )
{
  for(int i = 0; i < DIM; i++) {
    this->m_NumberOfLevels[i] = levels;
  }
  this->SetNumberOfLevels( this->m_NumberOfLevels );
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetNumberOfLevels(unsigned int* levels )
{
  this->m_NumberOfLevels = levels;
  this->m_MaximumNumberOfLevels = 1;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_NumberOfLevels[i] == 0 )
      {
      itkExceptionMacro(
        "The number of levels in each dimension must be greater than 0");
      }
    if( this->m_NumberOfLevels[i] > this->m_MaximumNumberOfLevels )
      {
      this->m_MaximumNumberOfLevels = this->m_NumberOfLevels[i];
      }
    }

  itkDebugMacro( "Setting m_NumberOfLevels to "
    << this->m_NumberOfLevels );
  itkDebugMacro( "Setting m_MaximumNumberOfLevels to "
    << this->m_MaximumNumberOfLevels );

  if( this->m_MaximumNumberOfLevels > 1 )
    {
    this->m_DoMultilevel = true;
    }
  else
    {
    this->m_DoMultilevel = false;
    }
  this->SetSplineOrder( this->m_SplineOrder );
  this->Modified();
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetPointWeights(float* weights )
{
  this->m_UsePointWeights = true;
  for(int i = 0; i < m_PointCount; i++) {
    m_DataSet_Weight[i] = weights[i];
  }
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::GenerateData()
{
  TOutputImage *output = this->GetOutput();
  //const TInputPointSet *inputPointSet = this->GetInput();
  /**
   *  Create the output image
   */
  itkDebugMacro( "Size: " << this->m_Size );
  itkDebugMacro( "Origin: " << this->m_Origin );
  itkDebugMacro( "Spacing: " << this->m_Spacing );
  itkDebugMacro( "Direction: " << this->m_Direction );

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_Size[i] == 0 ) /* image size. width & height */
      {
      itkExceptionMacro("Size must be specified.");
      }
    }

  output->SetOrigin( this->m_Origin );
  output->SetSpacing( this->m_Spacing );
  output->SetDirection( this->m_Direction );
  output->SetRegions( this->m_Size );
  output->Allocate(); /* allocate output images. */

  /**
   * Perform some error checking on the input
   */
  //if( this->m_UsePointWeights &&
  //  ( this->m_PointWeights->Size() != inputPointSet->GetNumberOfPoints() ) )
  //  {
  //  itkExceptionMacro(
  //    "The number of weight points and input points must be equal." );
  //  }

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_NumberOfControlPoints[i] < this->m_SplineOrder[i] + 1 )
      {
      itkExceptionMacro(
        "The number of control points must be greater than the spline order." );
      }
    }

  /**
   * Calculate the appropriate epsilon value.
   */
  unsigned int maximumNumberOfSpans = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    unsigned int numberOfSpans = this->m_NumberOfControlPoints[d] -
      this->m_SplineOrder[d];
    numberOfSpans <<= ( this->m_NumberOfLevels[d] - 1 );
    if( numberOfSpans > maximumNumberOfSpans )
      {
      maximumNumberOfSpans = numberOfSpans;
      }
    }
  this->m_BSplineEpsilon = 100 * vcl_numeric_limits<RealType>::epsilon();
  while( static_cast<RealType>( maximumNumberOfSpans ) ==
    static_cast<RealType>( maximumNumberOfSpans ) - this->m_BSplineEpsilon )
    {
    this->m_BSplineEpsilon *= 10;
    }

#if 0
  this->m_InputPointData->Initialize();
  this->m_OutputPointData->Initialize();
  if( inputPointSet->GetNumberOfPoints() > 0 )
    {
    typename PointDataContainerType::ConstIterator It =
      inputPointSet->GetPointData()->Begin();
    while( It != inputPointSet->GetPointData()->End() ) /* Input PointSet's Iterator. */
      {
      if( !this->m_UsePointWeights )
        {
        this->m_PointWeights->InsertElement( It.Index(), 1.0 );
        }
      this->m_InputPointData->InsertElement( It.Index(), It.Value() );
      this->m_OutputPointData->InsertElement( It.Index(), It.Value() ); /* output point data also has input? */
      ++It;
      }
    }
#endif

  this->m_CurrentLevel = 0;
  this->m_CurrentNumberOfControlPoints = this->m_NumberOfControlPoints; /* at first, set NumberofCtl Pts these are array. */

  /**
   * Set up multithread processing to handle generating the
   * control point lattice.
   */
  typename ImageSource<TOutputImage>::ThreadStruct str1;
  str1.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str1 );

  /**
   * Multithread the generation of the control point lattice.
   */
  this->BeforeThreadedGenerateData();
  this->GetMultiThreader()->SingleMethodExecute();
  this->AfterThreadedGenerateData();

  this->UpdatePointSet();

#ifdef MULTI_LEVEL
  if( this->m_DoMultilevel )
    {
    this->m_PsiLattice->SetRegions(
      this->m_PhiLattice->GetLargestPossibleRegion() );
    this->m_PsiLattice->Allocate();
    PointDataType P(0.0);
    this->m_PsiLattice->FillBuffer(P);
    }
#endif

#ifdef MULTI_LEVEL
  for( this->m_CurrentLevel = 1; /* maybe, this loop is the biggest part! */
        this->m_CurrentLevel < this->m_MaximumNumberOfLevels;
        this->m_CurrentLevel++ )
    {
    ImageRegionIterator<PointDataImageType> ItPsi(
      this->m_PsiLattice, this->m_PsiLattice->GetLargestPossibleRegion() );
    ImageRegionIterator<PointDataImageType> ItPhi(
      this->m_PhiLattice, this->m_PhiLattice->GetLargestPossibleRegion() );
    for( ItPsi.GoToBegin(), ItPhi.GoToBegin(); !ItPsi.IsAtEnd();
      ++ItPsi, ++ItPhi )
      {
      ItPsi.Set( ItPhi.Get() + ItPsi.Get() );
      }
    this->RefineControlPointLattice();

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
        {
        this->m_CurrentNumberOfControlPoints[i] =
          2 * this->m_CurrentNumberOfControlPoints[i] - this->m_SplineOrder[i];
        }
      }

    itkDebugMacro( "Current Level = " << this->m_CurrentLevel );
    itkDebugMacro( "  Current number of control points = "
      << this->m_CurrentNumberOfControlPoints );

    RealType averageDifference = 0.0;
    RealType totalWeight = 0.0;

#if 0
    typename PointDataContainerType::Iterator ItIn =
      this->m_InputPointData->Begin();
    typename PointDataContainerType::Iterator ItOut =
      this->m_OutputPointData->Begin();
    while( ItIn != this->m_InputPointData->End() )
      {
      this->m_InputPointData->InsertElement(
        ItIn.Index(), ItIn.Value() - ItOut.Value() );

      if( this->GetDebug() )
        {
        RealType weight = this->m_PointWeights->GetElement( ItIn.Index() );
        averageDifference += ( ItIn.Value() - ItOut.Value() ).GetNorm() * weight;
        totalWeight += weight;
        }

      ++ItIn;
      ++ItOut;
      }
#endif
    for(int i = 0; i < m_PointCount; i++) {
      m_InputPointDataSet_[i].data[0] -= m_OutputPointDataSet_[i].data[0];
    }
    /**
     * Set up multithread processing to handle generating the
     * control point lattice.
     */
    typename ImageSource<ImageType>::ThreadStruct str2;
    str2.Filter = this;

    this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
    this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str2 );

    /**
     * Multithread the generation of the control point lattice.
     */
    this->BeforeThreadedGenerateData();
    this->GetMultiThreader()->SingleMethodExecute();
    this->AfterThreadedGenerateData();

    this->UpdatePointSet();
    }
#endif

#ifdef MULTI_LEVEL
  if( this->m_DoMultilevel ) //at first, I ignore.
    {
    ImageRegionIterator<PointDataImageType> ItPsi( this->m_PsiLattice,
      this->m_PsiLattice->GetLargestPossibleRegion() );
    ImageRegionIterator<PointDataImageType> ItPhi( this->m_PhiLattice,
      this->m_PhiLattice->GetLargestPossibleRegion() );
    for( ItPsi.GoToBegin(), ItPhi.GoToBegin(); !ItPsi.IsAtEnd();
      ++ItPsi, ++ItPhi )
      {
      ItPsi.Set( ItPhi.Get() + ItPsi.Get() );
      }

    typedef ImageDuplicator<PointDataImageType> ImageDuplicatorType;
    typename ImageDuplicatorType::Pointer duplicator =
      ImageDuplicatorType::New();
    duplicator->SetInputImage( this->m_PsiLattice );
    duplicator->Update();
    this->m_PhiLattice = duplicator->GetModifiableOutput();

    this->UpdatePointSet();
    }
#endif

  this->m_IsFittingComplete = true;

  if( this->m_GenerateOutputImage ) //I can ignore this codes.
    {
    typename ImageSource<ImageType>::ThreadStruct str3;
    str3.Filter = this;

    this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
    this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str3 );

//    this->BeforeThreadedGenerateData();
    this->GetMultiThreader()->SingleMethodExecute();
//    this->AfterThreadedGenerateData();
    }

  this->SetPhiLatticeParametricDomainParameters();
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::BeforeThreadedGenerateData()
{
  if( !this->m_IsFittingComplete )
    {
    m_ThreadNum = this->GetNumberOfThreads();
    this->m_rawDeltaLattice = new float*[m_ThreadNum];
    this->m_rawOmegaLattice = new float*[m_ThreadNum];

    this->m_DeltaLatticePerThread.resize( this->GetNumberOfThreads() );
    this->m_OmegaLatticePerThread.resize( this->GetNumberOfThreads() );

    typename RealImageType::SizeType size;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( this->m_CloseDimension[i] )
        {
        size[i] = this->m_CurrentNumberOfControlPoints[i] - this->m_SplineOrder[i];
        }
      else
        {
        size[i] = this->m_CurrentNumberOfControlPoints[i];
        }
      }

    for( unsigned int n = 0; n < this->GetNumberOfThreads(); n++ )
      {
      this->m_OmegaLatticePerThread[n] = RealImageType::New();
      this->m_OmegaLatticePerThread[n]->SetRegions( size );
      this->m_OmegaLatticePerThread[n]->Allocate();
      this->m_OmegaLatticePerThread[n]->FillBuffer( 0.0 );

      this->m_DeltaLatticePerThread[n] = PointDataImageType::New();
      this->m_DeltaLatticePerThread[n]->SetRegions( size );
      this->m_DeltaLatticePerThread[n]->Allocate();
      this->m_DeltaLatticePerThread[n]->FillBuffer( 0.0 );

      m_oneDSize = size[0] * size[1];
      this->m_rawDeltaLattice[n] = new float[m_oneDSize];
      this->m_rawOmegaLattice[n] = new float[m_oneDSize];
      for(int _k = 0; _k < m_oneDSize; _k++) {
        this->m_rawDeltaLattice[n][_k] = 0; 
        this->m_rawOmegaLattice[n][_k] = 0;
      }
      }
    }
}

template<class TInputPointSet, class TOutputImage>
unsigned int
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SplitRequestedRegion( unsigned int i, unsigned int num,
  RegionType &splitRegion )
{
  // For fitting, the image regions are not used so we always return a valid
  // number.
  if( !this->m_IsFittingComplete )
    {
    return this->GetNumberOfThreads();
    }
  else // we split on the output region for reconstruction
    {
    // Get the output pointer
    ImageType *outputPtr = this->GetOutput();

    const SizeType requestedRegionSize =
      outputPtr->GetRequestedRegion().GetSize();

    int splitAxis;
    typename TOutputImage::IndexType splitIndex;
    typename TOutputImage::SizeType splitSize;

    // Initialize the splitRegion to the output requested region
    splitRegion = outputPtr->GetRequestedRegion();
    splitIndex = splitRegion.GetIndex();
    splitSize = splitRegion.GetSize();

    // split on the outermost dimension
    splitAxis = outputPtr->GetImageDimension() - 1;

    // determine the actual number of pieces that will be generated
    typename SizeType::SizeValueType range = requestedRegionSize[splitAxis];
    unsigned int valuesPerThread = static_cast<unsigned int>( vcl_ceil(
      range / static_cast<double>( num ) ) );
    unsigned int maxThreadIdUsed = static_cast<unsigned int>( vcl_ceil(
      range / static_cast<double>( valuesPerThread ) ) - 1 );

    // Split the region
    if ( i < maxThreadIdUsed )
      {
      splitIndex[splitAxis] += i * valuesPerThread;
      splitSize[splitAxis] = valuesPerThread;
      }
    if ( i == maxThreadIdUsed )
      {
      splitIndex[splitAxis] += i * valuesPerThread;
      // last thread needs to process the "rest" dimension being split
      splitSize[splitAxis] = splitSize[splitAxis] - i * valuesPerThread;
      }

    // set the split region ivars
    splitRegion.SetIndex( splitIndex );
    splitRegion.SetSize( splitSize );

    itkDebugMacro( "Split piece: " << splitRegion );

    return maxThreadIdUsed + 1;
    }
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::ThreadedGenerateData( const RegionType &region, ThreadIdType threadId )
{
  if( !this->m_IsFittingComplete )
    {
    this->ThreadedGenerateDataForFitting( region, threadId );
    }
#if 0 //this function is implemented by BSplineControlPointImageFilter.
  else
    {
    this->ThreadedGenerateDataForReconstruction( region, threadId );
    }
#endif
}


template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::ThreadedGenerateDataForFitting(
  const RegionType & itkNotUsed( region ), ThreadIdType threadId )
{
  //const TInputPointSet *input = this->GetInput();
  /**
   * Ignore the output region as we're only interested in dividing the
   * points among the threads.
   */
  typename RealImageType::SizeType size;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = this->m_SplineOrder[i] + 1;
    }
  RealImagePointer neighborhoodWeightImage = RealImageType::New();
  neighborhoodWeightImage->SetRegions( size ); //make m_SplineOrder[i] + 1 images. what's the meaning of this size?
  neighborhoodWeightImage->Allocate();
  neighborhoodWeightImage->FillBuffer( 0.0 );

  ImageRegionIteratorWithIndex<RealImageType> ItW(
    neighborhoodWeightImage, neighborhoodWeightImage->GetRequestedRegion() );

  float* neighborRawImage = new float[size[0] * size[1]];
  float p[ ImageDimension ]; /* vector size is 2. */
  float r[ ImageDimension ]; /* vector size is 2. */
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    r[i] = static_cast<RealType>( this->m_CurrentNumberOfControlPoints[i] -
      this->m_SplineOrder[i] ) / ( static_cast<RealType>( this->m_Size[i] - 1 ) *
      this->m_Spacing[i] );
    }

  /**
   * Determine which points should be handled by this particular thread.
   */
  ThreadIdType numberOfThreads = this->GetNumberOfThreads();
  SizeValueType numberOfPointsPerThread = static_cast<SizeValueType>(
    //input->GetNumberOfPoints() / numberOfThreads );
    m_PointCount / numberOfThreads );

  unsigned int start = threadId * numberOfPointsPerThread;
  unsigned int end = start + numberOfPointsPerThread;
  if( threadId == this->GetNumberOfThreads() - 1 )
    {
    //end = input->GetNumberOfPoints();
    end = m_PointCount;
    }

  for( unsigned int n = start; n < end; n++ ) /* executes start to end */
    {
    PointType point;
    point.Fill( 0.0 );

    //input->GetPoint( n, &point );
    for(int i = 0; i < ImageDimension; i++) {
      point[i] = m_InputPointDataSet_[n].index[i];
    }

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      unsigned int totalNumberOfSpans =
        this->m_CurrentNumberOfControlPoints[i] - this->m_SplineOrder[i];

      p[i] = ( point[i] - this->m_Origin[i] ) * r[i];
      if( vnl_math_abs( p[i] - static_cast<RealType>( totalNumberOfSpans ) ) <=
        this->m_BSplineEpsilon )
        {
        p[i] = static_cast<RealType>( totalNumberOfSpans )
               - this->m_BSplineEpsilon;
        }
      if( p[i] >= static_cast<RealType>( totalNumberOfSpans ) )
        {
        itkExceptionMacro( "The reparameterized point component " << p[i]
          << " is outside the corresponding parametric domain of [0, "
          << totalNumberOfSpans << "]." );
        }
      }

    RealType w2Sum = 0.0;

//#define DEB 1
#ifdef DEB
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW ) /* ItW is m_SplineOrder[i] + 1 * dim? */
      {
      typename RealImageType::IndexType idx = ItW.GetIndex();
#else
    for(int i_ = 0; i_ < size[0]; i_++) {
      int idx[2];
      for(int j_ = 0; j_ < size[1]; j_++) {
      idx[0] = j_;
      idx[1] = i_;
#endif
      RealType B = 1.0;
      //printf("start:%d %d %d\n", n, idx[0], idx[1]);
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        RealType u = static_cast<RealType>( p[i] -
          static_cast<unsigned>( p[i] ) - idx[i] ) + 0.5 *
          static_cast<RealType>( this->m_SplineOrder[i] - 1 );

        switch( this->m_SplineOrder[i] )
          {
          case 0:
            {
            //B *= this->m_KernelOrder0->Evaluate( u );
            B *= spline0(u);
            break;
            }
          case 1:
            {
            //B *= this->m_KernelOrder1->Evaluate( u );
            B *= spline1(u);
            break;
            }
          case 2:
            {
            //B *= this->m_KernelOrder2->Evaluate( u );
            B *= spline2(u);
            break;
            }
          case 3:
            {
            //B *= this->m_KernelOrder3->Evaluate( u );
            B *= spline3(u);
            break;
            }
          default:
            {
            //B *= this->m_Kernel[i]->Evaluate( u );
            B *= splineError(u);
            break;
            }
          }
        }

#ifdef DEB
      ItW.Set( B );
      w2Sum += B * B;
      }
#else
      neighborRawImage[idx[0] * size[1] + idx[1]] = B;
      w2Sum += B * B;
      //printf("idx[0]:%d idx[1]:%d B:%f\n", idx[0], idx[1], B);
      }
    }
#endif

    RealImageType * currentThreadOmegaLattice = this->m_OmegaLatticePerThread[threadId];
    PointDataImageType * currentThreadDeltaLattice = this->m_DeltaLatticePerThread[threadId];

    float* currentThreadRawOmegaLattice = this->m_rawOmegaLattice[threadId];
    float* currentThreadRawDeltaLattice = this->m_rawDeltaLattice[threadId];
    float data_;

#define DEB1
#ifdef DEB1
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW ) {
      typename RealImageType::IndexType idx = ItW.GetIndex();
#else

    for(int i_ = 0; i_ < size[0]; i_++) {
      int idx[2];
      for(int j_ = 0; j_ < size[1]; j_++) {
      idx[0] = j_;
      idx[1] = i_;
#endif
      int j_ = idx[0];
      int i_ = idx[1];
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        idx[i] += static_cast<unsigned>( p[i] );

        if( this->m_CloseDimension[i] )
          {
          idx[i] %= size[i];
          }

        }
      //RealType wc = this->m_PointWeights->GetElement(n);
      RealType wc = this->m_DataSet_Weight[n]; // = 1.0.

      //RealType t = ItW.Get();
      RealType t = neighborRawImage[j_ * size[1] + i_];
      //printf("idx[0]:%d idx[1]:%d ref:%f yao:%f\n", idx[0], idx[1], t, t_);
      RealType _current = currentThreadOmegaLattice->GetPixel( idx );
      _current = _current + wc * t * t;
      currentThreadOmegaLattice->SetPixel( idx,
        _current);
        //currentThreadOmegaLattice->GetPixel( idx ) + wc * t * t );
      //PointDataType data = this->m_InputPointData->GetElement( n );
      PointDataType data = this->m_InputPointDataSet_[n].data;
      data_ = this->m_InputPointDataSet_[n].data[0];
      data_ *= ( t * t * t * wc / w2Sum );
      data *= ( t * t * t * wc / w2Sum );
      currentThreadDeltaLattice->SetPixel( idx,
        currentThreadDeltaLattice->GetPixel( idx ) + data );

      int _idx = idx[0] * size[0] + idx[1];
      currentThreadRawOmegaLattice[_idx] = _current; 
      currentThreadRawDeltaLattice[_idx] += data_;
      }
    }
#ifdef DEB1

#else
    }
#endif
}

#if 0
template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::ThreadedGenerateDataForReconstruction( const RegionType &region, ThreadIdType
  itkNotUsed( threadId ) )
{
  typename PointDataImageType::Pointer collapsedPhiLattices[ImageDimension + 1];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    collapsedPhiLattices[i] = PointDataImageType::New();
    collapsedPhiLattices[i]->CopyInformation( this->m_PhiLattice );

    typename PointDataImageType::SizeType size;
    size.Fill(1);
    for( unsigned int j = 0; j < i; j++ )
      {
      size[j] = this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[j];
      }
    collapsedPhiLattices[i]->SetRegions( size );
    collapsedPhiLattices[i]->Allocate();
    }
  typedef ImageDuplicator<PointDataImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage( this->m_PhiLattice );
  duplicator->Update();

  collapsedPhiLattices[ImageDimension] = duplicator->GetModifiableOutput();

  unsigned int totalNumberOfSpans[DIM];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_CloseDimension[i] )
      {
      totalNumberOfSpans[i] =
        this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i];
      }
    else
      {
      totalNumberOfSpans[i] =
        this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i] -
        this->m_SplineOrder[i];
      }
    }
  FixedArray<RealType, ImageDimension> U;
  FixedArray<RealType, ImageDimension> currentU;
  currentU.Fill( -1 );

  typename ImageType::IndexType startIndex =
    this->GetOutput()->GetRequestedRegion().GetIndex();
  typename PointDataImageType::IndexType startPhiIndex =
    this->m_PhiLattice->GetLargestPossibleRegion().GetIndex();

  ImageRegionIteratorWithIndex<ImageType> It( this->GetOutput(), region );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType idx = It.GetIndex();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      U[i] = static_cast<RealType>( totalNumberOfSpans[i] ) *
        static_cast<RealType>( idx[i] - startIndex[i] ) /
        static_cast<RealType>( this->m_Size[i] - 1 );
      if( vnl_math_abs( U[i] - static_cast<RealType>( totalNumberOfSpans[i] ) )
        <= this->m_BSplineEpsilon )
        {
        U[i] = static_cast<RealType>( totalNumberOfSpans[i] ) -
          this->m_BSplineEpsilon;
        }
      if( U[i] >= static_cast<RealType>( totalNumberOfSpans[i] ) )
        {
        itkExceptionMacro( "The collapse point component " << U[i]
          << " is outside the corresponding parametric domain of [0, "
          << totalNumberOfSpans[i] << "]." );
        }
      }
    for( int i = ImageDimension - 1; i >= 0; i-- )
      {
      if( U[i] != currentU[i] )
        {
        for( int j = i; j >= 0; j-- )
          {
          this->CollapsePhiLattice( collapsedPhiLattices[j + 1],
            collapsedPhiLattices[j], U[j], j );
          currentU[j] = U[j];
          }
        break;
        }
      }
    It.Set( collapsedPhiLattices[0]->GetPixel( startPhiIndex ) ); //set output ImageData!
    }
}
#endif

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::AfterThreadedGenerateData()
{
  if( !this->m_IsFittingComplete )
    {
    /**
     * Accumulate all the delta lattice and omega lattice values to
     * calculate the final phi lattice.
     */
    ImageRegionIterator< PointDataImageType > ItD(
      this->m_DeltaLatticePerThread[0],
      this->m_DeltaLatticePerThread[0]->GetLargestPossibleRegion() );
    ImageRegionIterator< RealImageType > ItO(
      this->m_OmegaLatticePerThread[0],
      this->m_OmegaLatticePerThread[0]->GetLargestPossibleRegion() );

    for( ThreadIdType n = 1; n < this->GetNumberOfThreads(); n++ )
      {
      ImageRegionIterator< PointDataImageType > Itd(
        this->m_DeltaLatticePerThread[n],
        this->m_DeltaLatticePerThread[n]->GetLargestPossibleRegion() );
      ImageRegionIterator< RealImageType > Ito(
        this->m_OmegaLatticePerThread[n],
        this->m_OmegaLatticePerThread[n]->GetLargestPossibleRegion() );

      ItD.GoToBegin();
      ItO.GoToBegin();
      Itd.GoToBegin();
      Ito.GoToBegin();
      while( !ItD.IsAtEnd() )
        {
        ItD.Set( ItD.Get() + Itd.Get() );
        ItO.Set( ItO.Get() + Ito.Get() );

        ++ItD;
        ++ItO;
        ++Itd;
        ++Ito;
        }
      }

    float tempO[64];
    float tempD[64];
    float tempD_[64];
    ItD.GoToBegin();
    ItO.GoToBegin();
    int _j = 0;
    while( !ItD.IsAtEnd() ) {
      tempD[_j] = ItD.Get()[0];
      tempD_[_j] = ItD.Get()[1];
      tempO[_j] = ItO.Get();
      _j++;
      ++ItD;
      ++ItO;
    }

    for(int n = 1; n < m_ThreadNum; n++) {
      for(int j = 0; j < m_oneDSize; j++) {
        m_rawDeltaLattice[0][j] += m_rawDeltaLattice[n][j];
        m_rawOmegaLattice[0][j] += m_rawOmegaLattice[n][j];
      }
    }

    /**
     * Generate the control point lattice
     */
    typename RealImageType::SizeType size;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( this->m_CloseDimension[i] )
        {
        size[i] = this->m_CurrentNumberOfControlPoints[i] - this->m_SplineOrder[i];
        }
      else
        {
        size[i] = this->m_CurrentNumberOfControlPoints[i];
        }
      }
    this->m_PhiLattice = PointDataImageType::New();
    this->m_PhiLattice->SetRegions( size );
    this->m_PhiLattice->Allocate();
    this->m_PhiLattice->FillBuffer( 0.0 );
    this->m_IplPhiLattice = cvCreateImage(cvSize(size[0], size[1]),
          IPL_DEPTH_32F, 1);


    ImageRegionIterator<PointDataImageType> ItP(
      this->m_PhiLattice, this->m_PhiLattice->GetLargestPossibleRegion() );

    int _i = 0;
    float _P = 0;
    for( ItP.GoToBegin(), ItO.GoToBegin(), ItD.GoToBegin(); !ItP.IsAtEnd();
      ++ItP, ++ItO, ++ItD )
      {
      PointDataType P;
      P.Fill( 0 );
      float temp = ItO.Get();
      float temp_ = m_rawOmegaLattice[0][_i];
      if(temp != 0)
      //if( ItO.Get() != 0 )
        {
        P = ItD.Get() / ItO.Get();
        //P = m_rawDeltaLattice[0][_i] / m_rawOmegaLattice[0][_i];
        for( unsigned int i = 0; i < P.Size(); i++ )
          {
#if 0
          if( vnl_math_isnan( P[i] ) || vnl_math_isinf( P[i] ) )
            {
            P[i] = 0;
            }
#endif
          }
        ItP.Set( P );
        }

      if(m_rawOmegaLattice[_i] != 0) {
        _P = m_rawDeltaLattice[0][_i] / m_rawOmegaLattice[0][_i];

      }
      _i++;
      }
    }
}

#ifdef MULTI_LEVEL
template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::RefineControlPointLattice()
{
  unsigned int NumberOfNewControlPoints[ImageDimension];
  for(int i = 0; i < ImageDimension; i++) {
    NumberOfNewControlPoints[i] = this->m_CurrentNumberOfControlPoints[i];
  }
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
      {
      NumberOfNewControlPoints[i] = 2 *
        NumberOfNewControlPoints[i] - this->m_SplineOrder[i];
      }
    }
  typename RealImageType::SizeType size;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_CloseDimension[i] )
      {
      size[i] = NumberOfNewControlPoints[i] - this->m_SplineOrder[i];
      }
    else
      {
      size[i] = NumberOfNewControlPoints[i];
      }
    }

  PointDataImagePointer refinedLattice = PointDataImageType::New();
  refinedLattice->SetRegions( size );
  refinedLattice->Allocate();

  PointDataType data;
  data.Fill( 0.0 );
  refinedLattice->FillBuffer( data );

  typename PointDataImageType::IndexType idx;
  typename PointDataImageType::IndexType idxPsi;
  typename PointDataImageType::IndexType tmp;
  typename PointDataImageType::IndexType tmpPsi;
  typename PointDataImageType::IndexType off;
  typename PointDataImageType::IndexType offPsi;
  typename PointDataImageType::RegionType::SizeType sizePsi;

  size.Fill(2);
  unsigned int N = 1;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    N *= ( this->m_SplineOrder[i] + 1 );
    sizePsi[i] = this->m_SplineOrder[i] + 1;
    }

  ImageRegionIteratorWithIndex< PointDataImageType >
  It( refinedLattice, refinedLattice->GetLargestPossibleRegion() );

  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    idx = It.GetIndex();
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if( this->m_CurrentLevel < this->m_NumberOfLevels[i] )
        {
        idxPsi[i] = static_cast<unsigned int>( 0.5 * idx[i] );
        }
      else
        {
        idxPsi[i] = static_cast<unsigned int>( idx[i] );
        }
      }

    for( unsigned int i = 0; i < ( 2 << ( ImageDimension - 1 ) ); i++ )
      {
      PointDataType sum( 0.0 );
      PointDataType val( 0.0 );
      off = this->NumberToIndex( i, size );

      bool outOfBoundary = false;
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        tmp[j] = idx[j] + off[j];
        if( tmp[j] >= static_cast<int>( NumberOfNewControlPoints[j] ) &&
          !this->m_CloseDimension[j] )
          {
          outOfBoundary = true;
          break;
          }
        if( this->m_CloseDimension[j] )
          {
          tmp[j] %= refinedLattice->GetLargestPossibleRegion().GetSize()[j];
          }
        }
      if( outOfBoundary )
        {
        continue;
        }

      for( unsigned int j = 0; j < N; j++ )
        {
        offPsi = this->NumberToIndex( j, sizePsi );

        bool isOutOfBoundary = false;
        for( unsigned int k = 0; k < ImageDimension; k++ )
          {
          tmpPsi[k] = idxPsi[k] + offPsi[k];
          if( tmpPsi[k] >=
            static_cast<int>( this->m_CurrentNumberOfControlPoints[k] ) &&
            !this->m_CloseDimension[k] )
            {
            isOutOfBoundary = true;
            break;
            }
          if( this->m_CloseDimension[k] )
            {
            tmpPsi[k] %=
              this->m_PsiLattice->GetLargestPossibleRegion().GetSize()[k];
            }
          }
        if( isOutOfBoundary )
          {
          continue;
          }
        RealType coeff = 1.0;
        for( unsigned int k = 0; k < ImageDimension; k++ )
          {
          coeff *= this->m_RefinedLatticeCoefficients[k](off[k], offPsi[k]);
          }
        val = this->m_PsiLattice->GetPixel( tmpPsi );
        val *= coeff;
        sum += val;
        }
      refinedLattice->SetPixel( tmp, sum );
      }

    bool isEvenIndex = false;
    while( !isEvenIndex && !It.IsAtEnd() )
      {
      ++It;
      idx = It.GetIndex();
      isEvenIndex = true;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        if( idx[i] % 2 )
          {
          isEvenIndex = false;
          }
        }
      }
    }

  typedef ImageDuplicator<PointDataImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage( refinedLattice );
  duplicator->Update();
  this->m_PsiLattice = duplicator->GetModifiableOutput();
}
#endif

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::UpdatePointSet()
{
  //const TInputPointSet *input = this->GetInput();
  PointDataImagePointer collapsedPhiLattices[ImageDimension + 1];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    collapsedPhiLattices[i] = PointDataImageType::New();
    collapsedPhiLattices[i]->SetOrigin( this->m_PhiLattice->GetOrigin() );
    collapsedPhiLattices[i]->SetSpacing( this->m_PhiLattice->GetSpacing() );
    collapsedPhiLattices[i]->SetDirection( this->m_PhiLattice->GetDirection() );

    typename PointDataImageType::SizeType size;
    size.Fill( 1 );
    for( unsigned int j = 0; j < i; j++ )
      {
      size[j] = this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[j];
      }
    collapsedPhiLattices[i]->SetRegions( size );
    collapsedPhiLattices[i]->Allocate();
    }
  collapsedPhiLattices[ImageDimension] = this->m_PhiLattice;
  unsigned int totalNumberOfSpans[ImageDimension];
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if( this->m_CloseDimension[i] )
      {
      totalNumberOfSpans[i] =
        this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i];
      }
    else
      {
      totalNumberOfSpans[i] =
        this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i] -
        this->m_SplineOrder[i];
      }
    }
  FixedArray<RealType, ImageDimension> U;
  FixedArray<RealType, ImageDimension> currentU;
  currentU.Fill( -1 );

  typename PointDataImageType::IndexType startPhiIndex =
    this->m_PhiLattice->GetLargestPossibleRegion().GetIndex();

#if 0
  typename PointDataContainerType::ConstIterator ItIn =
    this->m_InputPointData->Begin();
  while( ItIn != this->m_InputPointData->End() )
#endif
  for(int j = 0; j < m_PointCount; j++)
    {
    PointType point;
    point.Fill( 0.0 );

    //input->GetPoint( ItIn.Index(), &point );
    //input->GetPoint( j, &point );
    for(int i = 0; i < ImageDimension; i++) {
      point[i] = m_InputPointDataSet_[j].index[i];
    }

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      U[i] = static_cast<RealType>( totalNumberOfSpans[i] ) *
        static_cast<RealType>( point[i] - this->m_Origin[i] ) /
        ( static_cast<RealType>( this->m_Size[i] - 1 ) * this->m_Spacing[i] );
      if( vnl_math_abs( U[i] - static_cast<RealType>( totalNumberOfSpans[i] ) )
        <= this->m_BSplineEpsilon )
        {
        U[i] = static_cast<RealType>( totalNumberOfSpans[i] ) -
          this->m_BSplineEpsilon;
        }
      if( U[i] >= static_cast<RealType>( totalNumberOfSpans[i] ) )
        {
        itkExceptionMacro( "The collapse point component " << U[i]
          << " is outside the corresponding parametric domain of [0, "
          << totalNumberOfSpans[i] << "]." );
        }
      }
    for( int i = ImageDimension - 1; i >= 0; i-- )
      {
      if( U[i] != currentU[i] )
        {
        for( int j = i; j >= 0; j-- )
          {
          this->CollapsePhiLattice( collapsedPhiLattices[j + 1],
            collapsedPhiLattices[j], U[j], j );
          currentU[j] = U[j];
          }
        break;
        }
      }
    //this->m_OutputPointData->InsertElement( ItIn.Index(),
    //this->m_OutputPointData->InsertElement( i,
    m_OutputPointDataSet_[j].data[0] = 
      collapsedPhiLattices[0]->GetPixel( startPhiIndex )[0];
    }
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::CollapsePhiLattice( PointDataImageType *lattice,
  PointDataImageType *collapsedLattice,
  const RealType u, const unsigned int dimension )
{
  ImageRegionIteratorWithIndex< PointDataImageType > It(
    collapsedLattice, collapsedLattice->GetLargestPossibleRegion() );

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PointDataType data;
    data.Fill( 0.0 );
    typename PointDataImageType::IndexType idx = It.GetIndex();
    for( unsigned int i = 0; i < this->m_SplineOrder[dimension] + 1; i++ )
      {
      idx[dimension] = static_cast<unsigned int>( u ) + i;
      RealType v = u - idx[dimension] + 0.5 * static_cast<RealType>(
        this->m_SplineOrder[dimension] - 1 );

      RealType B = 0.0;

      switch( this->m_SplineOrder[dimension] )
        {
        case 0:
          {
          B = spline0(v);
          break;
          }
        case 1:
          {
          B = spline1(v);
          break;
          }
        case 2:
          {
          B = spline2(v);
          break;
          }
        case 3:
          {
          B = spline3(v);
          break;
          }
        default:
          {
          B = splineError(v);
          break;
          }
        }

      if( this->m_CloseDimension[dimension] )
        {
        idx[dimension] %=
          lattice->GetLargestPossibleRegion().GetSize()[dimension];
        }
      data += ( lattice->GetPixel( idx ) * B );
      }
    It.Set( data );
    }
}

template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::SetPhiLatticeParametricDomainParameters()
{
  typename PointDataImageType::PointType origin;
  typename PointDataImageType::SpacingType spacing;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType domain = this->m_Spacing[i] * static_cast<RealType>(
      this->m_Size[i] - 1 );

    unsigned int totalNumberOfSpans =
      this->m_PhiLattice->GetLargestPossibleRegion().GetSize()[i];
    if( !this->m_CloseDimension[i] )
      {
      totalNumberOfSpans -= this->m_SplineOrder[i];
      }

    spacing[i] = domain / static_cast<RealType>( totalNumberOfSpans );

    origin[i] = -0.5 * spacing[i] * ( this->m_SplineOrder[i] - 1 );
    }
  origin = this->m_Direction * origin;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] += this->m_Origin[i];
    }

  this->m_PhiLattice->SetOrigin( origin );
  this->m_PhiLattice->SetSpacing( spacing );
  this->m_PhiLattice->SetDirection( this->m_Direction );
}

template<class TInputPointSet, class TOutputImage>
typename BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::IndexType
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::NumberToIndex( const unsigned int number, const SizeType size )
{
  IndexType k;
  k[0] = 1;

  for ( unsigned int i = 1; i < ImageDimension; i++ )
    {
    k[i] = size[ImageDimension - i - 1] * k[i - 1];
    }

  unsigned int numberModulo = number;

  IndexType index;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    index[ImageDimension - i - 1] =
      static_cast<unsigned int>( numberModulo / k[ImageDimension - i - 1] );
    numberModulo %= k[ImageDimension - i - 1];
    }
  return index;
}

/**
 * Standard "PrintSelf" method
 */
template<class TInputPointSet, class TOutputImage>
void
BSplineScatteredDataPointSetToImageFilter<TInputPointSet, TOutputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
#if 0
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    this->m_Kernel[i]->Print( os, indent );
    }
#endif
  os << indent << "B-spline order: " << this->m_SplineOrder << std::endl;
  os << indent << "Number of control points: "
     << this->m_NumberOfControlPoints << std::endl;
  os << indent << "Close dimension: " << this->m_CloseDimension << std::endl;
  os << indent << "Number of levels " << this->m_NumberOfLevels << std::endl;
  os << indent << "Parametric domain" << std::endl;
  os << indent << "  Origin:    " << this->m_Origin << std::endl;
  os << indent << "  Spacing:   " << this->m_Spacing << std::endl;
  os << indent << "  Size:      " << this->m_Size << std::endl;
  os << indent << "  Direction: " << this->m_Direction << std::endl;
}
} // end namespace itk

#endif
