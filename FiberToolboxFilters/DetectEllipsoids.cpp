/* ============================================================================
* Copyright (c) 2009-2016 BlueQuartz Software, LLC
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright notice, this
* list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
* contributors may be used to endorse or promote products derived from this software
* without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, Data, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "DetectEllipsoids.h"

#define NOMINMAX

#include <QtCore/QDateTime>

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
#include <tbb/atomic.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>
#endif

#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/DataArrays/StringDataArray.hpp"

#include "FiberToolbox/Algorithms/ComputeGradient.h"
#include "FiberToolbox/FiberToolboxConstants.h"
#include "FiberToolbox/FiberToolboxVersion.h"

#include <cmath>
#include <limits>

// Include the MOC generated file for this class
#include "moc_DetectEllipsoids.cpp"

/**
 * @brief The DetectEllipsoidsImpl class implements a threaded algorithm that detects ellipsoids in a FeatureIds array
 */
class DetectEllipsoidsImpl
{
  DetectEllipsoids*               m_Filter;
  int*                            m_CellFeatureIdsPtr;
  Int8ArrayType::Pointer          m_EdgesArray;
  size_t                          m_CellFeatureIdsDims[3];
  UInt32ArrayType::Pointer        m_Corners;
  int32_t                         m_FeatureIdStart;
  int32_t                         m_FeatureIdEnd;
  size_t                          m_TotalNumOfFeatures;
  DE_ComplexDoubleVector          m_ConvCoords_X;
  DE_ComplexDoubleVector          m_ConvCoords_Y;
  DE_ComplexDoubleVector          m_ConvCoords_Z;
  QVector<size_t>                 m_ConvKernel_tDims;
  Int32ArrayType::Pointer         m_ConvOffsetArray;
  std::vector<double>             m_SmoothKernel;
  Int32ArrayType::Pointer         m_SmoothOffsetArray;

public:
	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
  DetectEllipsoidsImpl(DetectEllipsoids* filter, int* cellFeatureIdsPtr, Int8ArrayType::Pointer edgesArray, size_t cellFeatureIdsDims[3], UInt32ArrayType::Pointer corners, int32_t featureIdStart, int32_t featureIdEnd, size_t totalNumOfFeatures, DE_ComplexDoubleVector convCoords_X, DE_ComplexDoubleVector convCoords_Y, DE_ComplexDoubleVector convCoords_Z, QVector<size_t> kernel_tDims, Int32ArrayType::Pointer convOffsetArray, std::vector<double> smoothFil, Int32ArrayType::Pointer smoothOffsetArray) :
    m_Filter(filter),
    m_CellFeatureIdsPtr(cellFeatureIdsPtr),
    m_EdgesArray(edgesArray),
    m_Corners(corners),
    m_FeatureIdStart(featureIdStart),
    m_FeatureIdEnd(featureIdEnd),
    m_TotalNumOfFeatures(totalNumOfFeatures),
    m_ConvCoords_X(convCoords_X),
    m_ConvCoords_Y(convCoords_Y),
    m_ConvCoords_Z(convCoords_Z),
    m_ConvKernel_tDims(kernel_tDims),
    m_ConvOffsetArray(convOffsetArray),
    m_SmoothKernel(smoothFil),
    m_SmoothOffsetArray(smoothOffsetArray)
  {
    m_CellFeatureIdsDims[0] = cellFeatureIdsDims[0];
    m_CellFeatureIdsDims[1] = cellFeatureIdsDims[1];
    m_CellFeatureIdsDims[2] = cellFeatureIdsDims[2];
  }

	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
  virtual ~DetectEllipsoidsImpl()
  {
  }

	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
  void operator()() const
  {
    //std::cout << "Feature Start: " << m_FeatureIdStart << "\tFeature End: " << m_FeatureIdEnd << std::endl;
    for(int i = m_FeatureIdStart; i < m_FeatureIdEnd; i++)
    {
      size_t topL_X = m_Corners->getComponent(i, 0);
      size_t topL_Y = m_Corners->getComponent(i, 1);
      size_t topL_Z = m_Corners->getComponent(i, 2);
      size_t bottomR_X = m_Corners->getComponent(i, 3);
      size_t bottomR_Y = m_Corners->getComponent(i, 4);
      size_t bottomR_Z = m_Corners->getComponent(i, 5);

      size_t obj_xDim = bottomR_X - topL_X + 1;
      size_t obj_yDim = bottomR_Y - topL_Y + 1;
      size_t image_xDim = (bottomR_X+1) - (topL_X-1) + 1;
      size_t image_yDim = (bottomR_Y+1) - (topL_Y-1) + 1;
      size_t zDim = bottomR_Z - topL_Z + 1;

      QVector<size_t> image_tDims;
      image_tDims.push_back(image_xDim);
      image_tDims.push_back(image_yDim);
      image_tDims.push_back(zDim);

      QVector<size_t> obj_tDims;
      obj_tDims.push_back(obj_xDim);
      obj_tDims.push_back(obj_yDim);
      obj_tDims.push_back(zDim);

      QVector<size_t> cDims(1, 1);
      DoubleArrayType::Pointer featureObjArray = DoubleArrayType::CreateArray(image_tDims, cDims, "Feature Object");
      featureObjArray->initializeWithZeros();
      Int8ArrayType::Pointer edgeArray = Int8ArrayType::CreateArray(image_tDims, cDims, "Feature Object Edges");
      edgeArray->initializeWithZeros();

      for (size_t z = topL_Z; z <= bottomR_Z; z++)
      {
        for (size_t y = topL_Y; y <= bottomR_Y; y++)
        {
          for (size_t x = topL_X; x <= bottomR_X; x++)
          {
            size_t objX = x - topL_X;
            size_t objY = y - topL_Y;
            size_t objZ = z - topL_Z;
            size_t objIndex = (image_yDim * image_xDim * objZ) + (image_xDim * (objY+1)) + (objX+1);
            size_t imageIndex = (m_CellFeatureIdsDims[1] * m_CellFeatureIdsDims[0] * z) + (m_CellFeatureIdsDims[0] * y) + x;
            double featureValue = m_CellFeatureIdsPtr[imageIndex];
            int8_t edgesValue = m_EdgesArray->getValue(imageIndex);

            if (edgesValue > 1) { edgesValue = 1; }
            if (featureValue > 1.0) { featureValue = 1.0; }

            featureObjArray->setValue(objIndex, featureValue);
            edgeArray->setValue(objIndex, edgesValue);
          }
        }
      }

      // Convolute Gradient of object with filters
      ComputeGradient grad(featureObjArray, image_xDim, image_yDim);
      grad.compute();

      DoubleArrayType::Pointer gradX = grad.getGradX();
      DoubleArrayType::Pointer gradY = grad.getGradY();

      //std::cout << "Feature Id: " << i << "\tNumTuples = " << gradX->getNumberOfTuples() << std::endl;

      DE_ComplexDoubleVector gradX_conv = convoluteImage(gradX, m_ConvCoords_X, m_ConvOffsetArray, image_tDims);
      DE_ComplexDoubleVector gradY_conv = convoluteImage(gradY, m_ConvCoords_Y, m_ConvOffsetArray, image_tDims);

      DoubleArrayType::Pointer obj_conv_mag = DoubleArrayType::CreateArray(gradX_conv.size(), QVector<size_t>(1, 1), "obj_conv_mag");
      for (int i=0; i < gradX_conv.size(); i++)
      {
        std::complex<double> complexValue = gradX_conv[i] + gradY_conv[i];

        // Calculate magnitude of convolution
        double value = std::abs(complexValue);
        obj_conv_mag->setValue(i, value);
      }

      // Smooth Accumulator with Smoothing filter
      std::vector<double> obj_conv_mag_smooth = convoluteImage(obj_conv_mag, m_SmoothKernel, m_SmoothOffsetArray, image_tDims);
      double obj_conv_max = 0;
      for (int i=0; i<obj_conv_mag_smooth.size(); i++)
      {
        // Find max peak to set threshold
        if (obj_conv_mag_smooth[i] > obj_conv_max)
        {
          obj_conv_max = obj_conv_mag_smooth[i];
        }
        obj_conv_mag->setValue(i, obj_conv_mag_smooth[i]);
      }

      // Threshold convolution
      DoubleArrayType::Pointer obj_conv_thresh = DoubleArrayType::CreateArray(obj_conv_mag->getNumberOfTuples(), QVector<size_t>(1, 1), "obj_conv_thresh");
      obj_conv_thresh->initializeWithZeros();
      for (int i=0; i<obj_conv_thresh->getNumberOfTuples(); i++)
      {
        if (obj_conv_mag->getValue(i) > (0.7 * obj_conv_max))
        {
          double value = obj_conv_mag->getValue(i);
          obj_conv_thresh->setValue(i, value);
        }
      }

      if (m_Filter->getCancel())
      {
        return;
      }

      m_Filter->setFeaturesCompleted(m_Filter->getFeaturesCompleted() + 1);
      QString ss = QObject::tr("%1/%2").arg(m_Filter->getFeaturesCompleted()).arg(m_TotalNumOfFeatures);
      m_Filter->notifyStatusMessage(m_Filter->getMessagePrefix(), m_Filter->getHumanLabel(), ss);
    }
  }

	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
  template <typename T>
  std::vector<T> convoluteImage(DoubleArrayType::Pointer image, std::vector<T> kernel, Int32ArrayType::Pointer offsetArray, QVector<size_t> image_tDims) const
  {
    std::vector<T> convArray;

    std::vector<T> reverse_kernel = kernel;
    std::reverse(std::begin(reverse_kernel), std::end(reverse_kernel));

    int* offsetArrayPtr = offsetArray->getPointer(0);
    double* imageArray = image->getPointer(0);
    int offsetArrayNumOfComps = offsetArray->getNumberOfComponents();

    T accumulator;
    size_t xDim = image_tDims[0], yDim = image_tDims[1], zDim = image_tDims[2];
    int gradNumTuples = image->getNumberOfTuples();
    int reverseKernelCount = reverse_kernel.size();
    for (int i=0; i<gradNumTuples; i++)
    {
      if (m_Filter->getCancel())
      {
        return std::vector<T>();
      }

      int imageCurrentX = (i % xDim);
      int imageCurrentY = ((i / xDim) % yDim);
      int imageCurrentZ = (((i / xDim) / yDim) % zDim);
      for (int j = 0; j < reverseKernelCount; j++)
      {
        int currCoord_X = imageCurrentX + offsetArrayPtr[j*offsetArrayNumOfComps];
        int currCoord_Y = imageCurrentY + offsetArrayPtr[(j*offsetArrayNumOfComps) + 1];
        int currCoord_Z = imageCurrentZ + offsetArrayPtr[(j*offsetArrayNumOfComps) + 2];

        if (currCoord_X >= 0)
        {
          if (currCoord_X < xDim)
          {
            if (currCoord_Y >= 0)
            {
              if (currCoord_Y < yDim)
              {
                if (currCoord_Z >= 0)
                {
                  if (currCoord_Z < zDim)
                  {
                    int gradIndex = (yDim * xDim * currCoord_Z) + (xDim * currCoord_Y) + currCoord_X;

                    T kernelVal = reverse_kernel[j];
                    double imageVal = imageArray[gradIndex];
                    T value = kernelVal * imageVal;
                    accumulator += value;
                  }
                }
              }
            }
          }
        }
      }

      convArray.push_back(accumulator);
      accumulator = 0;
    }

    return convArray;
  }

	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
	// helper method - grabs index from matrix size
	int sub2ind(int matrixSize[2], int row, int col)
	{
		return row * matrixSize[1] + col;
	}

	// -----------------------------------------------------------------------------
	//
	// -----------------------------------------------------------------------------
	int* plotlineEdgeInter(int x0, int y0, int x1, int y1, DoubleArrayType::Pointer binImage, int imageDims[2])
	{
		int* edge = new int[2];

		// store initial point
		int xi = x0;
		int yi = y0;

		int dx = abs(x1 - x0);
		int dy = abs(y1 - y0);

		int sx, sy;

		if (x0 < x1)
			sx = 1;
		else
			sx = -1;

		if (y0 < y1)
			sy = 1;
		else
			sy = -1;

		int err = dx - dy;
		int e2;

		if (err == 0)
		{
			edge[0] = 0;
			edge[1] = 0;
			return edge;
		}

		//Search forward
		while (true)
		{
			if (binImage->getComponent(sub2ind(imageDims, x0, y0), 0) == 0)
			{
				if (x0 == x1 && y0 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 == x1 && y0 + 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 + 1 == x1 && y0 + 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 + 1 == x1 && y0 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 + 1 == x1 && y0 - 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 == x1 && y0 - 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 - 1 == x1 && y0 - 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 - 1 == x1 && y0 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else if (x0 - 1 == x1 && y0 + 1 == y1)
				{
					edge[0] = sub2ind(imageDims, x1, y1);
					break;
				}
				else
				{
					edge[0] = 0;
					edge[1] = 0;
					return edge;
				}
			}



			e2 = 2 * err;

			if (e2 > -dy)
			{
				err -= dy;
				x0 += sx;
			}

			if (e2 < dx)
			{
				err += dx;
				y0 += sy;
			}
		}

		// Reverse direction!
		x0 = xi;
		y0 = yi;

		if (x0 > x1)
			sx = 1;
		else
			sx = -1;

		if (y0 > y1)
			sy = 1;
		else
			sy = -1;

		err = dx - dy;

		while (true)
		{
			if (binImage->getComponent(sub2ind(imageDims, x0, y0), 0) == 0)
			{
				edge[1] = sub2ind(imageDims, x0, y0);
				break;
			}

			e2 = 2 * err;

			if (e2 > -dy)
			{
				err -= dy;
				x0 += sx;
			}

			if (e2 < dx)
			{
				err += dx;
				y0 += sy;
			}
		}

		return edge;
	}
};

double DetectEllipsoids::img_scale_length = 588.0;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectEllipsoids::DetectEllipsoids() :
  AbstractFilter(),
  m_MinFiberAxisLength(4),
  m_MaxFiberAxisLength(18),
  m_HoughTransformThreshold(0.5f),
  m_MinAspectRatio(0.4f),
  m_ImageScaleBarLength(100),
  m_ImageScaleBarUnits(ScaleBarUnits::MicronUnits),
  m_FeaturesCompleted(0)
{
  initialize();
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectEllipsoids::~DetectEllipsoids()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::initialize()
{
  setErrorCondition(0);
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::setupFilterParameters()
{
  FilterParameterVector parameters;

  parameters.push_back(SIMPL_NEW_INTEGER_FP("Min Fiber Axis Length (in units of image scale bar)", MinFiberAxisLength, FilterParameter::Parameter, DetectEllipsoids));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Max Fiber Axis Length (in units of image scale bar)", MaxFiberAxisLength, FilterParameter::Parameter, DetectEllipsoids));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Threshold for Hough Transform", HoughTransformThreshold, FilterParameter::Parameter, DetectEllipsoids));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Minimum Aspect Ratio", MinAspectRatio, FilterParameter::Parameter, DetectEllipsoids));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Length of Image Scale Bar (in units of image scale bar)", ImageScaleBarLength, FilterParameter::Parameter, DetectEllipsoids));

  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Units of Image Scale Bar");
    parameter->setPropertyName("ImageScaleBarUnits");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(DetectEllipsoids, this, ImageScaleBarUnits));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(DetectEllipsoids, this, ImageScaleBarUnits));

    QVector<QString> choices;
    choices.push_back("mm");
    choices.push_back("microns");
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, DetectEllipsoids, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int8, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Edge Array", EdgesArrayPath, FilterParameter::RequiredArray, DetectEllipsoids, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::CellFeature, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Active", ActiveArrayPath, FilterParameter::RequiredArray, DetectEllipsoids, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::dataCheck()
{
  setErrorCondition(0);

  getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType,AbstractFilter>(this, m_FeatureIdsArrayPath, QVector<size_t>(1, 1));
  getDataContainerArray()->getPrereqArrayFromPath<Int8ArrayType,AbstractFilter>(this, m_EdgesArrayPath, QVector<size_t>(1, 1));
  getDataContainerArray()->getPrereqArrayFromPath<BoolArrayType,AbstractFilter>(this, m_ActiveArrayPath, QVector<size_t>(1, 1));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true); // Set the fact that we are preflighting.
  emit preflightAboutToExecute(); // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck(); // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted(); // We are done preflighting this filter
  setInPreflight(false); // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::execute()
{
  initialize();
  dataCheck();
  if(getErrorCondition() < 0) { return; }

  /* Finding the top-left and bottom-right corners of each featureId  */
  Int32ArrayType::Pointer cellFeatureIds = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType,AbstractFilter>(this, m_FeatureIdsArrayPath, QVector<size_t>(1, 1));
  if (getErrorCondition() < 0) { return; }

  int* cellFeatureIdsPtr = cellFeatureIds->getPointer(0);
  if (cellFeatureIdsPtr != nullptr)
  {
    int featureId = 0;
    size_t numComps = 6;
    QVector<size_t> cDims(1, numComps);
    int err = 0;
    AttributeMatrix::Pointer activeAM = getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, m_ActiveArrayPath, err);
    UInt32ArrayType::Pointer corners = UInt32ArrayType::CreateArray(activeAM->getTupleDimensions(), cDims, "Corners of Feature");
    for (int i=0; i<corners->getNumberOfTuples(); i++)
    {
      corners->setComponent(i, 0, std::numeric_limits<uint32_t>::max());
      corners->setComponent(i, 1, std::numeric_limits<uint32_t>::max());
      corners->setComponent(i, 2, std::numeric_limits<uint32_t>::max());
      corners->setComponent(i, 3, std::numeric_limits<uint32_t>::min());
      corners->setComponent(i, 4, std::numeric_limits<uint32_t>::min());
      corners->setComponent(i, 5, std::numeric_limits<uint32_t>::min());
    }

    DataContainer::Pointer featureIdsDC = getDataContainerArray()->getDataContainer(m_FeatureIdsArrayPath.getDataContainerName());

    size_t xDim = 0, yDim = 0, zDim = 0;
    size_t dims[3];
    featureIdsDC->getGeometryAs<ImageGeom>()->getDimensions(xDim, yDim, zDim);
    featureIdsDC->getGeometryAs<ImageGeom>()->getDimensions(dims);

    size_t index = 0;
    for(size_t z = 0; z < zDim; z++)
    {
      for(size_t y = 0; y < yDim; y++)
      {
        for(size_t x = 0; x < xDim; x++)
        {
          index = (yDim * xDim * z) + (xDim * y) + x ; // Index into cellFeatureIds array

          featureId = cellFeatureIdsPtr[index];

          uint32_t* featureCorner = corners->getPointer(featureId * numComps);

          uint32_t val = featureCorner[0];
          if(x < featureCorner[0])
          {
            featureCorner[0] = x;
          }
          val = featureCorner[1];
          if(y < featureCorner[1])
          {
            featureCorner[1] = y;
          }
          val = featureCorner[2];
          if(z < featureCorner[2])
          {
            featureCorner[2] = z;
          }

          val = featureCorner[3];
          if(x > featureCorner[3])
          {
            featureCorner[3] = x;
          }
          val = featureCorner[4];
          if(y > featureCorner[4])
          {
            featureCorner[4] = y;
          }
          val = featureCorner[5];
          if(z > featureCorner[5])
          {
            featureCorner[5] = z;
          }
        }
      }
    }

    double img_pix_length = m_ImageScaleBarLength / img_scale_length;
    double axis_min = std::round( m_MinFiberAxisLength / img_pix_length );
    double axis_max = std::round( m_MaxFiberAxisLength / img_pix_length );

    QVector<size_t> orient_tDims;
    QVector<size_t> hc_tDims;
    DoubleArrayType::Pointer orientArray = orientationFilter(axis_min, axis_max, orient_tDims);
    DE_ComplexDoubleVector houghCircleVector = houghCircleFilter(axis_min, axis_max, hc_tDims);

    if (orientArray->getNumberOfTuples() != houghCircleVector.size())
    {
      setErrorCondition(-31000);
      QString ss = QObject::tr("There was an internal error.  Please ask the DREAM.3D developers for more information.");
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    // This function fills the xCoords, yCoords, and zCoords arrays with values
    DE_ComplexDoubleVector convCoords_X;
    DE_ComplexDoubleVector convCoords_Y;
    DE_ComplexDoubleVector convCoords_Z;
    convolutionFilter(orientArray, houghCircleVector, convCoords_X, convCoords_Y, convCoords_Z);

    Int8ArrayType::Pointer edgesArray = getDataContainerArray()->getPrereqArrayFromPath<Int8ArrayType,AbstractFilter>(this, m_EdgesArrayPath, QVector<size_t>(1, 1));

    Int32ArrayType::Pointer convOffsetArray = createOffsetArray(orient_tDims);

    int n_size = 3;
    std::vector<double> smoothFil = smoothingFilter(n_size);
    QVector<size_t> smooth_tDims;
    smooth_tDims.push_back(2*n_size+1);
    smooth_tDims.push_back(2*n_size+1);
    smooth_tDims.push_back(1);
    Int32ArrayType::Pointer smoothOffsetArray = createOffsetArray(smooth_tDims);

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    tbb::task_scheduler_init init;
    bool doParallel = false;
#endif

    size_t totalNumOfFeatures = activeAM->getNumberOfTuples();

    QString ss = QObject::tr("%1/%2").arg(getFeaturesCompleted()).arg(totalNumOfFeatures);
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

    qint64 millis = QDateTime::currentMSecsSinceEpoch();

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::task_group* g = new tbb::task_group;
      int threads = init.default_num_threads();
      unsigned int numOfTasks = (totalNumOfFeatures-1) / threads;

      int32_t start = 1;
      int32_t end = 0 + numOfTasks;
      for (int i=0; i<threads; i++)
      {
        g->run(DetectEllipsoidsImpl(this, cellFeatureIdsPtr, edgesArray, dims, corners, start, end, totalNumOfFeatures, convCoords_X, convCoords_Y, convCoords_Z, orient_tDims, convOffsetArray, smoothFil, smoothOffsetArray));
        start = end;
        end = end + numOfTasks;
        if(end >= totalNumOfFeatures)
        {
          end = totalNumOfFeatures;
        }
      }

      g->wait();
      delete g;
    }
    else
#endif
    {
      DetectEllipsoidsImpl impl(this, cellFeatureIdsPtr, edgesArray, dims, corners, 21, 22, totalNumOfFeatures, convCoords_X, convCoords_Y, convCoords_Z, orient_tDims, convOffsetArray, smoothFil, smoothOffsetArray);
      impl();
    }

    qint64 secs = (QDateTime::currentMSecsSinceEpoch() - millis) / 1000;
    std::cout << "Seconds: " << secs;
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DoubleArrayType::Pointer DetectEllipsoids::orientationFilter(int minAxisLength, int maxAxisLength, QVector<size_t> &tDims)
{
  double doubleMax = static_cast<double>(maxAxisLength);
  double doubleMin = static_cast<double>(minAxisLength);
  double doubleMin_squared = doubleMin*doubleMin;
  double doubleMax_squared = doubleMax*doubleMax;

  size_t xDim = 2*maxAxisLength+1;
  size_t yDim = 2*maxAxisLength+1;
  size_t zDim = 1;  // This can be changed later to handle 3-dimensions
  QVector<size_t> cDims(1, 3);
  tDims.clear();
  tDims.push_back(xDim);
  tDims.push_back(yDim);
  tDims.push_back(zDim);
  DoubleArrayType::Pointer orientationCoords = DoubleArrayType::CreateArray(tDims, cDims, "Orientation Coordinates");

  for (int z = 1; z <= zDim; z++)
  {
    for (int y = 1; y <= yDim; y++)
    {
      for (int x = 1; x <= xDim; x++)
      {
        int xIdx = x - 1;
        int yIdx = y - 1;
        int zIdx = z - 1;
        size_t index = (yDim * xDim * zIdx) + (xDim * yIdx) + xIdx;

        double m = static_cast<double>(y)-1.0-doubleMax;
        double n = static_cast<double>(x)-1.0-doubleMax;
        double theta = std::atan2(n,m);

        if( (m*m) + (n*n) >= doubleMin_squared && (m*m) + (n*n) <= doubleMax_squared)
        {
          orientationCoords->setComponent(index, 0, std::cos(theta));
          orientationCoords->setComponent(index, 1, std::sin(theta));
          orientationCoords->setComponent(index, 2, 0); // This can be changed later to handle 3-dimensions
        }
        else
        {
          orientationCoords->setComponent(index, 0, 0);
          orientationCoords->setComponent(index, 1, 0);
          orientationCoords->setComponent(index, 2, 0);
        }
      }
    }
  }

  return orientationCoords;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DE_ComplexDoubleVector DetectEllipsoids::houghCircleFilter(int minAxisLength, int maxAxisLength, QVector<size_t> &tDims)
{
  size_t xDim = 2*maxAxisLength+1;
  size_t yDim = 2*maxAxisLength+1;
  size_t zDim = 1;  // This can be changed later to handle 3-dimensions
  size_t totalElements = xDim * yDim * zDim;
  tDims.clear();
  tDims.push_back(xDim);
  tDims.push_back(yDim);
  tDims.push_back(zDim);
  DE_ComplexDoubleVector houghCircleCoords(totalElements);
  int minAxisLength_squared = minAxisLength*minAxisLength;
  int maxAxisLength_squared = maxAxisLength*maxAxisLength;

  for (int z = 1; z <= zDim; z++)
  {
    for (int y = 1; y <= yDim; y++)
    {
      for (int x = 1; x <= xDim; x++)
      {
        int xIdx = x - 1;
        int yIdx = y - 1;
        int zIdx = z - 1;
        size_t index = (yDim * xDim * zIdx) + (xDim * yIdx) + xIdx;

        double m = y-1-maxAxisLength;
        double n = x-1-maxAxisLength;
        double phi = ( std::sqrt( (m*m) + (n*n) ) - minAxisLength ) / ( maxAxisLength - minAxisLength );

        if( (m*m) + (n*n) >= minAxisLength_squared && (m*m) + (n*n) <= maxAxisLength_squared)
        {
          std::complex<double> complexVal(std::cos(2*M_PI*phi), std::sin(2*M_PI*phi));
          std::complex<double> value = 1.0/2.0/M_PI/std::sqrt((m*m) + (n*n)) * complexVal;
          houghCircleCoords[index] = value;
        }
        else
        {
          houghCircleCoords[index] = 0;
        }
      }
    }
  }

  return houghCircleCoords;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::convolutionFilter(DoubleArrayType::Pointer orientationFilter, DE_ComplexDoubleVector houghCircleFilter, DE_ComplexDoubleVector &convCoords_X, DE_ComplexDoubleVector &convCoords_Y, DE_ComplexDoubleVector &convCoords_Z)
{
  if (orientationFilter->getNumberOfTuples() != houghCircleFilter.size()
      || orientationFilter->getNumberOfComponents() != 3)
  { return; }

  for (int i=0; i<orientationFilter->getNumberOfTuples(); i++)
  {
    std::complex<double> hcValue = houghCircleFilter[i];

    double orientValue_X = orientationFilter->getComponent(i, 0);
    std::complex<double> valueX = orientValue_X * hcValue;
    convCoords_X.push_back(valueX);

    double orientValue_Y = orientationFilter->getComponent(i, 1);
    std::complex<double> valueY = orientValue_Y * hcValue;
    convCoords_Y.push_back(valueY);

    double orientValue_Z = orientationFilter->getComponent(i, 2);
    std::complex<double> valueZ = orientValue_Z * hcValue;
    convCoords_Z.push_back(valueZ);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::vector<double> DetectEllipsoids::smoothingFilter(int n_size)
{
  int xDim = 2*n_size+1;
  int yDim = 2*n_size+1;
  int zDim = 1;
  std::vector<double> smooth(xDim*yDim*zDim);
  int n_size_squared = n_size*n_size;

  for (int z = 0; z < zDim; z++)
  {
    for (int y = 0; y < yDim; y++)
    {
      for (int x = 0; x < xDim; x++)
      {
        int m = y-n_size;
        int n = x-n_size;
        int index = (yDim * xDim * z) + (xDim * y) + x;

        if( ((m*m) + (n*n)) <= n_size_squared)
        {
          smooth[index] = 1;
        }
        else
        {
          smooth[index] = 0;
        }
      }
    }
  }

  return smooth;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Int32ArrayType::Pointer DetectEllipsoids::createOffsetArray(QVector<size_t> kernel_tDims)
{
  QVector<size_t> cDims(1, 3);
  Int32ArrayType::Pointer offsetArray = Int32ArrayType::CreateArray(kernel_tDims, cDims, "Coordinate Array");
  size_t xDim = kernel_tDims[0], yDim = kernel_tDims[1], zDim = kernel_tDims[2];
  int index = 0;

  for (int z = 0; z < zDim; z++)
  {
    for (int y = 0; y < yDim; y++)
    {
      for (int x = 0; x < xDim; x++)
      {
        int xVal = x - xDim/2;
        offsetArray->setComponent(index, 0, xVal);
        int yVal = y - yDim/2;
        offsetArray->setComponent(index, 1, yVal);
        int zVal = z - zDim/2;
        offsetArray->setComponent(index, 2, zVal);
        index++;
      }
    }
  }

  return offsetArray;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer DetectEllipsoids::newFilterInstance(bool copyFilterParameters)
{
  DetectEllipsoids::Pointer filter = DetectEllipsoids::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getCompiledLibraryName()
{ return FiberToolboxConstants::FiberToolboxBaseName; }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getBrandingString()
{
  return "FiberToolbox";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getFilterVersion()
{
  QString version;
  QTextStream vStream(&version);
  vStream <<  FiberToolbox::Version::Major() << "." << FiberToolbox::Version::Minor() << "." << FiberToolbox::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getGroupName()
{ return SIMPL::FilterGroups::Unsupported; }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getSubGroupName()
{ return "FiberToolbox"; }

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString DetectEllipsoids::getHumanLabel()
{ return "Detect Ellipsoids"; }

