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

#define STORE_PIXEL_VALUES(array, count)\
  array->setComponent(count, 0, xc+x);\
  array->setComponent(count, 1, yc+y);\
  count++;\
  array->setComponent(count, 0, xc-x);\
  array->setComponent(count, 1, yc-y);\
  count++;\

/**
 * @brief The DetectEllipsoidsImpl class implements a threaded algorithm that detects ellipsoids in a FeatureIds array
 */
class DetectEllipsoidsImpl
{
  DetectEllipsoids*               m_Filter;
  int*                            m_CellFeatureIdsPtr;
  QVector<size_t>                 m_CellFeatureIdsDims;
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
  double                          m_Axis_Min;
  double                          m_Axis_Max;
  float                           m_TolEllipse;
  float                           m_Ba_Min;
  DoubleArrayType::Pointer        m_Cenx;
  DoubleArrayType::Pointer        m_Ceny;
  DoubleArrayType::Pointer        m_Majaxis;
  DoubleArrayType::Pointer        m_Minaxis;
  DoubleArrayType::Pointer        m_Rotangle;
  DoubleArrayType::Pointer        m_AdditionalEllipses;

public:
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  DetectEllipsoidsImpl(DetectEllipsoids* filter, int* cellFeatureIdsPtr, QVector<size_t> cellFeatureIdsDims, UInt32ArrayType::Pointer corners, int32_t featureIdStart, int32_t featureIdEnd, size_t m_TotalNumberOfFeatures, DE_ComplexDoubleVector convCoords_X, DE_ComplexDoubleVector convCoords_Y, DE_ComplexDoubleVector convCoords_Z, QVector<size_t> kernel_tDims, Int32ArrayType::Pointer convOffsetArray, std::vector<double> smoothFil, Int32ArrayType::Pointer smoothOffsetArray, double axis_min, double axis_max,
  float tol_ellipse, float ba_min, DoubleArrayType::Pointer cenx, DoubleArrayType::Pointer ceny, DoubleArrayType::Pointer majaxis, DoubleArrayType::Pointer minaxis, DoubleArrayType::Pointer rotangle, DoubleArrayType::Pointer additionalEllipses) :
    m_Filter(filter),
    m_CellFeatureIdsPtr(cellFeatureIdsPtr),
    m_CellFeatureIdsDims(cellFeatureIdsDims),
    m_Corners(corners),
    m_FeatureIdStart(featureIdStart),
    m_FeatureIdEnd(featureIdEnd),
    m_TotalNumOfFeatures(m_TotalNumberOfFeatures),
    m_ConvCoords_X(convCoords_X),
    m_ConvCoords_Y(convCoords_Y),
    m_ConvCoords_Z(convCoords_Z),
    m_ConvKernel_tDims(kernel_tDims),
    m_ConvOffsetArray(convOffsetArray),
    m_SmoothKernel(smoothFil),
    m_SmoothOffsetArray(smoothOffsetArray),
    m_Axis_Min(axis_min),
    m_Axis_Max(axis_max),
    m_TolEllipse(tol_ellipse),
    m_Ba_Min(ba_min),
    m_Cenx(cenx),
    m_Ceny(ceny),
    m_Majaxis(majaxis),
    m_Minaxis(minaxis),
    m_Rotangle(rotangle),
    m_AdditionalEllipses(additionalEllipses)
  {

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
    // Initialize temporary arrays for candidate ellipse and accumulation counter
    DoubleArrayType::Pointer cenx_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "cenx_can"); // x-coordinate of ellipse
    for (int i=0; i < cenx_can->getNumberOfTuples(); i++)
    {
      cenx_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer ceny_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "ceny_can"); // y-coordinate of ellipse
    for (int i=0; i < ceny_can->getNumberOfTuples(); i++)
    {
      ceny_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer maj_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "maj_can"); // major semi-axis
    for (int i=0; i < maj_can->getNumberOfTuples(); i++)
    {
      maj_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer min_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "min_can"); // minor semi-axis
    for (int i=0; i < min_can->getNumberOfTuples(); i++)
    {
      min_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer rot_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "rot_can"); // Counter clockwise rotation from x-axis
    for (int i=0; i < rot_can->getNumberOfTuples(); i++)
    {
      rot_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer accum_can = DoubleArrayType::CreateArray(300, QVector<size_t>(1, 1), "accum_can"); // Accumulation matrix
    for (int i=0; i < accum_can->getNumberOfTuples(); i++)
    {
      accum_can->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    for(size_t featureId = m_FeatureIdStart; featureId < m_FeatureIdEnd; featureId++)
    {
      size_t topL_X = m_Corners->getComponent(featureId, 0);
      size_t topL_Y = m_Corners->getComponent(featureId, 1);
      size_t topL_Z = m_Corners->getComponent(featureId, 2);
      size_t bottomR_X = m_Corners->getComponent(featureId, 3);
      size_t bottomR_Y = m_Corners->getComponent(featureId, 4);
      size_t bottomR_Z = m_Corners->getComponent(featureId, 5);

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

            if (featureValue == featureId)
            {
              featureValue = 1.0;
            }
            else
            {
              featureValue = 0.0;
            }

            featureObjArray->setValue(objIndex, featureValue);
          }
        }
      }

      double min_pix = std::round(M_PI * m_Axis_Min * m_Axis_Min / 2);
      SizeTArrayType::Pointer objPixelsArray = findNonZeroIndices<double>(featureObjArray, image_tDims);

      while (objPixelsArray->getNumberOfTuples() > min_pix)
      {
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

        QList<int> obj_ext_indices = findExtrema(obj_conv_thresh, image_tDims);
        int obj_ext_num = obj_ext_indices.size();

        if(obj_ext_num > 3)
        {
          obj_ext_num = 3;
        }

        int detectedObjIdx = 0;
        for (; detectedObjIdx < obj_ext_num; detectedObjIdx++)
        {
          size_t obj_ext_index = obj_ext_indices[detectedObjIdx];

          // Get indices of extrema value
          size_t obj_ext_y = obj_ext_index % image_xDim;
          size_t obj_ext_x = ((obj_ext_index / image_xDim) % image_yDim);

          size_t mask_rad = m_Axis_Max + m_Axis_Min;

          int mask_min_x = obj_ext_x - mask_rad;
          if ( mask_min_x < 1)
          {
            mask_min_x = 1;
          }

          int mask_min_y = obj_ext_y - mask_rad;
          if ( mask_min_y < 1)
          {
            mask_min_y = 1;
          }

          int mask_max_x = obj_ext_x + mask_rad;
          if ( mask_max_x > image_tDims[0])
          {
            mask_max_x = image_tDims[0];
          }

          int mask_max_y = obj_ext_y + mask_rad;
          if ( mask_max_y > image_tDims[1])
          {
            mask_max_y = image_tDims[1];
          }

          QVector<size_t> edge_tDims;
          edge_tDims.push_back(mask_max_x - mask_min_x + 1);
          edge_tDims.push_back(mask_max_y - mask_min_y + 1);

          DoubleArrayType::Pointer obj_mask = DoubleArrayType::CreateArray(edge_tDims, QVector<size_t>(1, 1), "obj_mask");
          obj_mask->initializeWithZeros();

          for (size_t y = mask_min_y; y < mask_max_y; y++)
          {
            for (size_t x = mask_min_x; x < mask_max_x; x++)
            {
              int featureObj_index = image_tDims[0] * y + x;
              int edge_index = edge_tDims[0] * (y - mask_min_y + 1) + (x - mask_min_x + 1);
              obj_mask->setValue(edge_index, featureObjArray->getValue(featureObj_index));
            }
          }

          // Compute the edge array
          Int8ArrayType::Pointer edgeArray = findEdges<double>(obj_mask, edge_tDims);

          SizeTArrayType::Pointer obj_edges = findNonZeroIndices<int8_t>(edgeArray, edge_tDims);

          SizeTArrayType::Pointer obj_edge_pair_a = SizeTArrayType::CreateArray(obj_edges->getNumberOfTuples(), QVector<size_t>(1, 1), "obj_edge_pair_a");
          SizeTArrayType::Pointer obj_edge_pair_b = SizeTArrayType::CreateArray(obj_edges->getNumberOfTuples(), QVector<size_t>(1, 1), "obj_edge_pair_b");
          SizeTArrayType::Pointer obj_edge_pair_a1 = SizeTArrayType::CreateArray(obj_edges->getNumberOfTuples(), QVector<size_t>(1, 2), "obj_edge_pair_a1");
          SizeTArrayType::Pointer obj_edge_pair_b1 = SizeTArrayType::CreateArray(obj_edges->getNumberOfTuples(), QVector<size_t>(1, 2), "obj_edge_pair_b1");
          int count = 0;
          for (int j=0; j<obj_edges->getNumberOfTuples(); j++)
          {
            QPair<size_t,size_t> edgeIndices = plotlineEdgeInter(obj_ext_x, obj_ext_y, obj_edges->getComponent(j, 1), obj_edges->getComponent(j, 0), featureObjArray, image_tDims);

            if (edgeIndices.first != 0 && edgeIndices.second != 0)
            {
              obj_edge_pair_a->setValue(count, edgeIndices.first);
              obj_edge_pair_b->setValue(count, edgeIndices.second);

              size_t x1 = (edgeIndices.first % image_tDims[0]);
              size_t y1 = ((edgeIndices.first / image_tDims[0]) % image_tDims[1]);
              size_t x2 = (edgeIndices.second % image_tDims[0]);
              size_t y2 = ((edgeIndices.second / image_tDims[0]) % image_tDims[1]);

              obj_edge_pair_a1->setComponent(count, 0, x1);
              obj_edge_pair_a1->setComponent(count, 1, y1);
              obj_edge_pair_b1->setComponent(count, 0, x2);
              obj_edge_pair_b1->setComponent(count, 1, y2);
              count++;
            }
          }
          obj_edge_pair_a->resize(count);
          obj_edge_pair_b->resize(count);
          obj_edge_pair_a1->resize(count);
          obj_edge_pair_b1->resize(count);

          // Search current object for ellipses
          size_t can_num = 0;
          for (int k = 0; k < obj_edge_pair_a1->getNumberOfTuples(); k++)
          {
            analyzeEdgePair(obj_edge_pair_a1, obj_edge_pair_b1, k, edge_tDims, edgeArray, can_num, cenx_can, ceny_can, maj_can, min_can, rot_can, accum_can);
          }

          if(can_num > 0) // Assume best match is the ellipse
          {
            // Increment the ellipse counter
            m_Filter->setEllipse_Count(m_Filter->getEllipse_Count() + 1);

            int accum_idx = getIdOfMax<double>(accum_can);

            /* If this is another ellipse in the same overall object,
             * create a new feature id and resize our output arrays */
            size_t objId = featureId;
            if (detectedObjIdx > 0)
            {
              objId = m_Filter->getUniqueFeatureId();
              m_Cenx->resize(objId + 1);
              m_Ceny->resize(objId + 1);
              m_Majaxis->resize(objId + 1);
              m_Minaxis->resize(objId + 1);
              m_Rotangle->resize(objId + 1);

              int extraEllipsesIndex = m_Filter->getAdditionalEllipsesIndex();
              m_AdditionalEllipses->setComponent(extraEllipsesIndex, 0, featureId);
              m_AdditionalEllipses->setComponent(extraEllipsesIndex, 1, objId);
  //            m_AdditionalEllipses->setComponent(extraEllipsesIndex, 2, featureId);
  //            m_AdditionalEllipses->setComponent(extraEllipsesIndex, 3, objId);
  //            m_AdditionalEllipses->setComponent(extraEllipsesIndex, 4, featureId);
  //            m_AdditionalEllipses->setComponent(extraEllipsesIndex, 5, objId);
            }

            // Store ellipse parameters
            double cenx_val = cenx_can->getValue(accum_idx);
            double ceny_val = ceny_can->getValue(accum_idx);
            double majaxis_val = maj_can->getValue(accum_idx);
            double minaxis_val = min_can->getValue(accum_idx);

            double rotangle_val = rot_can->getValue(accum_idx);
            while (rotangle_val > M_PI_2 || rotangle_val < -M_PI_2)
            {
              if (rotangle_val > M_PI_2)
              {
                rotangle_val = rotangle_val - M_PI;
              }
              else
              {
                rotangle_val = rotangle_val + M_PI;
              }
            }

            m_Cenx->setValue(objId, cenx_val);
            m_Ceny->setValue(objId, ceny_val);
            m_Majaxis->setValue(objId, majaxis_val);
            m_Minaxis->setValue(objId, minaxis_val);
            m_Rotangle->setValue(objId, rotangle_val);

            /* Remove pixels from object (including additional 1 pixel thick boarder)
             * and reassign remaining pixels to array */
            DoubleArrayType::Pointer featureObjOnesArray = DoubleArrayType::CreateArray(image_tDims, QVector<size_t>(1, 1), "featureObjOnesArray");
            featureObjOnesArray->initializeWithValue(1);

            DoubleArrayType::Pointer I_tmp = fillEllipse(featureObjOnesArray, image_tDims, cenx_val, ceny_val, majaxis_val+1, minaxis_val+1, rotangle_val, 0);

            for (int i = 0; i < I_tmp->getNumberOfTuples(); i++)
            {
              double value = featureObjArray->getValue(i) * I_tmp->getValue(i);
              featureObjArray->setValue(i, value);
            }

            // Translate center of ellipse into real space
            size_t obj_x_min = topL_Y;
            size_t obj_y_min = topL_X;
            m_Cenx->setValue(objId, m_Cenx->getValue(objId) + obj_x_min);
            m_Ceny->setValue(objId, m_Ceny->getValue(objId) + obj_y_min);

            // Clean Accumulator
            accum_can->initializeWithZeros();

            objPixelsArray = findNonZeroIndices<double>(featureObjArray, image_tDims);

            break;
          }
        }

        if(detectedObjIdx == obj_ext_num)
        {
          break;
        }
      }

      if (m_Filter->getCancel())
      {
        return;
      }

      m_Filter->notifyFeatureCompleted();
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  template <typename T>
  Int8ArrayType::Pointer findEdges(typename DataArray<T>::Pointer imagePtr, QVector<size_t> image_tDims) const
  {
    Int8ArrayType::Pointer edgeArray = Int8ArrayType::CreateArray(image_tDims, QVector<size_t>(1, 1), "Edge Array");
    edgeArray->initializeWithZeros();

    size_t xDim = image_tDims[0];
    size_t yDim = image_tDims[1];

    for (int y = 0; y < yDim; y++)
    {
      for (int x = 0; x < xDim; x++)
      {
        size_t index = (xDim * y) + x;
        T currentValue = imagePtr->getValue(index);
        if (currentValue != 0)
        {
          if ((y - 1) >= 0)
          {
            int topIdx = (xDim * (y - 1)) + x;
            if (imagePtr->getValue(topIdx) == 0)
            {
              edgeArray->setValue(index, 1);
              continue;
            }
          }
          if ((y + 1) < yDim)
          {
            int bottomIdx = (xDim * (y + 1)) + x;
            if (imagePtr->getValue(bottomIdx) == 0)
            {
              edgeArray->setValue(index, 1);
              continue;
            }
          }
          if ((x - 1) >= 0)
          {
            int leftIdx = (xDim * y) + (x - 1);
            if (imagePtr->getValue(leftIdx) == 0)
            {
              edgeArray->setValue(index, 1);
              continue;
            }
          }
          if ((x + 1) < xDim)
          {
            int rightIdx = (xDim * y) + (x + 1);
            if (imagePtr->getValue(rightIdx) == 0)
            {
              edgeArray->setValue(index, 1);
              continue;
            }
          }
        }
      }
    }

    return edgeArray;
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
  QList<int> findExtrema(DoubleArrayType::Pointer thresholdArray, QVector<size_t> tDims) const
  {
    QSet<int> extrema;
    size_t xDim = tDims[0];
    size_t yDim = tDims[1];

    QSet<int> extremaCols;
    QList<int> extremaCol_flatIndices;

    // Search peaks through columns
    for (int x = 0; x < xDim; x++)
    {
      double value = 0;
      int extrema_Y = -1;
      int flat_index = -1;
      for (int y = 0; y < yDim; y++)
      {
        int index = (xDim * y) + x;
        double threshValue = thresholdArray->getValue(index);
        if (threshValue > value)
        {
          extrema_Y = y;
          value = threshValue;
          flat_index = index;
        }
      }
      if (extrema_Y >= 0)
      {
        extremaCols.insert(extrema_Y);
        extremaCol_flatIndices.push_back(flat_index);
      }
    }

    // Search peaks through rows, on columns with extrema points:
    QList<int> extremaRow_flatIndices;
    QList<int> extremaCols_list = extremaCols.toList();
    int extremaColSize = extremaCols_list.size();
    for (int i=0; i<extremaColSize; i++)
    {
      int y = extremaCols_list[i];
      int extremaIndex = xDim * y + 0;  // Initialize to index of first value in extremaRow_flatIndices
      for (int x=1; x<xDim; x++)
      {
        int index = (xDim * y) + x;
        if (thresholdArray->getValue(index) > thresholdArray->getValue(extremaIndex))
        {
          extremaIndex = index;
        }
      }
      extremaRow_flatIndices.push_back(extremaIndex);
    }

    // Peaks in rows and in columns (intersection):
    QSet<int> rcIntersection;
    int extremaColFlatIndicesSize = extremaCol_flatIndices.size();
    for (int i=0; i<extremaColFlatIndicesSize; i++)
    {
      if (extremaRow_flatIndices.contains(extremaCol_flatIndices[i]))
      {
        rcIntersection.insert(extremaCol_flatIndices[i]);
        extremaRow_flatIndices.removeAll(extremaCol_flatIndices[i]);
      }
    }

    int extremaRowFlatIndicesSize = extremaRow_flatIndices.size();
    for (int i=0; i<extremaRowFlatIndicesSize; i++)
    {
      if (extremaCol_flatIndices.contains(extremaRow_flatIndices[i]))
      {
        rcIntersection.insert(extremaRow_flatIndices[i]);
        extremaRow_flatIndices.removeAll(extremaCol_flatIndices[i]);
      }
    }

    // Peaks through diagonals
    QList<int> rcIntersectionList = rcIntersection.toList();
    int rcIntersectionListSize = rcIntersectionList.size();
    for (int i = 0; i < rcIntersectionListSize; i++)
    {
      int extremaIndex = rcIntersectionList[i];
      int x = (extremaIndex % xDim);
      int y = ((extremaIndex / xDim) % yDim);

      int xDiag = x - 1;
      int yDiag = y - 1;
      while (xDiag >= 0 && yDiag >= 0)
      {
        int diagIndex = (xDim * yDiag) + xDiag;
        if (thresholdArray->getValue(diagIndex) > thresholdArray->getValue(extremaIndex))
        {
          extremaIndex = diagIndex;
        }
        xDiag--;
        yDiag--;
      }

      xDiag = x + 1;
      yDiag = y + 1;
      while (xDiag < xDim && yDiag < yDim)
      {
        int diagIndex = (xDim * yDiag) + xDiag;
        if (thresholdArray->getValue(diagIndex) > thresholdArray->getValue(extremaIndex))
        {
          extremaIndex = diagIndex;
        }
        xDiag++;
        yDiag++;
      }

      xDiag = x - 1;
      yDiag = y + 1;
      while (xDiag >= 0 && yDiag < yDim)
      {
        int diagIndex = (xDim * yDiag) + xDiag;
        if (thresholdArray->getValue(diagIndex) > thresholdArray->getValue(extremaIndex))
        {
          extremaIndex = diagIndex;
        }
        xDiag--;
        yDiag++;
      }

      xDiag = x + 1;
      yDiag = y - 1;
      while (xDiag < xDim && yDiag >= 0)
      {
        int diagIndex = (xDim * yDiag) + xDiag;
        if (thresholdArray->getValue(diagIndex) > thresholdArray->getValue(extremaIndex))
        {
          extremaIndex = diagIndex;
        }
        xDiag++;
        yDiag--;
      }

      extrema.insert(extremaIndex);
    }

    return extrema.toList();
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  // helper method - grabs index from matrix size
  size_t sub2ind(QVector<size_t> tDims, size_t row, size_t col) const
  {
    return row * tDims[0] + col;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  QPair<size_t,size_t> plotlineEdgeInter(int x0, int y0, int x1, int y1, DoubleArrayType::Pointer binImage, QVector<size_t> imageDims) const
  {
    QPair<size_t,size_t> edge;

    // Store Initial Point
    size_t xi = x0;
    size_t yi = y0;

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
      edge.first = 0;
      edge.second = 0;
      return edge;
    }

    //Search forward
    while (true)
    {
      if (binImage->getValue(sub2ind(imageDims, x0, y0)) == 0)
      {
        if (x0 == x1 && y0 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 == x1 && y0 + 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 + 1 == x1 && y0 + 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 + 1 == x1 && y0 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 + 1 == x1 && y0 - 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 == x1 && y0 - 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 - 1 == x1 && y0 - 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 - 1 == x1 && y0 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else if (x0 - 1 == x1 && y0 + 1 == y1)
        {
          edge.first = sub2ind(imageDims, x1, y1);
          break;
        }
        else
        {
          edge.first = 0;
          edge.second = 0;
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
    x0=xi;
    y0=yi;

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
        edge.second = sub2ind(imageDims, x0, y0);
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

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  template <typename T>
  int getIdOfMax(typename DataArray<T>::Pointer array) const
  {
    if (array.get() == nullptr)
      return -1;

    int arrayLength = array->getNumberOfTuples() * array->getNumberOfComponents();

    if (arrayLength <= 0)
      return -1;

    double maxId = 0;
    double maxValue = std::numeric_limits<T>::min();

    for (int i = 0; i < arrayLength; i++)
    {
      T value = array->getValue(i);
      if (value > maxValue)
      {
        maxValue = value;
        maxId = i;
      }
    }

    return maxId;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  template <typename T>
  SizeTArrayType::Pointer findNonZeroIndices(typename DataArray<T>::Pointer array, QVector<size_t> tDims) const
  {
    int size = array->getNumberOfTuples();
    size_t xDim = tDims[0];
    size_t yDim = tDims[1];
    typename DataArray<size_t>::Pointer indices = DataArray<size_t>::CreateArray(size, QVector<size_t>(1, 2), "Non-Zero Indices");

    if (array.get() == nullptr)
      return indices;

    size_t count = 0;
    for (size_t x = 0; x < xDim; x++)
    {
      for (size_t y = 0; y < yDim; y++)
      {
        size_t index = (xDim * y) + x;
        T value = array->getValue(index);
        if (value != 0)
        {
          indices->setComponent(count, 0, x);
          indices->setComponent(count, 1, y);
          count++;
        }
      }
    }

    indices->resize(count);
    return indices;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  SizeTArrayType::Pointer bitwiseMatrixCombination(SizeTArrayType::Pointer matrix1, Int8ArrayType::Pointer matrix2) const
  {
    SizeTArrayType::Pointer result;

    if (matrix1->getNumberOfTuples() != matrix2->getNumberOfTuples() ||
        matrix1->getNumberOfComponents() != matrix2->getNumberOfComponents())
    {
      result = nullptr;
    }
    else
    {
      result = SizeTArrayType::CreateArray(matrix1->getNumberOfTuples(), matrix1->getComponentDimensions(), "Bitwise Matrix Combination");

      for (int y = 0; y < matrix1->getNumberOfTuples(); y++)
      {
        for (int x = 0; x < matrix1->getNumberOfComponents(); x++)
        {
          result->setComponent(y, x, matrix1->getComponent(y, x) & matrix2->getComponent(y, x));
        }
      }
    }

    return result;
  }


  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void analyzeEdgePair(SizeTArrayType::Pointer obj_edge_pair_a1, SizeTArrayType::Pointer obj_edge_pair_b1, size_t index, QVector<size_t> obj_tDims,
                     Int8ArrayType::Pointer obj_mask_edge, size_t& can_num, DoubleArrayType::Pointer cenx_can, DoubleArrayType::Pointer ceny_can, DoubleArrayType::Pointer maj_can, DoubleArrayType::Pointer min_can, DoubleArrayType::Pointer rot_can, DoubleArrayType::Pointer accum_can) const
  {
    const int daxis = m_Axis_Max - m_Axis_Min;
    int axis_min_sq = m_Axis_Min * m_Axis_Min;

    // Matlab is column-based, so x & y are swapped
    double dobj_x = obj_tDims[1];
    double dobj_y = obj_tDims[0];

    double x1 = obj_edge_pair_a1->getComponent(index, 1);
    double y1 = obj_edge_pair_a1->getComponent(index, 0);

    double x2 = obj_edge_pair_b1->getComponent(index, 1);
    double y2 = obj_edge_pair_b1->getComponent(index, 0);

    double x0 = (x1 + x2) / 2;
    double y0 = (y1 + y2) / 2;
    double a = sqrt(std::pow((x2 - x1), 2) + std::pow((y2 - y1), 2)) / 2;
    double alpha = atan2(y2 - y1, x2 - x1);

    SizeTArrayType::Pointer accum = SizeTArrayType::CreateArray(daxis, QVector<size_t>(1, 1), "Accumulator");
    accum->initializeWithZeros();

    if (a >= m_Axis_Min && a <= m_Axis_Max)
    {
      for (int k = 0; k < obj_edge_pair_a1->getNumberOfTuples(); k++)
      {

        double x3 = obj_edge_pair_a1->getComponent(k, 1);
        double y3 = obj_edge_pair_a1->getComponent(k, 0);

        double dsq = std::pow((x3 - x0), 2) + std::pow((y3 - y0), 2);
        double asq = a * a;

        if (dsq > axis_min_sq && dsq < asq)
        {
          double d = sqrt(dsq);
          double fsq = std::pow((x3 - x2), 2) + std::pow((y3 - y2), 2);

          double costau = (asq + dsq - fsq) / (2 * a*d);

          double costau_sq = costau * costau;
          double sintau_sq = 1 - costau_sq;

          double bsq = (asq * dsq * sintau_sq) / (asq - dsq * costau_sq);

          double b = sqrt(bsq);

          //Add one to count from one
          int bidx = static_cast<int>(std::round(b) - m_Axis_Min + 1);

          if (bidx <= daxis && bidx > 0)
          {
            size_t value = accum->getValue(bidx);
            accum->setValue(bidx, value + 1);
          }
        }
      }

      int accum_idx = getIdOfMax<size_t>(accum);
      double accum_max = accum->getValue(accum_idx);

      if (accum_max > 5)
      {
        double b = accum_idx + m_Axis_Min - 1;

        if (b / a > m_Ba_Min)
        {
          //Draw ellipse and compare to object
          size_t count = 0;
          DoubleArrayType::Pointer ellipseCoords = m_Filter->plotEllipsev2(std::round(x0), std::round(y0), std::round(a), std::round(b), alpha, count);

          // create I_check as a 2D array and assign all values to zero
          QVector<size_t> I_check_dims;
          I_check_dims.push_back(dobj_y);
          I_check_dims.push_back(dobj_x);

          SizeTArrayType::Pointer I_check = SizeTArrayType::CreateArray(I_check_dims, QVector<size_t>(1, 1), "I_check");
          I_check->initializeWithZeros();

          for (int k = 0; k < ellipseCoords->getNumberOfTuples(); k++)
          {
            double x = ellipseCoords->getComponent(k, 0);
            double y = ellipseCoords->getComponent(k, 1);

            if (x >= 0 && x < dobj_x && y >= 0 && y < dobj_y)
            {
              size_t index = (I_check_dims[0] * x) + y;
              I_check->setValue(index, 1);

              if (x + 1 < dobj_x)
              {
                size_t index = (I_check_dims[0] * (x + 1)) + y;
                I_check->setValue(index, 1);
              }
              if (x - 1 >= 0)
              {
                size_t index = (I_check_dims[0] * (x - 1)) + y;
                I_check->setValue(index, 1);
              }
              if (y + 1 < dobj_y)
              {
                size_t index = (I_check_dims[0] * x) + (y + 1);
                I_check->setValue(index, 1);
              }
              if (y - 1 >= 0)
              {
                size_t index = (I_check_dims[0] * x) + (y - 1);
                I_check->setValue(index, 1);
              }
              if (x + 1 < dobj_x && y + 1 < dobj_y)
              {
                size_t index = (I_check_dims[0] * (x + 1)) + (y + 1);
                I_check->setValue(index, 1);
              }
              if (x - 1 >= 0 && y + 1 < dobj_y)
              {
                size_t index = (I_check_dims[0] * (x - 1)) + (y + 1);
                I_check->setValue(index, 1);
              }
              if (x + 1 < dobj_x && y - 1 >= 0)
              {
                size_t index = (I_check_dims[0] * (x + 1)) + (y - 1);
                I_check->setValue(index, 1);
              }
              if (x - 1 >= 0 && y - 1 >= 0)
              {
                size_t index = (I_check_dims[0] * (x - 1)) + (y - 1);
                I_check->setValue(index, 1);
              }
            }
          }

          SizeTArrayType::Pointer combinedMatrix = bitwiseMatrixCombination(I_check, obj_mask_edge);
          SizeTArrayType::Pointer overlap = findNonZeroIndices<size_t>(combinedMatrix, I_check_dims);

          // Estimate perimeter length using Ramanujan'a approximation.
          double perim = M_PI * (3 * (a + b) - sqrt((3 * a + b)*(a + 3 * b)));
          // Calculate pixel tolerance based on
          // the calculated perimeter
          double tol_pix = std::round(perim * m_TolEllipse);

          if (overlap->getNumberOfTuples() > tol_pix)
          {
            // Accept point as a new candidate
            cenx_can->setComponent(can_num, 0, std::round(x0)); //x - coordinate of ellipse
            ceny_can->setComponent(can_num, 0, std::round(y0)); //y - coordinate of ellipse
            maj_can->setComponent(can_num, 0, std::round(a)); //major semi - axis
            min_can->setComponent(can_num, 0, std::round(b)); //minor semi - axis
            rot_can->setComponent(can_num, 0, alpha); //Counter clockwise rotation from x - axis
            accum_can->setComponent(can_num, 0, accum_max); //Accumulation matrix
            can_num++;
          }
        }
      }
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  DoubleArrayType::Pointer fillEllipse(DoubleArrayType::Pointer I, QVector<size_t> I_tDims, double xc, double yc, double p, double q, double theta, double val) const
  {
    /* if(theta >= 0) %(xa,xb) is in 1st quadrant and (xb,yb) is in 2nd quadrant
     else (xa,xb) is in 4th quadrant and (xb,yb) is in 1nd quadrant */
    double xa = p * std::cos(theta);
    double ya = p * std::sin(theta);
    double xb = -q * std::sin(theta);
    double yb = q * std::cos(theta);

    double xa_sq = xa * xa;
    double ya_sq = ya * ya;
    double xb_sq = xb * xb;
    double yb_sq = yb * yb;

    double xbyb_sqsq = std::pow((xb_sq + yb_sq), 2);
    double xaya_sqsq = std::pow((xa_sq + ya_sq), 2);

    /* (xa,ya) and (xb,yb) are the points on the ellipse where the major and
     minor axis intersect with the ellipse boundary */

    double a = (xa_sq * xbyb_sqsq) + (xb_sq * xaya_sqsq);
    double b = (xa * ya * xbyb_sqsq) + (xb * yb * xaya_sqsq);
    double c = (ya_sq * xbyb_sqsq) + (yb_sq * xaya_sqsq);
    double d = xaya_sqsq * xbyb_sqsq;

    /* a,b,c,d are the parameters of an ellipse centered at the origin such that
     the ellipse is described by the equation:
     A*x^2 + 2*B*x*y + C*y^2 = D */

    int maxStack = 1000; // Initial stack size

    // Artifical stack
    DoubleArrayType::Pointer stackX = DoubleArrayType::CreateArray(maxStack, QVector<size_t>(1, 1), "stackX");
    stackX->initializeWithZeros();
    DoubleArrayType::Pointer stackY = DoubleArrayType::CreateArray(maxStack, QVector<size_t>(1, 1), "stackY");
    stackY->initializeWithZeros();

    // push back current point to begin search
    stackX->setValue(0, xc);
    stackY->setValue(0, yc);
    int stackSize = 0;

    size_t sizeX = I_tDims[1];
    size_t sizeY = I_tDims[0];

    if(xc <= 0 || xc > sizeX || yc <= 0 || yc > sizeY)
    {
      //Error: Fill ellipse initiated outside of image boundary!
      return DoubleArrayType::NullPointer();
    }

    DoubleArrayType::Pointer I_tmp = DoubleArrayType::CreateArray(I_tDims, QVector<size_t>(1, 1), "I_tmp");
    I_tmp->initializeWithZeros();

    bool copy = I->copyIntoArray(I_tmp);
    if (copy == false)
    {
      //Error: Could not copy array contents into new array!
      return DoubleArrayType::NullPointer();
    }

    while(stackSize > -1)
    {
      double x = stackX->getValue(stackSize);
      double y = stackY->getValue(stackSize);
      stackSize--;

      double xt = x - xc;
      double yt = y - yc;

      double sigma = a*xt*xt+2*b*xt*yt+c*yt*yt-d;

      size_t index = I_tDims[0] * x + y;
      if (sigma <= 0 && I_tmp->getValue(index) != val) // fill in pixel
      {
        I_tmp->setValue(index, val);

        if(stackSize > maxStack - 5)
        {
          // Increase stack size
          stackX->resize(maxStack+1000);
          stackX->initializeWithValue(0, maxStack);

          stackY->resize(maxStack+1000);
          stackY->initializeWithValue(0, maxStack);
          maxStack = maxStack + 1000;
        }

        // x + 1
        if( x + 1 < sizeX )
        {
          stackSize++;
          stackX->setValue(stackSize, x + 1);
          stackY->setValue(stackSize, y);
        }

        // x - 1
        if( x - 1 >= 0 )
        {
          stackSize++;
          stackX->setValue(stackSize, x - 1);
          stackY->setValue(stackSize, y);
        }

        // y + 1
        if( y + 1 < sizeY )
        {
          stackSize++;
          stackX->setValue(stackSize, x);
          stackY->setValue(stackSize, y + 1);
        }

        // y - 1
        if( y - 1 >= 0 )
        {
          stackSize++;
          stackX->setValue(stackSize, x);
          stackY->setValue(stackSize, y - 1);
        }
      }
    }

    return I_tmp;
  }
};

double DetectEllipsoids::m_img_scale_length = 588.0;

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
  m_Ellipse_Count(0),
  m_MaxFeatureId(0),
  m_TotalNumberOfFeatures(0),
  m_FeaturesCompleted(0),
  m_AdditionalEllipsesCount(0),
  m_MaxFeatureIdSem(1),
  m_FeaturesCompletedSem(1),
  m_AdditionalEllipsesCountSem(1)
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

  m_TotalNumberOfFeatures = 0;
  m_FeaturesCompleted = 0;
  m_MaxFeatureId = 0;

  // Initialize counter to track number of detected ellipses
  m_Ellipse_Count = 0;
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

    AttributeMatrix::Pointer featureIdsAM = getDataContainerArray()->getAttributeMatrix(m_FeatureIdsArrayPath);

    QVector<size_t> imageDims = featureIdsAM->getTupleDimensions();
    size_t xDim = imageDims[0], yDim = imageDims[1], zDim = imageDims[2];

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

    double img_pix_length = m_ImageScaleBarLength / m_img_scale_length;
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

    Int32ArrayType::Pointer convOffsetArray = createOffsetArray(orient_tDims);

    int n_size = 3;
    std::vector<double> smoothFil = smoothingFilter(n_size);
    QVector<size_t> smooth_tDims;
    smooth_tDims.push_back(2*n_size+1);
    smooth_tDims.push_back(2*n_size+1);
    smooth_tDims.push_back(1);
    Int32ArrayType::Pointer smoothOffsetArray = createOffsetArray(smooth_tDims);

    m_TotalNumberOfFeatures = activeAM->getNumberOfTuples();
    QString ss = QObject::tr("0/%2").arg(m_TotalNumberOfFeatures);
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

    qint64 millis = QDateTime::currentMSecsSinceEpoch();

    /* Initialize output arrays.  Each array is initialized to m_TotalNumberOfFeatures+1 because the output arrays
       will not use the 0 index so that each output array index matches up to a featureId. */
    DoubleArrayType::Pointer cenx = DoubleArrayType::CreateArray(m_TotalNumberOfFeatures+1, QVector<size_t>(1, 1), "cenx"); // x-coordinate of ellipse
    for (int i=0; i < cenx->getNumberOfTuples(); i++)
    {
      cenx->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer ceny = DoubleArrayType::CreateArray(m_TotalNumberOfFeatures+1, QVector<size_t>(1, 1), "ceny"); // y-coordinate of ellipse
    for (int i=0; i < ceny->getNumberOfTuples(); i++)
    {
      ceny->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer majaxis = DoubleArrayType::CreateArray(m_TotalNumberOfFeatures+1, QVector<size_t>(1, 1), "majaxis"); // major semi-axis
    for (int i=0; i < majaxis->getNumberOfTuples(); i++)
    {
      majaxis->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer minaxis = DoubleArrayType::CreateArray(m_TotalNumberOfFeatures+1, QVector<size_t>(1, 1), "minaxis"); // minor semi-axis
    for (int i=0; i < minaxis->getNumberOfTuples(); i++)
    {
      minaxis->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    DoubleArrayType::Pointer rotangle = DoubleArrayType::CreateArray(m_TotalNumberOfFeatures+1, QVector<size_t>(1, 1), "rotangle"); // Counter clockwise rotation from x-axis
    for (int i=0; i < rotangle->getNumberOfTuples(); i++)
    {
      rotangle->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    }

    m_MaxFeatureId = m_TotalNumberOfFeatures;

    DoubleArrayType::Pointer additionalEllipses = DoubleArrayType::CreateArray(10, QVector<size_t>(1, 6), "additionalEllipses"); // Counter clockwise rotation from x-axis
    for (int i=0; i < additionalEllipses->getNumberOfTuples(); i++)
    {
      additionalEllipses->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
      additionalEllipses->setComponent(i, 1, std::numeric_limits<double>::quiet_NaN());
      additionalEllipses->setComponent(i, 2, std::numeric_limits<double>::quiet_NaN());
      additionalEllipses->setComponent(i, 3, std::numeric_limits<double>::quiet_NaN());
      additionalEllipses->setComponent(i, 4, std::numeric_limits<double>::quiet_NaN());
      additionalEllipses->setComponent(i, 5, std::numeric_limits<double>::quiet_NaN());
    }

    // ***********
    getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, "Test_DC");
    // ***********

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    tbb::task_scheduler_init init;
    bool doParallel = true;

    if(doParallel == true)
    {
      tbb::task_group* g = new tbb::task_group;
      int threads = init.default_num_threads();
      unsigned int numOfTasks = (m_TotalNumberOfFeatures-1) / threads;

      int32_t start = 1;
      int32_t end = 0 + numOfTasks;
      for (int i=0; i<threads; i++)
      {
        g->run(DetectEllipsoidsImpl(this, cellFeatureIdsPtr, imageDims, corners, start, end, m_TotalNumberOfFeatures, convCoords_X, convCoords_Y, convCoords_Z, orient_tDims, convOffsetArray, smoothFil, smoothOffsetArray, axis_min, axis_max, m_HoughTransformThreshold, m_MinAspectRatio, cenx, ceny, majaxis, minaxis, rotangle, additionalEllipses));
        start = end;
        end = end + numOfTasks;
        if(end >= m_TotalNumberOfFeatures)
        {
          end = m_TotalNumberOfFeatures;
        }
      }

      g->wait();
      delete g;
    }
    else
#endif
    {
      DetectEllipsoidsImpl impl(this, cellFeatureIdsPtr, imageDims, corners, 1, m_TotalNumberOfFeatures, m_TotalNumberOfFeatures, convCoords_X, convCoords_Y, convCoords_Z, orient_tDims, convOffsetArray, smoothFil, smoothOffsetArray, axis_min, axis_max, m_HoughTransformThreshold, m_MinAspectRatio, cenx, ceny, majaxis, minaxis, rotangle, additionalEllipses);
      impl();
    }
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
DoubleArrayType::Pointer DetectEllipsoids::plotEllipsev2(double xc, double yc, double p, double q, double theta, size_t &count)
{
  // xc, yc = center of ellipse
  // p, q = length of semi-major and semi-minor axes, respectively
  // theta = angle of counterclockwise rotation of major axis from x-axis in radians

//    if(isreal(xc) == 0 || isreal(yc) == 0 || isreal(p) == 0 || isreal(q) == 0 || isreal(theta) == 0)
//    {
//      // Error: Input must be real valued!
//      return DoubleArrayType::NullPointer();
//    }

  // theta must statisfy: -M_PI/2 < theta <= M_PI/2 (can be rotated due to symmetry)
  while(theta > M_PI/2)
  {
    theta = theta - M_PI;
  }
  while(theta <= -M_PI/2)
  {
    theta = theta + M_PI;
  }

  // if(theta >= 0) %(xa,xb) is in 1st quadrant and (xb,yb) is in 2nd quadrant
  // else (xa,xb) is in 4th quadrant and (xb,yb) is in 1nd quadrant
  double xa = p*cos(theta);
  double ya = p*sin(theta);
  double xb = -q*sin(theta);
  double yb = q*cos(theta);

  double xa_sq = xa * xa;
  double ya_sq = ya * ya;
  double xb_sq = xb * xb;
  double yb_sq = yb * yb;

  double xbyb_sqsq = std::pow((xb_sq + yb_sq), 2);
  double xaya_sqsq = std::pow((xa_sq + ya_sq), 2);

  // (xa,ya) and (xb,yb) are the points on the ellipse where the major and minor
  // axis intersect with the ellipse boundary

  double a = (xa_sq * xbyb_sqsq) + (xb_sq * xaya_sqsq);
  double b = (xa * ya * xbyb_sqsq) + (xb * yb * xaya_sqsq);
  double c = (ya_sq * xbyb_sqsq) + (yb_sq * xaya_sqsq);
  double d = xaya_sqsq * xbyb_sqsq;

  // a,b,c,d are the parameters of an ellipse centered at the origin such that
  // the ellipse is described by the eqution:
  //   A*x^2 + 2*B*x*y + C*y^2 = D

  //Initialize Values
  double y = -ya;
  double x = -xa;
  double dy = -( (a * x) + (b * y) );
  double dx = (b * x) + (c * y);

  // Round values to nearest whole integer
  a=std::round(a);
  b=std::round(b);
  c=std::round(c);
  d=std::round(d);
  x=std::round(x);
  y=std::round(y);
  dx=std::round(dx);
  dy=std::round(dy);

  // estimate number of points on ellipse for array pre-allocation using arc length

  // Estimate perimeter using approximate formula
  // (Note this is a bad approximation if the eccentricity is high)
  size_t perim = static_cast<size_t>(std::ceil( (M_PI * sqrt( 2 * (p*p + q*q) - std::pow((p-q), 2) / 2)) ));
  // Preallocate array using estimated perimeter
  DoubleArrayType::Pointer ellipseCoords = DoubleArrayType::CreateArray(perim, QVector<size_t>(1, 2), "Ellipse Coordinates");
  for (int i=0; i < ellipseCoords->getNumberOfTuples(); i++)
  {
    ellipseCoords->setComponent(i, 0, std::numeric_limits<double>::quiet_NaN());
    ellipseCoords->setComponent(i, 1, std::numeric_limits<double>::quiet_NaN());
  }

  if( x <= 0 && y <= 0 ) // (-xa,-ya) is in the third quadrant or on the x axis
  {
    if( dx == 0 || std::abs(dy/dx) > 1 ) // 1. Slope at (-xa,-ya) is larger than 1 in magnitude. Five sub-arcs are drawn.
    {
      /* (a) Arc from (-xa, -ya) to a point (x0, y0) whose slope is
            infinite. For all points between, the ellipse has slope larger
            than 1 in magnitude, so y is always incremented at each step. */

      while (dx < 0) // loop until point with infinite slope occurs
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y++;
        dy = dy - b;
        dx = dx + c;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          x--;
          dy = dy + a;
          dx = dx - b;
        }
      }

      /* (b) Arc from (x0, y0) to a point (x1, y1) whose slope is 1. For
            all points between, the ellipse has slope larger than 1 in
            magnitude, so y is always incremented at each step. */
      while (dy > dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y++;
        dy = dy - b;
        dx = dx + c;
        double sigma = a*(x+1)*(x+1)+2*b*(x+1)*y+c*y*y-d;
        if(sigma >= 0)
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }

      /* (c) Arc from (x1, y1) to a point (x2, y2) whose slope is 0. For
            all points between, the ellipse has slope less than 1 in
            magnitude, so x is always incremented at each step. */
      while (dy > 0)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          y++;
          dy = dy - b;
          dx = dx + c;
        }
      }

      /* (d) Arc from (x2, y2) to a point (x3, y3) whose slope is -1. For
            all points between, the ellipse has slope less than 1 in
            magnitude, so x is always incremented at each step. */
      while (std::abs(dy) < dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*(y-1)+c*(y-1)*(y-1)-d;
        if(sigma >= 0)
        {
          y--;
          dy = dy + b;
          dx = dx - c;
        }
      }

      /* (e) Arc from (x3, y3) to (xa, ya). For all points between, the
            ellipse has slope larger than 1 in magnitude, so y is always
            decremented at each step. */
      while (y > ya)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y--;
        dy = dy + b;
        dx = dx - c;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }
    }
    else // 2. Slope at (-xa,-ya) is smaller than equal to 1 in magnitude. Five subarcs are drawn (i.e., abs(dy/dx) <= 1 ).
    {
      /* (a) Arc from (-xa, -ya) to a point (x0, y0) whose slope is -1. For
            all points between, the ellipse has slope less than 1 in
            magnitude, so x is always decremented at each step.
          while ( dy < abs(dx) ) % loop until point with infinite slope occurs */

      //Store pixel values
      STORE_PIXEL_VALUES(ellipseCoords, count);

      x++;
      dy = dy + a;
      dx = dx - b;
      double sigma = a*x*x+2*b*x*(y+1)+c*(y+1)*(y+1)-d;
      if(sigma >= 0)
      {
        y++;
        dy = dy - b;
        dx = dx + c;
      }

      /* (b) Arc from (x0, y0) to a point (x1, y1) whose slope is infinite. For
            all points between, the ellipse has slope larger than 1 in
            magnitude, so y is always incremented at each step. */
      while (dx < 0)
      {

        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y++;
        dy = dy - b;
        dx = dx + c;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          x--;
          dy = dy + a;
          dx = dx - b;
        }
      }

      /* (c) Arc from (x1, y1) to a point (x2, y2) whose slope is 1. For
          % all points between, the ellipse has slope larger than 1 in
          % magnitude, so y is always incremented at each step. */
      while (dy > dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y++;
        dy = dy - b;
        dx = dx + c;
        double sigma = a*(x+1)*(x+1)+2*b*(x+1)*y+c*y*y-d;
        if( sigma >= 0 )
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }

      /* (d) Arc from (x2, y2) to a point (x3, y3) whose slope is
            zero. For all points between, the ellipse has slope less
            than 1 in magnitude, so x is always incremented at each step. */

      while (dy > 0) // loop until point with infinite slope occurs
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          y++;
          dy = dy - b;
          dx = dx + c;
        }
      }

      /* (e) Arc from (x3, y3) to (xa, ya). For all points between, the
            ellipse has slope less than 1 in magnitude, so x is always
            incremented at each step. */

      while (x < xa)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*(y-1)+c*(y-1)*(y-1)-d;
        if(sigma >= 0)
        {
          y--;
          dy = dy + b;
          dx = dx - c;
        }
      }
    }
  }
  else // (-xa,-xb) is in the second quadrant
  {
    if ( std::abs(dy/dx) >= 1 ) // 1. Slope at (-xa,-ya) is greater than or equal to 1 in magnitude. Five subarcs are drawn.
    {
      /* (a) Arc from (-xa,-ya) to a point (x0, y0) whose slope is 1.
            For all points between, the ellipse has slope larger than 1 in
            magnitude, so y is incremented at each step. */
      while (dy > dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y++;
        dy = dy - b;
        dx = dx + c;
        double sigma = a*(x+1)*(x+1)+2*b*(x+1)*y+c*y*y-d;
        if(sigma >= 0)
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }

      /* (b) Arc from (x0, y0) to a point (x1, y1) whose slope is
            zero. For all points between, the ellipse has slope less
            than 1 in magnitude, so x is always incremented. */
      while (dy > 0)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          y++;
          dy = dy - b;
          dx = dx + c;
        }
      }

      /* (c) Arc from (x1, y1) to a point (x2, y2) whose slope is -1.
            For all points between, the ellipse has slope less than 1 in
            magnitude, so x is always incremented at each step. */
      while (std::abs(dy) < dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*(y-1)+c*(y-1)*(y-1)-d;
        if(sigma >= 0)
        {
          y--;
          dy = dy + b;
          dx = dx - c;
        }
      }

      /* (d) Arc from (x2, y2) to a point (x3, y3) whose slope is infinity.
          % For all points between, the ellipse has slope greater than 1 in
          % magnitude, so y is always decremented at each step. */
      while (dx > 0)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y--;
        dy = dy + b;
        dx = dx - c;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }

      /* (e) Arc from (x3, y3) to (xa, ya). For all points between, the
          % ellipse has slope greater than 1 in magnitude, so y is always
          % decremented at each step. */
      while (y > ya)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y--;
        dy = dy + b;
        dx = dx - c;
        double sigma = a*(x-1)*(x-1)+2*b*(x-1)*y+c*y*y-d;
        if(sigma >= 0)
        {
          x--;
          dy = dy + a;
          dx = dx - b;
        }
      }
    }

    else // 2. Slop at (-xa,-ya) is smaller than 1 in magnitude (i.e., dy/dx < 0)
    {
      /* (a) Arc from (-xa, -ya) to a point (x0, y0) whose slope is 0.
            For all points between, the ellipse has slope less than 1 in
            magnitude, so x is always incremented at each step. */
      while (dy > 0)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          y++;
          dy = dy - b;
          dx = dx + c;
        }
      }

      /* (b) Arc from (x0,y0) to a point (x1, y1) whose slope is -1.
            For all points between, the ellipse has slope less than 1 in
            magnitude, so x is always decremented. */
      while (std::abs(dy) < dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x++;
        dy = dy - a;
        dx = dx + b;
        double sigma = a*x*x+2*b*x*(y-1)+c*(y-1)*(y-1)-d;
        if(sigma >= 0)
        {
          y--;
          dy = dy + b;
          dx = dx - c;
        }
      }

      /* (c) Arc from (x1, y1) to a point (x2, y2) whose slope is
            infinite. For all points between, the ellipse has slope larger
            than 1, so y is always incremented. */
      while (dx > 0)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y--;
        dy = dy + b;
        dx = dx - c;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          x++;
          dy = dy - a;
          dx = dx + b;
        }
      }

      /* (d) Arc from (x2, y2) to a point (x3, y3) whose slope is 1.
            For all points between, the ellipse has slope larger than 1 in
            magnitude, so y is always decremented at each step. */
      while (dy < dx)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        y--;
        dy = dy + b;
        dx = dx - c;
        double sigma = a*(x-1)*(x-1)+2*b*(x-1)*y+c*y*y-d;
        if(sigma >= 0)
        {
          x--;
          dy = dy + a;
          dx = dx - b;
        }
      }

      /* (e) Arc from (x3, y3) to (xa, ya). For all points between, the
          % ellipse has slope less than 1 in magnitude, so x is always
          % incremented at each step. */
      while (x > xa)
      {
        //Store pixel values
        STORE_PIXEL_VALUES(ellipseCoords, count);

        x--;
        dy = dy + a;
        dx = dx - b;
        double sigma = a*x*x+2*b*x*y+c*y*y-d;
        if(sigma < 0)
        {
          y--;
          dy = dy + b;
          dx = dx - c;
        }
      }

    }
  }

  count--;

  return ellipseCoords;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
size_t DetectEllipsoids::getUniqueFeatureId()
{
  m_MaxFeatureIdSem.acquire();
  m_MaxFeatureId++;
  size_t id = m_MaxFeatureId;
  m_MaxFeatureIdSem.release();
  return id;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
size_t DetectEllipsoids::getAdditionalEllipsesIndex()
{
  m_AdditionalEllipsesCountSem.acquire();
  size_t index = m_AdditionalEllipsesCount;
  m_AdditionalEllipsesCount++;
  m_AdditionalEllipsesCountSem.release();
  return index;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoids::notifyFeatureCompleted()
{
  m_FeaturesCompletedSem.acquire();
  m_FeaturesCompleted++;
  QString ss = QObject::tr("%1/%2").arg(m_FeaturesCompleted).arg(m_TotalNumberOfFeatures);
  notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
  m_FeaturesCompletedSem.release();
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

