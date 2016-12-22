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

#include "DetectEllipsoidsImpl.h"

#include "FiberToolbox/HelperClasses/ComputeGradient.h"
#include "FiberToolbox/FiberToolboxFilters/DetectEllipsoids.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectEllipsoidsImpl::DetectEllipsoidsImpl(DetectEllipsoids* filter, int* cellFeatureIdsPtr, QVector<size_t> cellFeatureIdsDims, UInt32ArrayType::Pointer corners, int32_t featureIdStart, int32_t featureIdEnd, DE_ComplexDoubleVector convCoords_X, DE_ComplexDoubleVector convCoords_Y, DE_ComplexDoubleVector convCoords_Z, QVector<size_t> kernel_tDims, Int32ArrayType::Pointer convOffsetArray, std::vector<double> smoothFil, Int32ArrayType::Pointer smoothOffsetArray, double axis_min, double axis_max, float tol_ellipse, float ba_min, DoubleArrayType::Pointer center, DoubleArrayType::Pointer majaxis, DoubleArrayType::Pointer minaxis, DoubleArrayType::Pointer rotangle, AttributeMatrix::Pointer ellipseFeatureAM) :
  m_Filter(filter),
  m_CellFeatureIdsPtr(cellFeatureIdsPtr),
  m_CellFeatureIdsDims(cellFeatureIdsDims),
  m_Corners(corners),
  m_FeatureIdStart(featureIdStart),
  m_FeatureIdEnd(featureIdEnd),
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
  m_Center(center),
  m_Majaxis(majaxis),
  m_Minaxis(minaxis),
  m_Rotangle(rotangle),
  m_EllipseFeatureAM(ellipseFeatureAM)
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectEllipsoidsImpl::~DetectEllipsoidsImpl()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DetectEllipsoidsImpl::operator()() const
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

  // Run the ellipse detection algorithm on each object
  for(size_t featureId = m_FeatureIdStart; featureId < m_FeatureIdEnd; featureId++)
  {
    size_t topL_X = m_Corners->getComponent(featureId, 0);
    size_t topL_Y = m_Corners->getComponent(featureId, 1);
    size_t topL_Z = m_Corners->getComponent(featureId, 2);
    size_t bottomR_X = m_Corners->getComponent(featureId, 3);
    size_t bottomR_Y = m_Corners->getComponent(featureId, 4);
    size_t bottomR_Z = m_Corners->getComponent(featureId, 5);

    // Calculate the object's dimensions
    size_t obj_xDim = bottomR_X - topL_X + 1;
    size_t obj_yDim = bottomR_Y - topL_Y + 1;

    // Calculate the object's dimensions with a 1-pixel border around it
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

    // Copy object pixels from m_CellFeatureIdsPtr into featureObjArray
    QVector<size_t> cDims(1, 1);
    DoubleArrayType::Pointer featureObjArray = DoubleArrayType::CreateArray(image_tDims, cDims, "featureObjArray");
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

    size_t numberOfDetectedEllipses = 0;
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
        size_t obj_ext_y = (obj_ext_index % image_xDim);
        size_t obj_ext_x = (((obj_ext_index / image_xDim) % image_yDim));

        size_t mask_rad = m_Axis_Max + m_Axis_Min;

        int mask_min_x = obj_ext_x - mask_rad + 1;
        if ( mask_min_x < 1)
        {
          mask_min_x = 1;
        }

        int mask_min_y = obj_ext_y - mask_rad + 1;
        if ( mask_min_y < 1)
        {
          mask_min_y = 1;
        }

        int mask_max_x = obj_ext_x + mask_rad + 1;
        if ( mask_max_x > image_tDims[1])
        {
          mask_max_x = image_tDims[1];
        }

        int mask_max_y = obj_ext_y + mask_rad + 1;
        if ( mask_max_y > image_tDims[0])
        {
          mask_max_y = image_tDims[0];
        }

        DoubleArrayType::Pointer obj_mask = DoubleArrayType::CreateArray(image_tDims, QVector<size_t>(1, 1), "obj_mask");
        obj_mask->initializeWithZeros();

        for (size_t y = mask_min_y - 1; y < mask_max_y; y++)
        {
          for (size_t x = mask_min_x - 1; x < mask_max_x; x++)
          {
            int index = image_tDims[0] * x + y;
            obj_mask->setValue(index, featureObjArray->getValue(index));
          }
        }

        // Compute the edge array
        Int8ArrayType::Pointer edgeArray = findEdges<double>(obj_mask, image_tDims);

        SizeTArrayType::Pointer obj_edges = findNonZeroIndices<int8_t>(edgeArray, image_tDims);

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
          analyzeEdgePair(obj_edge_pair_a1, obj_edge_pair_b1, k, image_tDims, edgeArray, can_num, cenx_can, ceny_can, maj_can, min_can, rot_can, accum_can);
        }

        if(can_num > 0) // Assume best match is the ellipse
        {
          // Increment the ellipse counter
          m_Filter->setEllipse_Count(m_Filter->getEllipse_Count() + 1);

          int accum_idx = getIdOfMax<double>(accum_can);

          /* If this is another ellipse in the same overall object,
             * create a new feature id and resize our output arrays */
          size_t objId = featureId;
          if (numberOfDetectedEllipses > 0)
          {
            objId = m_Filter->getUniqueFeatureId();
            m_EllipseFeatureAM->resizeAttributeArrays(QVector<size_t>(1, objId + 1));
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

          m_Center->setComponent(objId, 0, cenx_val);
          m_Center->setComponent(objId, 1, ceny_val);
          m_Majaxis->setValue(objId, majaxis_val);
          m_Minaxis->setValue(objId, minaxis_val);
          m_Rotangle->setValue(objId, rotangle_val);

          /* Remove pixels from object (including additional 1 pixel thick boarder)
             * and reassign remaining pixels to array */
          Int32ArrayType::Pointer featureObjOnesArray = Int32ArrayType::CreateArray(image_tDims, QVector<size_t>(1, 1), "featureObjOnesArray");
          featureObjOnesArray->initializeWithValue(1);

          Int32ArrayType::Pointer I_tmp = m_Filter->fillEllipse(featureObjOnesArray, image_tDims, cenx_val, ceny_val, majaxis_val+1, minaxis_val+1, rotangle_val, 0);

          for (int i = 0; i < I_tmp->getNumberOfTuples(); i++)
          {
            double value = featureObjArray->getValue(i) * I_tmp->getValue(i);
            featureObjArray->setValue(i, value);
          }

          // Translate center of ellipse into real space
          size_t obj_x_min = topL_Y;
          size_t obj_y_min = topL_X;

          m_Center->setComponent(objId, 0, m_Center->getComponent(objId, 0) + obj_x_min);
          m_Center->setComponent(objId, 1, m_Center->getComponent(objId, 1) + obj_y_min);

          // Clean Accumulator
          accum_can->initializeWithZeros();

          objPixelsArray = findNonZeroIndices<double>(featureObjArray, image_tDims);

          numberOfDetectedEllipses++;
          break;
        }
      }

      if(detectedObjIdx == obj_ext_num)
      {
        break;
      }

      if (m_Filter->getCancel())
      {
        return;
      }
    }

    m_Filter->notifyFeatureCompleted();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QList<int> DetectEllipsoidsImpl::findExtrema(DoubleArrayType::Pointer thresholdArray, QVector<size_t> tDims) const
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

  QList<int> extremaList = extrema.toList();
  std::sort(extremaList.begin(), extremaList.end());

  return extremaList;
}

// -----------------------------------------------------------------------------
// helper method - grabs index from matrix size
// -----------------------------------------------------------------------------
size_t DetectEllipsoidsImpl::sub2ind(QVector<size_t> tDims, size_t row, size_t col) const
{
  return row * tDims[0] + col;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QPair<size_t,size_t> DetectEllipsoidsImpl::plotlineEdgeInter(int x0, int y0, int x1, int y1, DoubleArrayType::Pointer binImage, QVector<size_t> imageDims) const
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
SizeTArrayType::Pointer DetectEllipsoidsImpl::bitwiseMatrixCombination(SizeTArrayType::Pointer matrix1, Int8ArrayType::Pointer matrix2) const
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
void DetectEllipsoidsImpl::analyzeEdgePair(SizeTArrayType::Pointer obj_edge_pair_a1, SizeTArrayType::Pointer obj_edge_pair_b1, size_t index, QVector<size_t> obj_tDims, Int8ArrayType::Pointer obj_mask_edge, size_t& can_num, DoubleArrayType::Pointer cenx_can, DoubleArrayType::Pointer ceny_can, DoubleArrayType::Pointer maj_can, DoubleArrayType::Pointer min_can, DoubleArrayType::Pointer rot_can, DoubleArrayType::Pointer accum_can) const
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
        int bidx = static_cast<int>(std::round(b) - m_Axis_Min);

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
      double b = accum_idx + m_Axis_Min;

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
