/*
 * Your License or Copyright can go here
 */

#include "DetectEllipsoids.h"

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
#include <tbb/atomic.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>
#endif

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "FiberToolbox/FiberToolboxConstants.h"
#include "FiberToolbox/FiberToolboxVersion.h"

// Include the MOC generated file for this class
#include "moc_DetectEllipsoids.cpp"

/**
 * @brief The DetectEllipsoidsImpl class implements a threaded algorithm that detects ellipsoids in a FeatureIds array
 */
class DetectEllipsoidsImpl
{
  int*                            m_CellFeatureIdsPtr;
  size_t                          m_CellFeatureIdsDims[3];
  UInt32ArrayType::Pointer        m_Corners;
  int32_t                         m_FeatureIdStart;
  int32_t                         m_FeatureIdEnd;

public:
  DetectEllipsoidsImpl(int* cellFeatureIdsPtr, size_t cellFeatureIdsDims[3], UInt32ArrayType::Pointer corners, int32_t featureIdStart, int32_t featureIdEnd) :
    m_CellFeatureIdsPtr(cellFeatureIdsPtr),
    m_Corners(corners),
    m_FeatureIdStart(featureIdStart),
    m_FeatureIdEnd(featureIdEnd)
  {
    m_CellFeatureIdsDims[0] = cellFeatureIdsDims[0];
    m_CellFeatureIdsDims[1] = cellFeatureIdsDims[1];
    m_CellFeatureIdsDims[2] = cellFeatureIdsDims[2];
  }

  virtual ~DetectEllipsoidsImpl()
  {
  }

  void operator()() const
  {
    for(int i = m_FeatureIdStart; i < m_FeatureIdEnd; i++)
    {
      size_t topL_X = m_Corners->getComponent(i, 0);
      size_t topL_Y = m_Corners->getComponent(i, 1);
      size_t topL_Z = m_Corners->getComponent(i, 2);
      size_t bottomR_X = m_Corners->getComponent(i, 3);
      size_t bottomR_Y = m_Corners->getComponent(i, 4);
      size_t bottomR_Z = m_Corners->getComponent(i, 5);

      size_t zDim = bottomR_Z - topL_Z + 1;
      size_t yDim = bottomR_Y - topL_Y + 1;
      size_t xDim = bottomR_X - topL_X + 1;

      QVector<size_t> tDims(zDim, xDim * yDim);
      QVector<size_t> cDims(1, 1);
      Int32ArrayType::Pointer featureObjArray = Int32ArrayType::CreateArray(tDims, cDims, "Feature Object");

      for (size_t z = topL_Z; z <= bottomR_Z; z++)
      {
        for (size_t y = topL_Y; y <= bottomR_Y; y++)
        {
          for (size_t x = topL_X; x <= bottomR_X; x++)
          {
            size_t featureX = x - topL_X;
            size_t featureY = y - topL_Y;
            size_t featureZ = z - topL_Z;
            size_t featureObjIndex = (yDim * xDim * featureZ) + (xDim * featureY) + featureX;
            size_t cellFeatureIdsIndex = (m_CellFeatureIdsDims[1] * m_CellFeatureIdsDims[0] * z) + (m_CellFeatureIdsDims[0] * y) + x;
            int32_t value = m_CellFeatureIdsPtr[cellFeatureIdsIndex];
            featureObjArray->setValue(featureObjIndex, value);
          }
        }
      }

      //Now send featureObjArray to craig's code
    }
  }
};


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
  m_ImageScaleBarUnits(ScaleBarUnits::MicronUnits)
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
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, DetectEllipsoids, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, SIMPL::AttributeMatrixType::CellFeature, SIMPL::GeometryType::ImageGeometry);
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

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    tbb::task_scheduler_init init;
    bool doParallel = true;
#endif

    size_t totalNumOfFeatures = activeAM->getNumberOfTuples();

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::task_group* g = new tbb::task_group;
      int threads = init.default_num_threads();
      unsigned int numOfTasks = totalNumOfFeatures / threads;

      int32_t start = 0;
      int32_t end = 0 + numOfTasks;
      for (int i=0; i<threads; i++)
      {
        g->run(DetectEllipsoidsImpl(cellFeatureIdsPtr, dims, corners, start, end));
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
      DetectEllipsoidsImpl impl(cellFeatureIdsPtr, dims, corners, 0, totalNumOfFeatures);
      impl();
    }
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
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

