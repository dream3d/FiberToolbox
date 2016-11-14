/*
 * Your License or Copyright can go here
 */

#include "DetectEllipsoids.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"


#include "FiberToolbox/FiberToolboxConstants.h"
#include "FiberToolbox/FiberToolboxVersion.h"

// Include the MOC generated file for this class
#include "moc_DetectEllipsoids.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DetectEllipsoids::DetectEllipsoids() :
  AbstractFilter()
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

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, DetectEllipsoids, req));
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

  if (getCancel() == true) { return; }

  if (getErrorCondition() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
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

