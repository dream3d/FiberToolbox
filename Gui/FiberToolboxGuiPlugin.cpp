

#include "FiberToolboxGuiPlugin.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FiberToolboxGuiPlugin::FiberToolboxGuiPlugin()
: FiberToolboxPlugin()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FiberToolboxGuiPlugin::~FiberToolboxGuiPlugin() = default;

#include "FiberToolbox/Gui/FilterParameterWidgets/RegisterKnownFilterParameterWidgets.cpp"
