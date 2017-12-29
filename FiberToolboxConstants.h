/*
 * Your License should go here
 */
#ifndef _fibertoolboxconstants_h_
#define _fibertoolboxconstants_h_

#include <QtCore/QString>

/**
* @brief This namespace is used to define some Constants for the plugin itself.
*/
namespace FiberToolboxConstants
{
  const QString FiberToolboxPluginFile("FiberToolboxPlugin");
  const QString FiberToolboxPluginDisplayName("FiberToolboxPlugin");
  const QString FiberToolboxBaseName("FiberToolbox");

  namespace FilterGroups
  {
  	const QString FiberToolboxFilters("FiberToolbox");
  }
}

/**
* @brief Use this namespace to define any custom GUI widgets that collect FilterParameters
* for a filter. Do NOT define general reusable widgets here.
*/
namespace FilterParameterWidgetType
{
/* const QString SomeCustomWidget("SomeCustomWidget"); */
}
#endif
