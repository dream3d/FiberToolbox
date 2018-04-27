#pragma once

#include "FiberToolbox/FiberToolboxPlugin.h"

class FiberToolboxGuiPlugin : public FiberToolboxPlugin
{
  Q_OBJECT
  Q_INTERFACES(ISIMPLibPlugin)
  Q_PLUGIN_METADATA(IID "net.bluequartz.dream3d.FiberToolboxGuiPlugin")

public:
  FiberToolboxGuiPlugin();
   ~FiberToolboxGuiPlugin() override;
  
  /**
   * @brief Register all the filters with the FilterWidgetFactory
   */
  void registerFilterWidgets(FilterWidgetManager* fwm) override;
  

public:
  FiberToolboxGuiPlugin(const FiberToolboxGuiPlugin&) = delete;            // Copy Constructor Not Implemented
  FiberToolboxGuiPlugin(FiberToolboxGuiPlugin&&) = delete;                 // Move Constructor
  FiberToolboxGuiPlugin& operator=(const FiberToolboxGuiPlugin&) = delete; // Copy Assignment Not Implemented
  FiberToolboxGuiPlugin& operator=(FiberToolboxGuiPlugin&&) = delete;      // Move Assignment Not Implemented
};
