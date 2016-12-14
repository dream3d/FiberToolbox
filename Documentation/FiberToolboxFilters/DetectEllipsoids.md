Detect Ellipsoids {#detectellipsoids}
=============

## Group (Subgroup) ##
FiberToolbox (FiberToolbox)

## Description ##
This **Filter** detects ellipsoids in an existing *FeatureIds* array, and outputs a new FeatureId array that contains all the ellipses that were found.

The input *FeatureIds* array 

## Parameters ##
| Name | Type | Description |
|------|------|------|
| MinFiberAxisLength | Integer | The minimum length of the fiber axis |
| MaxFiberAxisLength | Integer | The maximum length of the fiber axis |
| HoughTransformThreshold | Double | Threshold for Hough Transform |
| MinAspectRatio | Double | Minimum Aspect Ratio |
| ImageScaleBarLength | Integer | Length of the Image Scale Bar |
| ImageScaleBarUnits | Choice | Units used in calculations |

## Required Geometry ##
Image

## Required Objects ##
| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| Cell **Attribute Array** | FeatureIds | int32_t | (1) | The Feature Ids array to analyze for ellipsoids |
| CellFeature **Attribute Array** | Active | bool | (1) | The Feature Ids that are active in the Feature Ids array |

## Created Objects ##
| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| Cell **Attribute Array** | DetectedEllipsoidsFeatureIdsArrayPath | int32_t | (1) | The path to the Feature Ids array that contains detected ellipsoids |

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users