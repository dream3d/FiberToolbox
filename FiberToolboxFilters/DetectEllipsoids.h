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

#ifndef _fibertoolboxfilter_h_
#define _fibertoolboxfilter_h_

#include <complex>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"

typedef std::vector<std::complex<double> > DE_ComplexDoubleVector;

/**
 * @brief The DetectEllipsoids class. See [Filter documentation](@ref fibertoolboxfilter) for details.
 */
class DetectEllipsoids : public AbstractFilter
{
  Q_OBJECT

  public:
    SIMPL_SHARED_POINTERS(DetectEllipsoids)
    SIMPL_STATIC_NEW_MACRO(DetectEllipsoids)
    SIMPL_TYPE_MACRO_SUPER(DetectEllipsoids, AbstractFilter)

    virtual ~DetectEllipsoids();

    enum ScaleBarUnits
    {
      MillimeterUnits = 0,
      MicronUnits = 1
    };

    SIMPL_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
    Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

    SIMPL_FILTER_PARAMETER(DataArrayPath, ActiveArrayPath)
    Q_PROPERTY(DataArrayPath ActiveArrayPath READ getActiveArrayPath WRITE setActiveArrayPath)

    SIMPL_FILTER_PARAMETER(int, MinFiberAxisLength)
    Q_PROPERTY(int MinFiberAxisLength READ getMinFiberAxisLength WRITE setMinFiberAxisLength)

    SIMPL_FILTER_PARAMETER(int, MaxFiberAxisLength)
    Q_PROPERTY(int MaxFiberAxisLength READ getMaxFiberAxisLength WRITE setMaxFiberAxisLength)

    SIMPL_FILTER_PARAMETER(float, HoughTransformThreshold)
    Q_PROPERTY(float HoughTransformThreshold READ getHoughTransformThreshold WRITE setHoughTransformThreshold)

    SIMPL_FILTER_PARAMETER(float, MinAspectRatio)
    Q_PROPERTY(float MinAspectRatio READ getMinAspectRatio WRITE setMinAspectRatio)

    SIMPL_FILTER_PARAMETER(int, ImageScaleBarLength)
    Q_PROPERTY(int ImageScaleBarLength READ getImageScaleBarLength WRITE setImageScaleBarLength)

    SIMPL_FILTER_PARAMETER(int, ImageScaleBarUnits)
    Q_PROPERTY(int ImageScaleBarUnits READ getImageScaleBarUnits WRITE setImageScaleBarUnits)

    /**
     * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getCompiledLibraryName();

    /**
     * @brief getBrandingString Returns the branding string for the filter, which is a tag
     * used to denote the filter's association with specific plugins
     * @return Branding string
    */
    virtual const QString getBrandingString();

    /**
     * @brief getFilterVersion Returns a version string for this filter. Default
     * value is an empty string.
     * @return
     */
    virtual const QString getFilterVersion();

    /**
     * @brief newFilterInstance Reimplemented from @see AbstractFilter class
     */
    virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);

    /**
     * @brief getGroupName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getGroupName();

    /**
     * @brief getSubGroupName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getSubGroupName();

    /**
     * @brief getHumanLabel Reimplemented from @see AbstractFilter class
     */
    virtual const QString getHumanLabel();

    /**
     * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
     */
    virtual void setupFilterParameters();

    /**
     * @brief execute Reimplemented from @see AbstractFilter class
     */
    virtual void execute();

    /**
    * @brief preflight Reimplemented from @see AbstractFilter class
    */
    virtual void preflight();

  signals:
    /**
     * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
     * be pushed from a user-facing control (such as a widget)
     * @param filter Filter instance pointer 
     */
    void updateFilterParameters(AbstractFilter* filter);

    /**
     * @brief parametersChanged Emitted when any Filter parameter is changed internally
     */
    void parametersChanged();

    /**
     * @brief preflightAboutToExecute Emitted just before calling dataCheck()
     */
    void preflightAboutToExecute();

    /**
     * @brief preflightExecuted Emitted just after calling dataCheck()
     */
    void preflightExecuted();

  protected:
    DetectEllipsoids();

    /**
    * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
    */
    void dataCheck();

    /**
    * @brief Initializes all the private instance variables.
    */
    void initialize();

  private:
    static double img_scale_length;

    /**
     * @brief orientationFilter
     * @return
     */
    DoubleArrayType::Pointer orientationFilter(int minAxisLength, int maxAxisLength);

    /**
     * @brief houghCircleFilter
     * @return
     */
    DE_ComplexDoubleVector houghCircleFilter(int minAxisLength, int maxAxisLength);

    DetectEllipsoids(const DetectEllipsoids&); // Copy Constructor Not Implemented
    void operator=(const DetectEllipsoids&); // Operator '=' Not Implemented
};

#endif /* _DetectEllipsoids_H_ */
