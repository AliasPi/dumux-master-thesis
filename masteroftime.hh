// masteroftime.hh
#include <iostream>
#include <optional>
#include <vector>
#include <algorithm> // for std::clamp

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

namespace Dumux {

template <class TypeTag, class Scalar>
class MasterOfTime {
    using Problem = GetPropType<TypeTag, Properties::Problem>;
public:
    MasterOfTime(TimeLoop<Scalar>& timeLoop, Problem& problem)
        : timeLoop_(timeLoop)
        , problem_(problem)
    {
        std::cout << "MasterOfTime Constructor" << std::endl;

        stepTimeIncreaseFactor_ = getParam<Scalar>("TimeLoop.StepTimeIncreaseFactor");
        maxTimeStepSize_ = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

        useWell1_            = getParam<bool>("Well1.UseWell");
        startTimeWell1_      = getParam<Scalar>("Well1.StartTime");
        endTimeWell1_        = getParam<Scalar>("Well1.EndTime");
        dtInitialStartWell1_ = getParam<Scalar>("Well1.DtInitialStart");
        dtInitialEndWell1_   = getParam<Scalar>("Well1.DtInitialEnd");

        useWell2_            = getParam<bool>("Well2.UseWell");
        startTimeWell2_      = getParam<Scalar>("Well2.StartTime");
        endTimeWell2_        = getParam<Scalar>("Well2.EndTime");
        dtInitialStartWell2_ = getParam<Scalar>("Well2.DtInitialStart");
        dtInitialEndWell2_   = getParam<Scalar>("Well2.DtInitialEnd");

        // --- New Parameters for overriding max time step ---
        // Read the new max time step range as two percentages (e.g. 50 60)
        std::vector<Scalar> newMaxRange = getParam<std::vector<Scalar>>("TimeLoop.NewMaxTimeStepSizeRange");
        if(newMaxRange.size() == 2) {
            newMaxTimeStepRangeLower_ = newMaxRange[0];
            newMaxTimeStepRangeUpper_ = newMaxRange[1];
        }
        else {
            // If not properly provided, use the full range [0,100]
            newMaxTimeStepRangeLower_ = 0;
            newMaxTimeStepRangeUpper_ = 100;
        }
        newMaxTimeStepSize_ = getParam<Scalar>("TimeLoop.NewMaxTimeStepSize");
        // We also need the total simulation time (TEnd) to compute the percent of elapsed time.
        simulationEndTime_ = getParam<Scalar>("TimeLoop.TEnd");
        // --------------------------------------------------------
    };

    void advanceTimeStep() {
        // Advance the time loop and report step statistics.
        timeLoop_.advanceTimeStep();
        timeLoop_.reportTimeStep();

        // Compute an effective maximum time step value based on the current simulation time.
        // (We assume that the simulation time is between 0 and simulationEndTime_.)
        Scalar currentTime = time();
        Scalar currentPercent = (currentTime / simulationEndTime_) * 100;
        Scalar effectiveMax = (currentPercent >= newMaxTimeStepRangeLower_ && currentPercent < newMaxTimeStepRangeUpper_) 
                                  ? newMaxTimeStepSize_ 
                                  : maxTimeStepSize_;

        Scalar nextTimeStepSize;

        std::cout << "Current Time: " << currentTime << std::endl;

        if (cachedTimeStepSize_.has_value()){
            std::cout << "------------------------ Loading Cached Step Size ------------------------" << std::endl;
            nextTimeStepSize = *cachedTimeStepSize_;
            cachedTimeStepSize_.reset();
            std::cout << "Time Step Size: " << nextTimeStepSize << std::endl;
            std::cout << "Next Cached Time: " << currentTime + nextTimeStepSize << std::endl; 
        }
        else {
            nextTimeStepSize = timeLoop_.timeStepSize() * stepTimeIncreaseFactor_;
        }

        Scalar nextTime = currentTime + nextTimeStepSize;

        // Check if we are approaching or leaving a well injection interval.
        // In each case use std::clamp with the effective maximum (which may be reduced in the specified percent range).
        if (nextTime >= startTimeWell1_ && nextTime < endTimeWell1_ && !isWellInjecting1() && useWell1_){
            std::cout << "------------------------ Time Step Size Adapted / Injection Well 1 Starting ------------------------" << std::endl;         
            cachedTimeStepSize_ = dtInitialStartWell1_;
            nextTimeStepSize = std::clamp(startTimeWell1_ - currentTime, 1e-5, effectiveMax);
            std::cout << "Time Step Size: " << nextTimeStepSize << std::endl;
            std::cout << "Injection Start Time: " << startTimeWell1_ << std::endl;
            problem_.setWellInjecting1(true);
        }
        if (nextTime >= startTimeWell2_ && nextTime < endTimeWell2_ && !isWellInjecting2() && useWell2_){
            std::cout << "------------------------ Time Step Size Adapted / Injection Well 2 Starting ------------------------" << std::endl;         
            cachedTimeStepSize_ = dtInitialStartWell2_;
            nextTimeStepSize = std::clamp(startTimeWell2_ - currentTime, 1e-5, effectiveMax);
            std::cout << "Time Step Size: " << nextTimeStepSize << std::endl;
            std::cout << "Injection Start Time: " << startTimeWell2_ << std::endl;
            problem_.setWellInjecting2(true);
        }
        if (nextTime >= endTimeWell1_ && isWellInjecting1() && useWell1_){
            std::cout << "------------------------ Time Step Size Adapted / Injection Well 1 Stopping ------------------------" << std::endl;         
            cachedTimeStepSize_ = dtInitialEndWell1_;
            nextTimeStepSize = std::clamp(endTimeWell1_ - currentTime, 1e-5, effectiveMax);
            std::cout << "Time Step Size: " << nextTimeStepSize << std::endl;
            std::cout << "Injection End Time: " << endTimeWell1_ << std::endl;
            problem_.setWellInjecting1(false);
        }
        if (nextTime >= endTimeWell2_ && isWellInjecting2() && useWell2_){
            std::cout << "------------------------ Time Step Size Adapted / Injection Well 2 Stopping ------------------------" << std::endl;         
            cachedTimeStepSize_ = dtInitialEndWell2_;
            nextTimeStepSize = std::clamp(endTimeWell2_ - currentTime, 1e-5, effectiveMax);
            std::cout << "Time Step Size: " << nextTimeStepSize << std::endl;
            std::cout << "Injection End Time: " << endTimeWell2_ << std::endl;
            problem_.setWellInjecting2(false);
        }

        // If no injection override is active, ensure the new time step does not exceed the effective maximum.
        if (!cachedTimeStepSize_.has_value()){
            nextTimeStepSize = std::min(nextTimeStepSize, effectiveMax);
        }

        timeLoop_.setTimeStepSize(nextTimeStepSize);
    }

    bool isWellInjecting1()
    { return time() >= startTimeWell1_ && time() < endTimeWell1_; }

    bool isWellInjecting2()
    { return time() >= startTimeWell2_ && time() < endTimeWell2_; }

    Scalar time()
    { return timeLoop_.time(); }

private:
    Scalar stepTimeIncreaseFactor_;
    Scalar maxTimeStepSize_;

    bool useWell1_, useWell2_;
    Scalar startTimeWell1_, endTimeWell1_, startTimeWell2_, endTimeWell2_;
    Scalar dtInitialStartWell1_, dtInitialEndWell1_, dtInitialStartWell2_, dtInitialEndWell2_;

    // --- New members for alternate maximum time step during a given simulation percentage range ---
    Scalar newMaxTimeStepSize_;
    Scalar newMaxTimeStepRangeLower_;
    Scalar newMaxTimeStepRangeUpper_;
    Scalar simulationEndTime_;
    // ---------------------------------------------------------------------------------------------

    std::optional<Scalar> cachedTimeStepSize_;

    TimeLoop<Scalar>& timeLoop_;
    Problem& problem_;
};

} // end namespace Dumux
