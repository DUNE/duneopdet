#ifndef DUNEOPDET_OPTICALDETECTOR_ANA_OPDETWAVEFORMTIMESTAMPHELPER_H
#define DUNEOPDET_OPTICALDETECTOR_ANA_OPDETWAVEFORMTIMESTAMPHELPER_H

#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include <cmath>
#include <cstdint>
#include <limits>

namespace opdet {

  class OpDetWaveformTimeStampHelper {
  public:
    using TriggerTimeStamp_t = std::uint64_t;
    using TriggerTimeDiff_t = std::int64_t;
    using WaveformTimeStamp_t = raw::TimeStamp_t;

    static constexpr unsigned int TrimmedBits = 40U;
    static constexpr TriggerTimeStamp_t TrimmedMask =
      (TriggerTimeStamp_t{1} << TrimmedBits) - TriggerTimeStamp_t{1};
    static constexpr TriggerTimeStamp_t TrimmedModulus =
      TriggerTimeStamp_t{1} << TrimmedBits;

    static TriggerTimeStamp_t WaveformToTriggerTicks(
      raw::OpDetWaveform const& waveform,
      double opticalClockTickPeriod)
    {
      return WaveformToTriggerTicks(waveform.TimeStamp(), opticalClockTickPeriod);
    }

    static TriggerTimeStamp_t WaveformToTriggerTicks(
      raw::OpDetWaveform const& waveform,
      detinfo::DetectorClocksData const& clockData)
    {
      return WaveformToTriggerTicks(waveform, clockData.OpticalClock().TickPeriod());
    }

    static TriggerTimeStamp_t WaveformToTriggerTicks(
      WaveformTimeStamp_t waveformTimeStamp,
      double opticalClockTickPeriod)
    {
      return RoundToTicks(waveformTimeStamp, opticalClockTickPeriod) & TrimmedMask;
    }

    static TriggerTimeStamp_t WaveformToTriggerTicks(
      WaveformTimeStamp_t waveformTimeStamp,
      detinfo::DetectorClocksData const& clockData)
    {
      return WaveformToTriggerTicks(
        waveformTimeStamp, clockData.OpticalClock().TickPeriod());
    }

    static WaveformTimeStamp_t TriggerToWaveformTimeStamp(
      TriggerTimeStamp_t triggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return static_cast<WaveformTimeStamp_t>(triggerTimeStamp & TrimmedMask) *
             opticalClockTickPeriod;
    }

    static WaveformTimeStamp_t TriggerToWaveformTimeStamp(
      TriggerTimeStamp_t triggerTimeStamp,
      detinfo::DetectorClocksData const& clockData)
    {
      return TriggerToWaveformTimeStamp(
        triggerTimeStamp, clockData.OpticalClock().TickPeriod());
    }

    static TriggerTimeDiff_t DifferenceTicks(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t triggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return DifferenceTicks(
        WaveformToTriggerTicks(waveform, opticalClockTickPeriod),
        triggerTimeStamp);
    }

    static TriggerTimeDiff_t DifferenceTicks(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t triggerTimeStamp,
      detinfo::DetectorClocksData const& clockData)
    {
      return DifferenceTicks(
        waveform, triggerTimeStamp, clockData.OpticalClock().TickPeriod());
    }

    static TriggerTimeDiff_t DifferenceTicks(
      WaveformTimeStamp_t waveformTimeStamp,
      TriggerTimeStamp_t triggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return DifferenceTicks(
        WaveformToTriggerTicks(waveformTimeStamp, opticalClockTickPeriod),
        triggerTimeStamp);
    }

    static TriggerTimeDiff_t DifferenceTicks(
      TriggerTimeStamp_t waveformTrimmedTicks,
      TriggerTimeStamp_t triggerTimeStamp)
    {
      const auto waveformLow = waveformTrimmedTicks & TrimmedMask;
      const auto triggerLow = triggerTimeStamp & TrimmedMask;

      auto difference =
        static_cast<TriggerTimeDiff_t>(waveformLow) -
        static_cast<TriggerTimeDiff_t>(triggerLow);
      const auto halfRange =
        static_cast<TriggerTimeDiff_t>(TrimmedModulus >> 1U);

      if (difference >= halfRange) {
        difference -= static_cast<TriggerTimeDiff_t>(TrimmedModulus);
      }
      else if (difference < -halfRange) {
        difference += static_cast<TriggerTimeDiff_t>(TrimmedModulus);
      }

      return difference;
    }

    static double Difference(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t triggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return static_cast<double>(
               DifferenceTicks(waveform, triggerTimeStamp, opticalClockTickPeriod)) *
             opticalClockTickPeriod;
    }

    static double Difference(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t triggerTimeStamp,
      detinfo::DetectorClocksData const& clockData)
    {
      return Difference(
        waveform, triggerTimeStamp, clockData.OpticalClock().TickPeriod());
    }

    static double Difference(
      WaveformTimeStamp_t waveformTimeStamp,
      TriggerTimeStamp_t triggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return static_cast<double>(
               DifferenceTicks(waveformTimeStamp, triggerTimeStamp, opticalClockTickPeriod)) *
             opticalClockTickPeriod;
    }

    static TriggerTimeStamp_t WaveformToTriggerTicksNear(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t referenceTriggerTimeStamp,
      double opticalClockTickPeriod)
    {
      return WaveformToTriggerTicksNear(
        WaveformToTriggerTicks(waveform, opticalClockTickPeriod),
        referenceTriggerTimeStamp);
    }

    static TriggerTimeStamp_t WaveformToTriggerTicksNear(
      raw::OpDetWaveform const& waveform,
      TriggerTimeStamp_t referenceTriggerTimeStamp,
      detinfo::DetectorClocksData const& clockData)
    {
      return WaveformToTriggerTicksNear(
        waveform, referenceTriggerTimeStamp, clockData.OpticalClock().TickPeriod());
    }

    static TriggerTimeStamp_t WaveformToTriggerTicksNear(
      TriggerTimeStamp_t waveformTrimmedTicks,
      TriggerTimeStamp_t referenceTriggerTimeStamp)
    {
      const auto difference =
        DifferenceTicks(waveformTrimmedTicks, referenceTriggerTimeStamp);

      if (difference < 0) {
        const auto magnitude = static_cast<TriggerTimeStamp_t>(-difference);
        return (referenceTriggerTimeStamp < magnitude) ?
          TriggerTimeStamp_t{0} : referenceTriggerTimeStamp - magnitude;
      }

      const auto magnitude = static_cast<TriggerTimeStamp_t>(difference);
      const auto max = std::numeric_limits<TriggerTimeStamp_t>::max();
      return (max - referenceTriggerTimeStamp < magnitude) ?
        max : referenceTriggerTimeStamp + magnitude;
    }

  private:
    static TriggerTimeStamp_t RoundToTicks(
      WaveformTimeStamp_t waveformTimeStamp,
      double opticalClockTickPeriod)
    {
      if (opticalClockTickPeriod <= 0.0 || waveformTimeStamp <= 0.0) {
        return TriggerTimeStamp_t{0};
      }

      const auto ticks = std::llround(
        static_cast<long double>(waveformTimeStamp) /
        static_cast<long double>(opticalClockTickPeriod));

      return (ticks <= 0) ?
        TriggerTimeStamp_t{0} : static_cast<TriggerTimeStamp_t>(ticks);
    }
  };

} // namespace opdet

#endif // DUNEOPDET_OPTICALDETECTOR_ANA_OPDETWAVEFORMTIMESTAMPHELPER_H
