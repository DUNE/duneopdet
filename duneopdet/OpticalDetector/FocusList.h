#ifndef duneopdet_OpticalDetector_FocusList_h
#define duneopdet_OpticalDetector_FocusList_h

#include <algorithm>
#include <utility>
#include <vector>

namespace opdet {

  // Sample ranges of a waveform that need to be read out. Call Finalize()
  // after the last AddRange() and before reading ranges.
  class FocusList {
  public:
    FocusList(int nSamples, int padding) : fNSamples(nSamples), fPadding(padding) {}

    // Pad the range, clamp it to the waveform, and store it
    void AddRange(int from, int to)
    {
      from = std::max(from - fPadding, 0);
      to = std::min(to + fPadding, fNSamples - 1);
      ranges.emplace_back(from, to);
    }

    // Sort the ranges and merge overlapping ones
    void Finalize()
    {
      if (ranges.empty()) return;
      std::sort(ranges.begin(), ranges.end());
      std::vector<std::pair<int, int>> merged{ranges.front()};
      for (size_t i = 1; i < ranges.size(); ++i) {
        if (ranges[i].first <= merged.back().second)
          merged.back().second = std::max(merged.back().second, ranges[i].second);
        else
          merged.push_back(ranges[i]);
      }
      ranges = std::move(merged);
    }

    // Replace the ranges with the whole waveform
    void Reset() { ranges = {{0, fNSamples - 1}}; }

    std::vector<std::pair<int, int>> ranges;

  private:
    int fNSamples;
    int fPadding;
  };

}

#endif
