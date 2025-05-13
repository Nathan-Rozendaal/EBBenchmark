#include <iostream>
#include <chrono>

#include <EveryBeam/facade/load.h>
#include <EveryBeam/options.h>
#include <EveryBeam/load.h>
#include <EveryBeam/facade/response.h>
#include <EveryBeam/coords/itrfconverter.h>
#include <EveryBeam/telescope/telescope.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#include <casacore/ms/MeasurementSets.h>
#include <xtensor.hpp>
#include <aocommon/matrix2x2.h>
#include <aocommon/banddata.h>

void ResponseNewTelescope(const std::vector<int> &stations,
                          const std::unique_ptr<everybeam::facade::Telescope> &newtelescope,
                          const std::vector<everybeam::vector3r_t> &directionsXtime,
                          const std::vector<double> &timeXdirections,
                          const std::vector<double> &freqs);
void ResponseOldTelescope(std::vector<int> &stations,
                          const std::unique_ptr<everybeam::telescope::Telescope> &telescope,
                          std::vector<everybeam::vector3r_t> &directions,
                          const std::vector<double> &times,
                          const std::vector<double> &freqs);
/// Normalize a 3‐vector
static everybeam::vector3r_t normalize(everybeam::vector3r_t v) {
  double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return { v[0]/n, v[1]/n, v[2]/n };
}

/// Generate directions
std::vector<everybeam::vector3r_t> sampleNear(
    everybeam::vector3r_t center,
    size_t N,
    double max_angle_rad = 5.0 * M_PI / 180.0
) {
  std::mt19937_64 rng{std::random_device{}()};
  std::uniform_real_distribution<double> unif01(0.0, 1.0);
  // rough offset magnitude scale = sin(max_angle)
  double offset_magnitude = std::sin(max_angle_rad);

  std::vector<everybeam::vector3r_t> out;
  out.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    // small random offset vector in [-δ,δ]^3
    everybeam::vector3r_t offset = {
        (unif01(rng)*2.0 - 1.0) * offset_magnitude,
        (unif01(rng)*2.0 - 1.0) * offset_magnitude,
        (unif01(rng)*2.0 - 1.0) * offset_magnitude
    };
    everybeam::vector3r_t modified = normalize({ center[0] + offset[0],
                                       center[1] + offset[1],
                                       center[2] + offset[2] });
    out.push_back(modified);
  }
  return out;
}

int main() {
  constexpr int num_stations = 10; //max 31
  constexpr int num_freqs = 5; // max 900
  constexpr int num_directions = 150;
  constexpr int num_times = 3; // max 14
  using namespace std::chrono;

  //Setup
  std::vector<int> stations;
  for (int i = 0; i < num_stations; ++i)
  {
    stations.push_back(i);
  }
  const casacore::MeasurementSet ms("../measurementsets/TC00_first30.MS");

  everybeam::Options options;
  options.use_channel_frequency = true;
  std::unique_ptr<everybeam::telescope::Telescope> telescope = everybeam::Load(ms, options);
  std::unique_ptr<everybeam::facade::Telescope> newtelescope = everybeam::facade::Load(ms);
  casacore::ROScalarColumn<double> time_column(ms, "TIME");

  std::set<double> timeset;
  for (casacore::uInt i = 0; i < time_column.nrow(); ++i) {
    timeset.insert(time_column(i));
  }
  std::vector<double> time(timeset.begin(), timeset.end());

  everybeam::coords::ItrfConverter itrf_converter(time.front());
  std::pair<double,double> radecpointing = newtelescope->antennas[0].pointing_direction;
  auto itrfpointing = itrf_converter.RaDecToItrf(radecpointing.first, radecpointing.second);

  auto directions = sampleNear(itrfpointing, num_directions);
  std::vector<everybeam::vector3r_t> directionsXtime;
  directionsXtime.reserve(directions.size() * time.size());
  for (int i = 0; i < time.size(); ++i) {
    directionsXtime.insert(directionsXtime.end(), directions.begin(), directions.end());
  }

  std::vector<double> timeXdirections;
  timeXdirections.reserve(directionsXtime.size());
  for (double i : time) {
    for (int j = 0; j < directions.size(); ++j) {
      timeXdirections.push_back(i);
    }
  }
  assert(timeXdirections.size() == directionsXtime.size());

  aocommon::BandData band(ms.spectralWindow());
  std::vector<double> freqs;
  for (int i = 0; i < num_freqs; ++i) {
    freqs.push_back(band.ChannelFrequency(i));
  }


  //Benchmark
  auto startA = high_resolution_clock::now();
  ResponseOldTelescope(stations, telescope, directions, time, freqs);
  auto endA = high_resolution_clock::now();
  auto durationA = duration_cast<milliseconds >(endA - startA).count();
  std::cout << "Old response took " << durationA << " milliseconds\n";

  auto startB = high_resolution_clock::now();
  ResponseNewTelescope(stations, newtelescope, directionsXtime, timeXdirections, freqs);
  auto endB = high_resolution_clock::now();
  auto durationB = duration_cast<milliseconds >(endB - startB).count();
  std::cout << "New response took " << durationB << " milliseconds\n";

  return 0;
}
void ResponseOldTelescope(std::vector<int> &stations,
                          const std::unique_ptr<everybeam::telescope::Telescope> &telescope,
                          std::vector<everybeam::vector3r_t> &directions,
                          const std::vector<double> &times,
                          const std::vector<double> &freqs) {
  int count = 0;
  for (double time: times) {
    std::unique_ptr<everybeam::pointresponse::PointResponse> response =
        telescope->GetPointResponse(time);
    for (int station : stations) {
      for (double freq : freqs){
        for (const auto & direction : directions) {
          count ++;
          response->Response( everybeam::BeamMode::kFull,station,freq, direction);
        }
      }
    }
  }
  std::cout << "Old response count: " << count << "\n";
}
void ResponseNewTelescope(const std::vector<int> &stations,
                          const std::unique_ptr<everybeam::facade::Telescope> &newtelescope,
                          const std::vector<everybeam::vector3r_t> &directionsXtime,
                          const std::vector<double> &timeXdirections,
                          const std::vector<double> &freqs) {
  auto responses  = Response(*newtelescope, stations, freqs, directionsXtime,
                              timeXdirections);
  std::cout << "New response count: " << responses.size() << "\n";
}