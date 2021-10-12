#include <Well.hpp>


Well::Well(const WellConfig & config)
    :
    name(config.name),
    radius(config.radius),
    perforated(config.perforated),
    _force_frac_connect(config.force_connect_fractures)
{
  assert(!config.coordinates.empty());
  assert(config.radius > 0 );

  if ( config.coordinates.size() == 1 ) {
    // simple well
    coordinate = config.coordinates[0];
    reference_depth = -config.coordinates[0][2];
  }
  else
  {
    segments.reserve( perforated.size() );
    segments.emplace_back();
    bool segment_open = false;
    std::size_t i = 0;
    while (i < config.coordinates.size())
    {
      auto & segment = segments.back();
      if (!segment_open)
      {
        segment.first = config.coordinates[i];
        segment_open = true;
        i++;
      }
      else
      {
        segment.second = config.coordinates[i];
        segment_open = false;

        if (i < config.coordinates.size() - 1)
          segments.emplace_back();
        else
          break;
      }
    }
  }
}
