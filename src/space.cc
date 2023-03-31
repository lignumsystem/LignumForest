//About space colonization etc.

#include <Space.h>
namespace LignumForest{
bool SphericalSector::in_SphericalSector(const Point& p) {
  Point ap = p - apex;
  double apl = ap.getLength();
  if(apl > height) {
    return false;
  }

  double proj = Dot(ap, direction);
  if(proj <= 0.0) {
    return false;
  }

  if(sqrt(apl*apl - proj*proj) > tan(half_angle)*proj) {
    return false;
  }

  return true;
}
}
