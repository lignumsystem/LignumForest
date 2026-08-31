/// \file CalculateLight.h
/// \brief Radiation calculations
///
/// Two functors for Scots pine
/// \arg \c EvaluateRadiationForCfTreeSegmentInVoxelSpace
/// \arg \c AccumulateOpticalDepth
#ifndef CALCULATE_LIGHT_H
#define CALCULATE_LIGHT_H
#include <Lignum.h>
#include <ScotsPine.h>
#include <CopyDumpCfTree.h>
#include <VoxelSpace.h>
#include <SomeFunctors.h>
#include <vector>
#include <utility>
#include <BorderForest.h>
#include<limits>
using namespace std;

using namespace Lignum;
using namespace sky;


#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1

//===============================================================================
// VoxelSpace calculation
namespace LignumForest{
  
  template <class TS, class BUD>
  ///\brief Calculate tree segment radiation climate
  ///
  ///\anchor radiative_properties
  ///\par Radiative properties: transmittance and optical depth
  ///
  ///Transmittance (\f$ T \f$) is the fraction of light that passes through a material,
  ///while optical depth (\f$ \tau \f$) is a measure of how much light is blocked
  ///or scattered by that material.
  ///
  ///Transmittance has the range \f$[0,1]\f$, from completely opaque to completely transparent on a linear scale.
  ///Optical depth has the range \f$[0,\inf]\f$, from no absorption or scattering to infinity on an additive, logarithmic scale.
  ///
  ///More precisly, transmittance is the ratio of transmitted light intensity (\f$I\f$) to the initial incident
  ///light intensity (\f$ I_0 \f$), expressed as a fraction (\f$ T = I / I_0 \f$) or percentage.
  ///The Beer-Lambert law links transmittance and optical depth by \f$ T = \exp(-\tau) \f$, 
  ///meaning optical depth is the negative natural logarithm of transmittance (\f$ \tau = -\ln(T) \f$).
  ///
  ///Due to their mathematical relationship optical depths are additive, while transmittances are multiplicative.
  ///For example, consider two consecutive voxels (\f$\mathrm{v1,v2}\f$) that each pass 50% (0.5) of incoming light. The first voxel
  ///passes 0.5 of the initial light and the second one passes 0.5 of the \e remaining light of the first voxel.
  ///The total transmittance is \f$ 0.5 \times 0.5 = 0.25\f$. Applying the Beer–Lambert law clarifies this behavior:
  ///\f[
  /// T_{\mathrm{tot}} = \exp(-\tau_{\mathrm{v1}}) \times \exp(-\tau_{\mathrm{v2}}) = \exp(-(\tau_{\mathrm{v1}}+\tau_{\mathrm{v2}}))
  ///\f]
  ///The additive nature of optical depth is clearly seen in the exponent.
  ///
  ///Optical depth measures the cumulative attenuation along the light path. Consequently, its value is always non-negative.
  ///In the example both voxels exhibit optical depth of  \f$ \tau = -\ln(0.5) \approx 0.69315 \f$ units of attenuation,
  ///and the total optical depth is \f$ \approx 1.3863 \f$.
  class EvaluateRadiationForCfTreeSegmentInVoxelSpace {
  public:
    EvaluateRadiationForCfTreeSegmentInVoxelSpace(const ParametricCurve& k) : K(k),
									      evaluate_border(false) {}
    EvaluateRadiationForCfTreeSegmentInVoxelSpace(const ParametricCurve& k,
						  VoxelSpace* vs,BorderForest* bf, bool border,
						  LGMdouble kbc, bool wd, bool ps):
      K(k),voxel_space(vs), border_forest(bf), evaluate_border(border),
      k_border_conifer(kbc), wood(wd), pairwise_self(ps) {}
    ///\brief Calculate the incoming and absorbed radiation
    ///
    ///To calculate a tree segment's incoming and absorbed radiation, these sequential operations are applied:
    /// -# Pairwise tree segment comparison for self-shading:
    ///  - Contorolled by boolean \p pairwise_self
    ///    \snippet{lineno} CalculateLightI.h PairwiseComparison
    ///  - Self-shading optical depths for each sky sector are stored in the <b>\p v_pairwise_self </b> vector.
    /// -# Calculate the radiative properties along the voxel path between the target tree segment and each sky sector:
    ///   -# The incident radiant intensity, the original intensity of the radiation \p Iop from a sky sector:
    ///     \snippet{lineno} CalculateLightI.h InitialI0
    ///   -# Optical depths based on the initial incident light intensity from one sky sector:
    ///     \snippet{lineno} CalculateLightI.h OpticalDepth
    ///     - The default values in <b>\p v_pairwise_self </b> are 0.
    ///     - AccumulateOpticalDepth  \p AOD summates the independent optical depth values of individual voxels.
    ///     - It is correct to subtract self-shading by segment's own foliage to not over estimate shading by \f$\overline{\mathrm{STAR}}\f$,
    ///       but make sure the logic matches calculations in \p AOD.
    ///     .
    ///    \note  AccumulateOpticalDepth includes the wood part of the segment.
    ///    \sa AccumulateOpticalDepth
    ///    \todo  \htmlonly <span style="color:red;">&#9679;</span>\endhtmlonly Check the logic how the self-shading is treated.
    ///
    ///   -# Transmittance from one sky sector:
    ///      \snippet{lineno} CalculateLightI.h TransmittanceForest
    ///      \todo \htmlonly <span style="color:red;">&#9679;</span>\endhtmlonly Ensure the logic correctly flags or handles optical depths greater than 20.
    ///
    ///   -# Transmittance from the border forest:
    ///      \snippet{lineno} CalculateLightI.h TransmittanceBorderForest
    ///   \sa BorderForest::getBorderForestExtinction()
    ///   -# Total transmittance from one sector:
    ///      \snippet{lineno} CalculateLightI.h TransmittanceFromOneSector
    ///      - Transmittances for each sky sector are stored in the <b>\p s </b> vector.
    /// -# Calculate the total incoming radiation from all sectors:
    ///    \snippet{lineno} CalculateLightI.h QinTotal
    /// -# Calculate the segment's absorbed radiation for each sky sector:
    ///    - The absorption is realized according to Oker-Blom and Smolander. 1988. Forest Science,34(4):894–906.
    ///    - Calculate the radiation absorption coeffient, denoted \p Ask.
    ///    - Segment angle of incidence and consequent shadow area are used as the measure of absorption efficiency
    ///      \snippet{lineno} CalculateLightI.h QabsSectors
    ///    - The absorbed radiation is in the <b>\p s</b> vector.
    /// -# Total absorbed radiation from all sky sectors:
    ///      \snippet{lineno}  CalculateLightI.h QabsTotal
    /// -# Save the calculated values for incoming and absorbed radiation in the tree segment:
    ///    \snippet{lineno} CalculateLightI.h SaveQinQabs
    ///.
    ///\par Mathematical details for the radiation absorption coeffient
    ///
    ///According to the Beer-Lambert law for radiation transmittance (\f$ T = \exp(-\tau) \f$)
    ///the cumulative optical depth or opacity (\f$ \tau \f$) is a function of 
    ///extinction coefficient (\f$ k \f$) and the localized leaf area index or area density:
    ///\f[
    /// \tau = k\frac{A_f}{A_{\mathrm{ref}}}
    ///\f]
    ///where \f$ A_f \f$ is the foliage area, and \f$ A_{\mathrm{ref}} \f$ is the silhoutte - that is, cross-sectional projection - of the object.
    ///Substituting \f$ \tau \f$ into the transmittance formula yields:
    ///\f[
    /// T = \exp\left(-k\frac{A_f}{A_{\mathrm{ref}}}\right)
    ///\f]
    ///which expresses the probability light escaping absorption entirely.
    ///
    ///Assuming no scattering, the light can either be transmitted or absorbed. The fraction absorbed
    ///(\f$ \alpha \f$) by a tree segment is \f$ 1.0-T \f$:
    ///\f[
    /// \alpha = 1.0 - \exp\left(-k\frac{A_f}{A_{\mathrm{ref}}}\right)
    ///\f]
    ///To convert the relative absorption fraction into spatial metrics, that fraction
    ///must be scaled by the object's silhouette footprint or projection area (\f$ A_{\mathrm{geom}} \f$).
    ///Consequently, the total radiation absorption area - or the radiation absorption
    ///coefficient (\f$A_{\mathrm{abs}}\f$) - is expressed as:
    ///\f[
    /// A_{\mathrm{abs}} = \alpha A_{\mathrm{geom}}
    ///\f]
    ///In many applications, such as with tree segment cylinders, \f$ A_{\mathrm{ref}} = A_{\mathrm{geom}} \f$,
    ///yet they serve different roles in the calculation. \f$ A_{\mathrm{ref}} \f$ acts as a normalizing
    ///factor that scales the foliage area to a local leaf area index, ensuring the exponent is unitless.
    ///Conversely, \f$ A_{\mathrm{geom}} \f$ is an absolute multiplier that scales the the relative absorption
    ///fraction back into physical geometric area.
    ///\sa \ref radiative_properties "Radiative properties"
    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;

  private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_border;
    LGMdouble k_border_conifer;
    bool wood;         //If the needleless woody parts are included into calculation
    bool pairwise_self; //If self-shading is evaluated by ray casting
  };

  ///\brief Calculate the optical depth of a single voxel
  ///
  ///The total optical depth along a ray's path is then computed
  ///by summing these individual voxel depths using std::accumulate().
  class AccumulateOpticalDepth{
  public:
    ///Constructor
    ///\param side Voxel side length
    ///\param wd Calculate the wood optical depth
    AccumulateOpticalDepth(LGMdouble side, bool wd) :
      box_side_length(side), wood(wd)
    {
      box_volume = pow(box_side_length,3.0);
    }
    ///\brief Calculate the optical depth of a single voxel.
    ///
    ///Calculate the optical depth of a single voxel, accounting for both foliage and woody material.
    ///Evaluation of the wood area is conditional on the \p wood flag.
    ///\param o_d Optical depth
    ///\param vm Voxel data
    ///\retval o_d The sum of foliage and wood optical depths 
    ///
    ///\sa AccumulateOpticalDepth()
    ///\sa AccumulateOpticalDepth::wood
    ///\sa voxelspace::VoxelMovement
    ///\sa voxelspace::VoxelSpace::getRoute()
    ///\sa \ref radiative_properties "Radiative properties"
    ///\par Formulas for optical depth
    ///
    ///Calculate the single-voxel foliage optical depth (\f$ \tau_f \f$) based on voxel's \f$ \overline{\mathrm{STAR}}_\mathrm{voxel} \f$:
    ///\f[
    /// \tau_f = \frac{\overline{\mathrm{STAR}}_\mathrm{voxel}A_fL}{V_{\mathrm{voxel}}}
    ///\f]
    ///where \f$ A_f \f$ is the total foliage area within the voxel, \f$ L \f$ is the ray's path length through it,
    ///and \f$ V_{\mathrm{voxel}} \f$ is the volume of the voxel.
    ///
    ///Calculate the single-voxel wood optical depth (\f$ \tau_w \f$) based on the mean cylinder projection area (\f$\overline{A}_p\f$):
    ///\f[
    /// \tau_w = \frac{\overline{A}_pA_wL}{V_{\mathrm{voxel}}}
    ///\f]
    ///where \f$ A_w \f$ is the voxel wood area.
    ///
    ///The constant \f$ \overline{A}_p = \frac{1}{4} \f$ is based on Cauchy's surface area theorem.
    ///
    ///\par Single-voxel optical depth implementation
    ///
    ///The variable \p o_d represents the total optical depth, calculated as the sum of the optical depths along the entire ray path. 
    ///\snippet{lineno} CalculateLight.h VoxelOpticalDepth
    ///
    ///\sa voxelspace::VoxelSpace::getRoute()
    ///\sa voxelspace::VoxelBox::updateValues()
    ///\par Mathematical details for Cauchy's surface area theorem
    ///
    ///Define  \f$ \overline{A}_p \f$ as the weighted integral of the projected area divided
    ///by the total integrated weight over all possible 3D orientations. Due to the symmetry of the right cylinder the integration is performed
    ///from \f$ 0 \f$ to \f$ \frac{\pi}{2} \f$:
    ///\f[
    /// \overline{A}_p = \frac{\int_0^{\pi/2}A_p(\Theta)\sin\Theta d\Theta}{\int_0^{\pi/2}\sin\Theta d\Theta}
    ///\f]
    ///where  \f$ \sin\Theta \f$ is the standard geometric weighting factor; cylinder orientations are equally probable but
    ///sweeping the sphere with \f$ \Theta \f$ does not cover equal areas. Denominator is 1 and substituting the angled view formula:
    ///\f[
    /// \overline{A}_p = \int_0^{\pi/2}(\pi r^2\ cos\Theta + 2rh\sin\Theta)\sin\Theta d\Theta
    ///\f]
    ///
    ///The integration (using the trigonometric identities \f$\cos\theta \sin\theta = \frac{1}{2}\sin(2\theta)\f$ and
    ///\f$\sin^2\theta = \frac{1 - \cos(2\theta)}{2}\f$) followed by a simple algebraic manipulation reveals the total surface area formula:
    ///\f[
    ///  \overline{A}_p = \frac{1}{2}\pi r^2 + \frac{1}{2}\pi rh = \frac{2}{4}\pi r^2 + \frac{2}{4}\pi rh = \frac{1}{4}(2\pi r^2 + 2\pi rh) = \frac{1}{4}A
    ///\f]
    ///This completes Cauchy's surface area theorem for the flat faced right cylinder.
    double operator()(double o_d, VoxelMovement& vm){
      if(vm.af > R_EPSILON || (wood && vm.wood_area > R_EPSILON)) {
	// [VoxelOpticalDepth]
	//Foliage
	LGMdouble k = vm.STAR_mean;
	o_d += k * vm.af * vm.l / box_volume;
	if(wood) {
	  //Wood
	  //Cauchy's theorem: Mean projection area of surface of a circular cylinder
	  //(excluding end disks) is 1/4 of its area
	  o_d += 0.25 * vm.wood_area * vm.l / box_volume;
	}
	// [VoxelOpticalDepth]
      }
      return o_d;
    }
  private:
    LGMdouble box_side_length;
    LGMdouble box_volume;
    bool wood;                       //  If woody parts are considered
  };

#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD
}//End namespace LignumForest
#include <CalculateLightI.h>

#endif
