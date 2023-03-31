#ifndef CALCULATE_LIGHTI_H
#define CALCULATE_LIGHTI_H
/* #include <Lignum.h> */
/* #include <ScotsPine.h> */
/* #include <CopyDumpCfTree.h> */
/* #include <VoxelSpace.h> */
/* #include <SomeFunctors.h> */
/* #include <vector> */
/* #include <utility> */
/* #include <BorderForest.h> */
/* using namespace std; */

/* using namespace Lignum; */
/* using namespace sky; */


//=================================================================================================
//
// VoxelSpace calculation

namespace LignumForest{
template <class TS, class BUD>
TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegmentInVoxelSpace<TS,BUD>::operator()
  (TreeCompartment<TS, BUD>* tc)const
{
    if (TS* ts = dynamic_cast<TS*>(tc)){
        SetValue(*ts, LGAQin, 0.0);
        SetValue(*ts, LGAQabs, 0.0);
        //Radiation  conditions are not  evaluated if  the segment  has no
        //foliage (in practice  there would be division by  0 in computing
        //absorbed radiation)
        if (GetValue(*ts, LGAWf) < R_EPSILON){
            return tc;
        }

        Tree<TS,BUD>& tt = GetTree(*ts);
        FirmamentWithMask& firmament = GetFirmament(tt);
        int number_of_sectors = firmament.numberOfRegions();
        double a_dot_b = 0.0;
        vector<double> radiation_direction(3);
        Point middle = GetMidPoint(*ts);

        vector<double> v_pairwise_self(number_of_sectors,0.0);
	if(pairwise_self) {
	  ShadingEffectOfCfTreeSegment<TS,BUD> s_e(ts,K,v_pairwise_self);
	  //This  goes  through  the  tree  and computes  shading  based  on
	  //1)distance  light beam traverses  in foliage,  2)foliage density
	  //and 3) inclination light beam hits the segment.
	  //vector v_pairwise_self now contains optical depths (== -ln(transmittance))
	  //in the directions of sky sectors
	  ForEach(tt,s_e);
	}
        vector<double> s(number_of_sectors, 0.0);
 
        AccumulateOpticalDepth AOD(voxel_space->getXSideLength(), wood);
        for (int i = 0; i < number_of_sectors; i++){
            MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction);

            //first attenuation in the voxel space
            vector<VoxelMovement> vm;
            PositionVector dir(radiation_direction);

            voxel_space->getRoute(vm, middle, dir, K, false,false);
            LGMdouble optical_depth = accumulate(vm.begin(),vm.end(),0.0,AOD);

	    optical_depth += v_pairwise_self[i]; //Possible contribution of own crown by ray casting

            //Subtract the effect of own foliage in the case of not pairwise_self
            //NOTE Assumes that all own foliage is in the first voxel box
	    //(= the one that contains the middle point of the segment)
	    if(!pairwise_self) {
	      LGMdouble k =  vm[0].STAR_mean;
	      optical_depth -= k * GetValue(*ts,LGAAf) * vm[0].l / voxel_space->getBoxVolume();
	      if(optical_depth < 0.0) {
		optical_depth = 0.0;
	      }
	    }
	     
	    LGMdouble transmittance = 1.0;
            if(optical_depth > R_EPSILON){
	      if(optical_depth < 20.0) {
                    transmittance = exp(-optical_depth);
	      }
	      else {
                    transmittance = 0.0;
	      }
	    }
            Iop *= transmittance;

            //then attenuation in the BorderForest
            if(evaluate_border) {
                Iop *= border_forest->getBorderForestExtinction(middle, dir,k_border_conifer);
	    }
            s[i] = Iop;
        } //End of no_sectors ...

        MJ Q_in = accumulate(s.begin(),s.end(),0.0);

        //s contains now incoming radiation from each sector. Evaluate how
        //much segment absorbs from incoming radation.
	//The interception is realized acconding to P. Oker-Blom and H. Smolander.
	//Forest Science, 34(4):894–906, 1988

        LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk;
        Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0;
        Lk = GetValue(*ts, LGAL);   //length is > 0.0, otherwise we would not bee here
        Rfk = GetValue(*ts, LGARf);  //Radius to foliage limit
        Wfk = GetValue(*ts, LGAWf); //Foliage mass
        //sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree
        sfk  = GetValue(*ts, LGAsf); //Foliage m2/kg from segment!!!

        for (int i = 0; i < number_of_sectors; i++){
            firmament.diffuseRegionRadiationSum(i,radiation_direction);
            a_dot_b = Dot(GetDirection(*ts), PositionVector(radiation_direction));
            inclination = PI_DIV_2 - acos(fabs(a_dot_b));

            Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination);
            extinction = (double)K(inclination);

            if (Ack == 0.0){
                cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)"
                     << endl;
            }

            Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack;
            s[i] *= Ask;
        }
        MJ Q_abs = accumulate(s.begin(),s.end(),0.0);
        SetValue(*ts, LGAQabs, Q_abs);
        SetValue(*ts, LGAQin, Q_in);
    }   
    return tc;
}  //end of EvaluateRadiationForCfTreeSegmentInVoxelSpace  { ...
}
#endif
