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


#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1


//This functor ShadingEffectOfCfTreeSegment<TS,BUD> evaluates shading caused
//by a conifer segment on this conifer segment (shaded_s)

template <class TS,class BUD>
TreeCompartment<TS,BUD>* ShadingEffectOfCfTreeSegment_1<TS,BUD>::operator()(TreeCompartment<TS,BUD>* tc)const {

    //  int beamShading(Point& p0, PositionVector& v,

    if (CfTreeSegment<TS,BUD>* ts = dynamic_cast<CfTreeSegment<TS,BUD>*>(tc)) {
        //Don't compare to yourself
        if (ts == shaded_s)
            return tc;

        //Now go on computing shading
        int i = 0, number_of_sectors = 0, result = NO_HIT;
        double distance = 0.0;
        vector<double> radiation_direction(3);

        Tree<TS,BUD>& tt = GetTree(*ts);

        FirmamentWithMask& firmament = GetFirmament(tt);

        number_of_sectors = firmament.numberOfRegions();

        //Foliage density: Foliage area divided by  volume. Perhaps a good idea to
        //implement it as GetValue?
        LGMdouble af = GetValue(*ts,LGAAf);
        LGMdouble fol_dens;
        //NOTE: foliage density for cylinder containing also wood
        if(af > R_EPSILON)
            //      fol_dens = af/(PI_VALUE*(pow(GetValue(*ts,LGARf),2.0)-pow(GetValue(*ts,LGAR),2.0))
            fol_dens = af/(PI_VALUE*(pow(GetValue(*ts,LGARf),2.0))
                           *GetValue(*ts,LGAL));
        else
            fol_dens = 0.0;


        for (i = 0; i < number_of_sectors; i++) {
            //If the sector is blocked by another shoot
            //do not make computations, check the next sector instead
            if (S[i] == HIT_THE_WOOD) {
                continue;
            }
            //The radiation and its direction of sector i. We need the direction
            firmament.diffuseRegionRadiationSum(i,radiation_direction);

            Point r_0 =  GetPoint(*shaded_s)+0.5*GetValue(*shaded_s,LGAL)*
                    (Point)GetDirection(*shaded_s);        //Midpoint of shaded seg

            //If foliage, wood radius = 0
            if(fol_dens == 0.0) {
                LGMdouble rw = GetValue(*ts, LGAR);
		//                LGMdouble rf = rw;
		//old one
/*                 result = CylinderBeamShading(r_0, */
/*                                              radiation_direction, */
/*                                              GetPoint(*ts), */
/*                                              GetDirection(*ts), */
/*                                              rf, */
/*                                              rw, */
/*                                              GetValue(*ts, LGAL), */
/*                                              distance); */
                result = CylinderBeamShading(r_0,
                                             radiation_direction,
                                             GetPoint(*ts),
                                             GetDirection(*ts),
                                             rw,
                                             GetValue(*ts, LGAL));

            }
            else {
		//Old one:
/*                 result = CylinderBeamShading(r_0, */
/*                                              radiation_direction, */
/*                                              GetPoint(*ts), */
/*                                              GetDirection(*ts), */
/*                                              GetValue(*ts, LGARf), */
/*                                              rw, */
/*                                              GetValue(*ts, LGAL), */
/*                                              distance); */
                result = CylinderBeamShading(r_0,
                                             radiation_direction,
                                             GetPoint(*ts),
                                             GetDirection(*ts),
                                             GetValue(*ts, LGARf),
                                             GetValue(*ts, LGAL), distance);

            }

            if (result == HIT_THE_WOOD){
                //mark the sector blocked
                S[i] = HIT_THE_WOOD;
            }
            else if (result == HIT_THE_FOLIAGE){
                //otherwise compute Vp (the shadiness):
                //1. compute the inclination of light beam and the segment
                //1a. angle between segment and light beam
                double a_dot_b = Dot(GetDirection(*ts),
                                     PositionVector(radiation_direction));
                //1b.  inclination: Perpendicular  (PI_DIV_2) to segment minus
                //angle between segment and light beam
                double inclination = PI_DIV_2 - acos(fabs(a_dot_b));
                //2.the light extinction coefficient according to inclination
                double extinction = K(inclination);
                //3.Vp = extinction*distance*foliage_density
                double Vp = extinction *distance*fol_dens;
                S[i] += Vp;
            }

        }
    }
    return tc;
}


//=======================================================================================================
//This version of radiation evaluates radiation conditions for subject tree by pairwise
// comparison


// 2) EvaluateRadiationForCfTreeSegment_2 evaluates shading by all other segments by paiwise comparison (segments
// in own crown & other trees). It uses ShadingEffectOfCfTreeSegment_1<TS,BUD> to evaluate shading.
// This functor evaluates shading by own crown and shading by other trees (stand) separately and updates
// Qin_stand in TreeSegment.


template <class TS, class BUD, class TREE>
TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegment_2<TS,BUD,TREE>::operator() (TreeCompartment<TS, BUD>* tc)const
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

        vector<double> v(number_of_sectors,0.0);
        ShadingEffectOfCfTreeSegment_1<TS,BUD> s_e(ts,K,v);


        //This  goes  through  the  tree  and computes  shading  based  on
        //1)distance  light beam traverses  in foliage,  2)foliage density
        //and 3) inclination light beam hits the segment.

        //just go through all trees: first the others (stand) and then the tree itself.
        //This is to evaluate the self-shading and stand-shading components.

        //The target tree is the first tree in the vector; do it last in order
        //to get the effect of surrounding stand

        //In the case of only_self it is the only tree in the vector
        if(only_self) {
            TREE* t = vtree[0];
            ForEach(*t,s_e);
        }
        else {
            if(vtree.size() > 1)
                for(unsigned int k = 1; k < vtree.size(); k++) {
                    TREE* t = vtree[k];
                    ForEach(*t,s_e);
                }
        }

        //Now the Qin after shading of others
        vector<double> qis(number_of_sectors,0.0);
        vector<double>& ss = s_e.getS();
        for (int i = 0; i < number_of_sectors; i++){
            if (ss[i] != HIT_THE_WOOD){
                MJ Io = firmament.diffuseRegionRadiationSum(i,radiation_direction);
                qis[i] = Io*exp(-ss[i]);
            }
        }
        LGMdouble Qin_stand = 0.0;
        Qin_stand = accumulate(qis.begin(),qis.end(),0.0);
        ts->setQinStand(Qin_stand);

        //Now the last tree i.e tree itself
        //    ForEach(*(vtree[0]), s_e);

        //implement  "Ip  =  Iope^(-Vp)",  s[i] =  radiation  coming  from
        //direction i after this
        vector<double>& s = s_e.getS();
        for (int i = 0; i < number_of_sectors; i++){
            if (s[i] == HIT_THE_WOOD){
                s[i] = 0.0;
            }
            else {
                MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction);
                s[i] = Iop*exp(-s[i]);
            }
        }
        //Total incoming radiation
        MJ Q_in = accumulate(s.begin(),s.end(),0.0);

        //s contains now incoming radiation from each sector. Evaluate how
        //much segment absorbs from incoming radation.
        LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk;
        Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0;
        Lk = GetValue(*ts, LGAL);   //length is > 0.0, otherwise we would not bee here
        Rfk = GetValue(*ts, LGARf);  //Radius to foliage limit
        Wfk = GetValue(*ts, LGAWf); //Foliage mass
        //   sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree
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

            //implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it
            //should if there is any foliage)
            Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack;
            s[i] *= Ask;
        }

        MJ Q_abs = accumulate(s.begin(),s.end(),0.0);
        SetValue(*ts, LGAQabs, Q_abs);
        SetValue(*ts, LGAQin, Q_in);
    }
    return tc;
}


//=================================================================================================
//
// VoxelSpace calculation

template <class TS, class BUD>
TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegment_3<TS,BUD>::operator() (TreeCompartment<TS, BUD>* tc)const
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

        vector<double> v(number_of_sectors,0.0);
        //    ShadingEffectOfCfTreeSegment_1<TS,BUD> s_e(ts,K,v);
        //This  goes  through  the  tree  and computes  shading  based  on
        //1)distance  light beam traverses  in foliage,  2)foliage density
        //and 3) inclination light beam hits the segment.
        //HUOMMMM   !!!!!!!!
        //ForEach(tt,s_e);

        //implement  "Ip  =  Iope^(-Vp)",  s[i] =  radiation  coming  from
        //direction i after this
        //    vector<double>& s = s_e.getS();
        vector<double> s(number_of_sectors, 0.0);
        vector<double> qis(number_of_sectors, 0.0);
        // std::cout << std::noboolalpha << calculateDirectionalStar << " == " << std::boolalpha << calculateDirectionalStar << std::endl;

        AccumulateOpticalDepth AOD(voxel_space->getXSideLength(), par_a, par_b, middle,K,
                                   dir_effect, wood, correct_star, constant_star, calculateDirectionalStar);

        for (int i = 0; i < number_of_sectors; i++){
            MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction);

            //first attenuation in the voxel space
            LGMdouble transmission_voxel = 1.0;
            vector<VoxelMovement> vm;
            PositionVector dir(radiation_direction);

            voxel_space->getRoute(vm, middle, dir, K, false,calculateDirectionalStar);  //this shoud return only the "box route"


            //with traveled lengths
            //calculate the extinction coeffient
            //Consider also the mean direction of shoots in box
            AOD.beam_dir = dir;

            LGMdouble optical_depth = accumulate(vm.begin(),vm.end(),0.0,AOD);
            //If tree itself is in calculation (dump_self = true) subtract the effect of own foliage
            //NOTE Assumes that all own foliage is in the first (= the one that contains the middle point
            //of the segment) voxel box

            if(dump_self) {
                LGMdouble k =  max(0.0,-0.014+1.056*vm[0].STAR_mean);
                optical_depth -= k * GetValue(*ts,LGAAf) * vm[0].l / voxel_space->getBoxVolume();

                if(optical_depth < 0.0)
                    optical_depth = 0.0;
            }
            if(optical_depth > R_EPSILON){
	      if(optical_depth < 20.0) {
                    transmission_voxel = exp(-optical_depth);
	      }
	      else {
                    transmission_voxel = 0.0;
	      }
	    }
            Iop *= transmission_voxel;

            //then attenuation in the BorderForest
            if(evaluate_border)
                Iop *= border_forest->getBorderForestExtinction(middle, dir,k_border_conifer);

            qis[i] = Iop;
            s[i] = Iop;


            /*       if (s[i] == HIT_THE_WOOD){ */
            /* 	s[i] = 0.0; */
            /*       } */
            /*       else */
            /* 	s[i] = Iop*exp(-s[i]); */

        } //End of no_sectors ...


        //Total incoming radiation and radiation after stand
        LGMdouble Qin_stand = accumulate(qis.begin(),qis.end(),0.0);
        ts->setQinStand(Qin_stand);

        MJ Q_in = accumulate(s.begin(),s.end(),0.0);

        //s contains now incoming radiation from each sector. Evaluate how
        //much segment absorbs from incoming radation.
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

            //implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it
            //should if there is any foliage)
            Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack;
            s[i] *= Ask;
        }
        MJ Q_abs = accumulate(s.begin(),s.end(),0.0);
        SetValue(*ts, LGAQabs, Q_abs);
        SetValue(*ts, LGAQin, Q_in);
    }
    return tc;
}  //end of EvaluateRadiationForCfTreeSegment()  { ...



//=======================================================================================================

//This functor is a modification of ShadingEffectOfCfTreeSegment_1<TS,BUD> and calculates transmission
//through a conifer segment to a point in space from different directions.
//The point and the directions are data members of this functor.

template <class TS,class BUD>
  TreeCompartment<TS,BUD>* ShadingEffectOfCfTreeSegmentToPoint<TS,BUD>::operator()(TreeCompartment<TS,BUD>* tc)const {

  if (CfTreeSegment<TS,BUD>* ts = dynamic_cast<CfTreeSegment<TS,BUD>*>(tc)) {

    //Foliage density: Foliage area divided by volume. Perhaps a good idea to
    //implement it as GetValue?
    LGMdouble af = GetValue(*ts,LGAAf);
    LGMdouble fol_dens;
    //NOTE: foliage density for cylinder containing also wood
    if(af > R_EPSILON) {
      //      fol_dens = af/(PI_VALUE*(pow(GetValue(*ts,LGARf),2.0)-pow(GetValue(*ts,LGAR),2.0))
      fol_dens = af/(PI_VALUE*(pow(GetValue(*ts,LGARf),2.0))*GetValue(*ts,LGAL));
    }
    else {
      fol_dens = 0.0;
    }

    for(int i = 0; i < number_of_directions; i++) {

      if(ellipsoid_calculation && (ellipsoid_hits[i] == 0)) {   //In case of ellipsoid calculation, don't
	                                             //calculate if has not intercepted an ellipsoid crown
	continue;
      }

      //If the sector is blocked by another shoot
      //do not make computations, check the next sector instead
      if (S[i] == HIT_THE_WOOD) {
	continue;
      }
      //The radiation and its direction of sector i. We need the direction
      //            firmament.diffuseRegionRadiationSum(i,radiation_direction);

      vector<double> radiation_direction(3);
      radiation_direction = directions[i];

      int result = NO_HIT;
      double distance = 0.0;


      Point r_0 = p0;

      //If foliage, wood radius = 0
      if(fol_dens == 0.0) {
	LGMdouble rw = GetValue(*ts, LGAR);
	//old version
/* 	LGMdouble rf = rw; */
/* 	result = CylinderBeamShading(r_0, */
/* 				     radiation_direction, */
/* 				     GetPoint(*ts), */
/* 				     GetDirection(*ts), */
/* 				     rf, */
/* 				     rw, */
/* 				     GetValue(*ts, LGAL), */
/* 				     distance) */
	result = CylinderBeamShading(r_0,
				     radiation_direction,
				     GetPoint(*ts),
				     GetDirection(*ts),
				     rw,
				     GetValue(*ts, LGAL));
      }
      else {
	//old version
/* 	LGMdouble rw = 0.0; */
/* 	result = CylinderBeamShading(r_0, */
/* 				     radiation_direction, */
/* 				     GetPoint(*ts), */
/* 				     GetDirection(*ts), */
/* 				     GetValue(*ts, LGARf), */
/* 				     rw, */
/* 				     GetValue(*ts, LGAL), */
/* 				     distance); */
	result = CylinderBeamShading(r_0,
				     radiation_direction,
				     GetPoint(*ts),
				     GetDirection(*ts),
				     GetValue(*ts, LGARf),
				     GetValue(*ts, LGAL),
				     distance);
      }
      if (result == HIT_THE_WOOD){
	//mark the sector blocked
	S[i] = HIT_THE_WOOD;
      }
      else if (result == HIT_THE_FOLIAGE){
	//1a. angle between segment and light beam
	double a_dot_b = Dot(GetDirection(*ts),
			     PositionVector(radiation_direction));
	//1b.  inclination: Perpendicular  (PI_DIV_2) to segment minus
	//angle between segment and light beam
	double inclination = PI_DIV_2 - acos(fabs(a_dot_b));
	//2.the light extinction coefficient according to inclination
	double extinction = K(inclination);
	//3.Vp = extinction*distance*foliage_density
	double Vp = extinction *distance*fol_dens;
	S[i] += Vp;
      }

    }  //for(int i = 0; i < number_of_directions; ...
  }
  return tc;
}

#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD

#endif
