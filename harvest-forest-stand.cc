/// \file harvest-forest-stand.cc
/// \brief HarvestForestStand and RemoveDeadTrees implementation
#include <sstream>
#include <ScotsPine.h>
#include <HarvestForestStand.h>


namespace LignumForest{ 
  ///Stand harvest when the trees are individuals: remove a location with
  ///probability 'r' to  get the wanted stand density,  remove and delete
  ////the tree, remove and delete its L-system
  int HarvestForestStand(vector<pair<double,double> >& v,
			 vector<Tree<ScotsPineSegment,ScotsPineBud>*>& tv,
			 vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
			 PBNAME,PineBudData>*>& lv, double r)
  {
    int r_seed = 1;
    vector<pair<double,double> >::iterator current = v.begin();
    vector<Tree<ScotsPineSegment,ScotsPineBud>*>::iterator ct = tv.begin();
    vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>*>::iterator lt = lv.begin();

    //scan the location vector and remove  position if p <= r. Note that
    //the  tree to  be considered  is the  last in  the tree  vector and
    //v.size == tv.size-1, so it we won't touch the tree we study.
    while (current != v.end()){
      double p = ran3(&r_seed);
      if (p <= r){
	v.erase(current);
	//delete the tree and its vector element
	delete *ct;
	tv.erase(ct);
	//delete the L-system of the tree and its vector element
	Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>* lsys = *lt;
	delete lsys;
	lv.erase(lt);
      }
      else{
	advance(current,1);
	advance(ct,1);
	advance(lt,1);
      }
    }      
    //the final density
    return v.size();
  }


  ///Foliage dump  and carbon  allocation may fail  (out of  voxel space,
  ///bracket failure, allocation failure or nothing to allocate). In that
  ///case the  tree is  deleted and marked  NULL. Remove such  trees from
  ///simulation
  void RemoveDeadTrees(vector<pair<double,double> >& v,
		       vector<Tree<ScotsPineSegment,ScotsPineBud>*>& pinev,
		       vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
		       PBNAME,PineBudData>*>& plv)
  {
    vector<pair<double,double> >::iterator v_it = v.begin();
    vector<Tree<ScotsPineSegment,ScotsPineBud>*>::iterator pinev_it = pinev.begin();
    vector<Pine::LSystem<ScotsPineSegment,ScotsPineBud,
			 PBNAME,PineBudData>*>::iterator plv_it = plv.begin();
    while(pinev_it != pinev.end()){
      //if the tree is deleted
      if (GetValue(**pinev_it,LGAage) < 0){
	ostringstream error;
	Point p = GetPoint(**pinev_it);
	error << "Removing Tree at Point " << p.getX() << " " << p.getY() << " " << p.getZ() <<endl;
	LGMError(error.str());
	//delete the tree
	delete *pinev_it;
	//remove the tree
	pinev.erase(pinev_it);
	//delete the L-system
	Pine::LSystem<ScotsPineSegment,ScotsPineBud,PBNAME,PineBudData>* lsys = *plv_it;      
	delete lsys;
	//remove the L-system
	plv.erase(plv_it);
	//remove the position
	v.erase(v_it);
      }
      //After deletion iterator points  to the element immediately after
      //the deleted  one, possibly end(). Otherwise  move explicitely to
      //the next element.
      else{
	advance(pinev_it,1);
	advance(plv_it,1);
	advance(v_it,1);
      }
    }
    if (pinev.size() == 0){
      LGMError("No trees left, expect program exit");
    }
  }
}//End namespace LignmForest    
