#ifndef CREATE_VOXELSPACE_DATA_H
#define CREATE_VOXELSPACE_DATA_H
#include <TMatrix3D.h>
#include <VoxelSpace.h>
///\file CreateVoxelSpaceData.h
///\brief Collect VoxelSpace size data during the simulation.

namespace LignumForest{
  ///\brief Collect VoxelSpace lower left and upper right
  ///corner point data to a 2D matrix.
  ///
  ///The 2D matrix can be saved to a HDF5 file.
  ///\sa LignumForest::LGMHDF5File
  class CreateVoxelSpaceData{
  public:
    ///Create empty VoxelSpaceData. The number of data items
    ///is not predefined but the content of a data row, the number
    ///columns, can be decided later.
    ///\post The `current_row` = -1.
    ///\note Matrix indexing starts from 0 and `current_row` is updated
    ///after each row insertion
    ///\sa CreateVoxelSpaceData::insertData
  CreateVoxelSpaceData():
    current_row(-1){}
    ///Create predefined data matrix for collection years,
    ///lower left and upper right hand corners.
    ///\param rows Number of years to be collected
    ///\param cols Data columns for the simulation year and for two corner points
    ///\post The `current_row` = -1.
    ///\note Matrix indexing starts from 0 and `current_row` is updated
    ///after each row insertion
  CreateVoxelSpaceData(int rows, int cols=7):
    vspacedata(rows,cols),current_row(-1){};
    ///Insert the lower left and the upper right corner points from VoxelSpace
    ///for the given simulation year.
    ///\param vs VoxelSpace
    ///\param year Simulation yesr
    ///\param interval Frequency to collect VoxelSpace data
    ///\pre Data collected if \f$ \mathit{year} \bmod \mathit{interval} = 0 \f$
    ///\post The `current_row` incremented by 1.
    void insertData(const VoxelSpace& vs,int year, int interval);
    ///Return VoxelSpace data
    ///\retval vspacedata VoxelSpace data
    const TMatrix2D<double>& getData()const{return vspacedata;}
  private:
    TMatrix2D<double> vspacedata;///< VoxelSpace data
    int current_row;///< The current row when inserting VoxelSapce data
  };
}
#endif
