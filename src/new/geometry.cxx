/*
  This file is part of the CWIPI library.

  Copyright (C) 2012  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <map>
#include <geometry.hxx>
#include <mpi.h>

#include <pdm_mesh_nodal.h>
#include <pdm_dist_cloud_surf.h>
#include <pdm_gnum.h>
#include <pdm_timer.h>
#include <pdm_gnum_location.h>
#include <bftc_error.h>
#include <bftc_printf.h>
#include "cwp.h"
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>

namespace cwipi {

  Geometry::Geometry()
  {
  }


  Geometry::~Geometry()
  {
  }
  
  

/***************************************************************************/
/***************************************************************************/



  void Geometry::_IAlltoallIndexSend(void* send_buffer,
                                     int* send_count,
                                     int* send_disp,
                                     MPI_Datatype type, 
                                     MPI_Comm comm,
                                     std::vector<int> connectableRanks
                                    ){
      int comm_size;
      MPI_Comm_size(comm,&comm_size);
      int nranks = connectableRanks.size();

      _send_requests.resize(nranks,0);

      int tagsend = -1;
      if(localName == _codeVector[0]) {
        tagsend =0;
      }
      else {
       tagsend =1; 
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        
        if(type == MPI_BYTE){
           int* sendptr = (int*)send_buffer;
          //  target_data* sendptr = (target_data*)send_buffer;
            MPI_Issend(&(  sendptr [ send_disp[distant_rank]/sizeof(int) /*recv_size[distant_rank]*/ ] ), send_count[distant_rank] * 1, type, distant_rank, tagsend,
                   comm,
                   &(_send_requests[i_rank]));        
        }   
      }//end for on i_rank
  }




 void Geometry::_IAlltoallIndexRecv(void* recv_buffer,
                int* recv_count,
                int* recv_disp,
                MPI_Datatype type, 
                MPI_Comm comm,
                std::vector<int> connectableRanks
                ){

      int comm_size;
      MPI_Comm_size(comm,&comm_size);
      int nranks = connectableRanks.size();

      _recv_requests.resize(nranks,0);

      int tagrecv = -1;
      if(localName == _codeVector[0]) {
        tagrecv = 1;
      }
      else {
        tagrecv = 0;
      }

      for(int i_rank=0;i_rank<nranks;i_rank++) {
        int distant_rank = connectableRanks[i_rank];
        if(type == MPI_BYTE){
          int* recvptr = (int*)recv_buffer;
          //target_data* recvptr = (target_data*)recv_buffer;
          MPI_Irecv(&(  recvptr [ recv_disp[distant_rank]/sizeof(int)  /*recv_size[distant_rank]*/ ] ), recv_count[distant_rank] * 1,type, distant_rank, tagrecv,
                    comm,
                    &(_recv_requests[i_rank])); 
        }   
      }//end for on i_rank

  }

/***************************************************************************/
/***************************************************************************/



} // end namespace cwipi



