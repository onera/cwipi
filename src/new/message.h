#ifndef __MESSAGE_H__
#define __MESSAGE_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

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

/*
  This file is inspired from OpenPALM.
  OpenPALM is a free software under the GNU Lesser General Public License.
  See: https://www.cerfacs.fr/globc/PALM_WEB/
*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_MSG_MAXMSGSIZE 8192

/* Request numbering for CWIPI mirror operations */
#define CWP_MSG_DIE                                 0

#define CWP_MSG_CWP_INIT                            1
#define CWP_MSG_CWP_FINALIZE                        2

#define CWP_MSG_CWP_PARAM_LOCK                      3
#define CWP_MSG_CWP_PARAM_UNLOCK                    4
#define CWP_MSG_CWP_PARAM_ADD                       5
#define CWP_MSG_CWP_PARAM_GET                       6
#define CWP_MSG_CWP_PARAM_SET                       7
#define CWP_MSG_CWP_PARAM_DEL                       8
#define CWP_MSG_CWP_PARAM_N_GET                     9
#define CWP_MSG_CWP_PARAM_LIST_GET                  10
#define CWP_MSG_CWP_PARAM_IS                        11
#define CWP_MSG_CWP_PARAM_REDUCE                    12

#define CWP_MSG_CWP_CPL_CREATE                      13
#define CWP_MSG_CWP_CPL_DEL                         14

#define CWP_MSG_CWP_PROPERTIES_DUMP                 15
#define CWP_MSG_CWP_VISU_SET                        16
#define CWP_MSG_CWP_STATE_UPDATE                    17
#define CWP_MSG_CWP_TIME_UPDATE                     18
#define CWP_MSG_CWP_STATE_GET                       19
#define CWP_MSG_CWP_CODES_NB_GET                    20
#define CWP_MSG_CWP_CODES_LIST_GET                  21
#define CWP_MSG_CWP_LOC_CODES_NB_GET                22
#define CWP_MSG_CWP_LOC_CODES_LIST_GET              23
#define CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET           24
#define CWP_MSG_CWP_UNCOMPUTED_TGTS_GET             25
#define CWP_MSG_CWP_N_COMPUTED_TGTS_GET             26
#define CWP_MSG_CWP_COMPUTED_TGTS_GET               27
#define CWP_MSG_CWP_N_INVOLVED_SRCS_GET             28
#define CWP_MSG_CWP_INVOLVED_SRCS_GET               29
#define CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE  30
#define CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET     31
#define CWP_MSG_CWP_USER_TGT_PTS_SET                32

#define CWP_MSG_CWP_MESH_INTERF_FINALIZE            33
#define CWP_MSG_CWP_MESH_INTERF_VTX_SET             34
#define CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD           35
#define CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET       36
#define CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET       37
#define CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET    38
#define CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET    39
#define CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET    40
#define CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET    41
#define CWP_MSG_CWP_MESH_INTERF_DEL                 42
#define CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET   43
#define CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET   44

/* Init an request */
#define NEWMESSAGE(msg,msg_type) {memset(&msg,0,sizeof(t_message));msg.message_type=msg_type;}

/*============================================================================
 * Types definition
 *============================================================================*/

/* Data structure used for each client-server request */

typedef struct t_message
{
  int message_type;
  int flag;
  int msg_tag;
  int msg_time;
  int data_size;
        int data1;
        int data2;
        int data3;
  char space[256];
  char object[256];
}t_message,*p_message;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MESSAGE_H__ */
