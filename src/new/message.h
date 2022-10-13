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
#define CWP_MSG_DIE                             0

#define CWP_MSG_CWP_INIT                        1
#define CWP_MSG_CWP_FINALIZE                    2

#define CWP_MSG_CWP_PARAM_LOCK                  3
#define CWP_MSG_CWP_PARAM_UNLOCK                4
#define CWP_MSG_CWP_PARAM_ADD                   5
#define CWP_MSG_CWP_PARAM_GET                   6

#define CWP_MSG_CWP_CPL_CREATE                  7
#define CWP_MSG_CWP_CPL_DEL                     8

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
