
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.
*********************************************************************************************/

/****************************************  IMPORTANT NOTE  **********************************

    Comments in this file that start with / * ! or / / ! are being used by Doxygen to
    document the software.  Dashes in these comment blocks are used to create bullet lists.
    The lack of blank lines after a block of dash preceeded comments means that the next
    block of dash preceeded comments is a new, indented bullet list.  I've tried to keep the
    Doxygen formatting to a minimum but there are some other items (like <br> and <pre>)
    that need to be left alone.  If you see a comment that starts with / * ! or / / ! and
    there is something that looks a bit weird it is probably due to some arcane Doxygen
    syntax.  Be very careful modifying blocks of Doxygen comments.

*****************************************  IMPORTANT NOTE  **********************************/



#ifndef MASKTHREAD_H
#define MASKTHREAD_H


#include <QtCore>
#include <QtGui>
#if QT_VERSION >= 0x050000
#include <QtWidgets>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#include "nvutility.h"
#include "nvutility.hpp"


class maskThread:public QThread
{
  Q_OBJECT 


public:

  maskThread (QObject *parent = 0);
  ~maskThread ();

  void mask (uint8_t *bl = NULL, int32_t np = 0, int32_t *pc = NULL, double **py = NULL, double **px = NULL, double slt = 0.0, double sln = 0.0,
             double xr = 0.0, double yr = 0.0, int32_t dx = 0, int32_t dy = 0, int32_t sx = 0, int32_t sy = 0, int32_t w = 0, uint8_t wtr = NVFalse,
             uint8_t *c = NULL, int32_t p = -1);


signals:


protected:


  QMutex           mutex;

  uint8_t          *l_block, *l_complete, l_water;

  int32_t          l_resolution, l_num_poly, *l_poly_count, l_num_threads, l_pass;

  double           **l_poly_y, **l_poly_x, l_sw_lat, l_sw_lon, l_x_res, l_y_res, l_x_dim, l_y_dim, l_start_x, l_start_y, l_width;


  void             run ();


protected slots:

private:
};

#endif
