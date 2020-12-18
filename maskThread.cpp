
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



#include "maskThread.hpp"

maskThread::maskThread (QObject *parent)
  : QThread(parent)
{
}



maskThread::~maskThread ()
{
}



void maskThread::mask (uint8_t *bl, int32_t np, int32_t *pc, double **py, double **px, double slt, double sln, double xr, double yr,
                       int32_t dx, int32_t dy, int32_t sx, int32_t sy, int32_t w, uint8_t wtr, uint8_t *c, int32_t p)
{
  QMutexLocker locker (&mutex);

  l_block = bl;
  l_num_poly = np;
  l_poly_count = pc;
  l_poly_y = py;
  l_poly_x = px;
  l_sw_lat = slt;
  l_sw_lon = sln;
  l_x_res = xr;
  l_y_res = yr;
  l_x_dim = dx;
  l_y_dim = dy;
  l_start_x = sx;
  l_start_y = sy;
  l_width = w;
  l_water = wtr;
  l_complete = c;
  l_pass = p;

  if (!isRunning ()) start ();
}



void maskThread::run ()
{
  mutex.lock ();

  uint8_t *block = l_block;
  int32_t num_poly = l_num_poly;
  int32_t *poly_count = l_poly_count;
  double **poly_y = l_poly_y;
  double **poly_x = l_poly_x;
  double sw_lat = l_sw_lat;
  double sw_lon = l_sw_lon;
  double x_res = l_x_res;
  double y_res = l_y_res;
  int32_t x_dim = l_x_dim;
  int32_t y_dim = l_y_dim;
  int32_t start_x = l_start_x;
  int32_t start_y = l_start_y;
  int32_t width = l_width;
  uint8_t water = l_water;
  uint8_t *complete = l_complete;
  int32_t pass = l_pass;

  mutex.unlock ();


  //qDebug () << __LINE__ << pass;

  int32_t percent = 0, old_percent = -1;


  int32_t end_x = start_x + x_dim;
  int32_t end_y = start_y + y_dim;


  //  Latitude loop.

  for (int32_t i = start_y ; i < end_y ; i++)
    {
      //  Compute the latitude of the center of the "y_res" sized bin (that's why we add 0.5).

      double slat = (double) sw_lat + ((double) i + 0.5) * y_res;


      //  Longitude loop.

      for (int32_t j = start_x ; j < end_x ; j++)
        {
          //  Compute the longitude of the center of the "x_res" sized bin (that's why we add 0.5).

          double slon = (double) sw_lon + ((double) j + 0.5) * x_res;


          int32_t inside_count = 0;


          //  Check against all polygons.

          for (int32_t k = 0 ; k < num_poly ; k++)
            {
              if (inside_polygon2 (poly_x[k], poly_y[k], poly_count[k], slon, slat)) inside_count++;
            }


          //  Set the flag (NVTrue for land, NVFalse for water).

          if (inside_count % 2)
            {
              if (water)
                {
                  block[i * width + j] = NVFalse;
                }
              else
                {
                  block[i * width + j] = NVTrue;
                }
            }
          else
            {
              if (water)
                {
                  block[i * width + j] = NVTrue;
                }
              else
                {
                  block[i * width + j] = NVFalse;
                }
            }
        }


      percent = (int32_t) (((double) (i - start_y) / (double) y_dim) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Pass %d - %03d%% processed\r", pass, percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  //qDebug () << __LINE__ << pass;

  complete[pass] = NVTrue;
}
