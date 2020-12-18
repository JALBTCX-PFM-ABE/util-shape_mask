
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


/*******************************************************************************************/
/*!

  - Module Name:        shape_mask

  - Programmer(s):      Jan C. Depner (PFM Software)

  - Date Written:       February 2014

  - Purpose:            Reads a shape file containing land or water polygons and creates a
                        land mask at the specified resolution.

  - Arguments:          argv[1]         -   resolution (in integer meters - 1, 2, 3...)

*********************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#include "nvutility.hpp"
#include "nvutility.h"

#include "shapefil.h"
#include "version.h"

#include "maskThread.hpp"


void usage (char *string)
{
  fprintf (stderr, "Program: %s\n", string);
  fprintf (stderr, "Purpose: Reads a shape file containing land or water polygons\n");
  fprintf (stderr, "and a surrounding generic area file and creates a land mask at\n");
  fprintf (stderr, "the specified resolution.  The output file is ALWAYS a land mask\n");
  fprintf (stderr, "regardless of whether the shape file contains land or water\n");
  fprintf (stderr, "polygons.  The output file can be used in cpfFilter to invalidate\n");
  fprintf (stderr, "land processed shots over water or water processed shots over land.\n\n");
  fprintf (stderr, "Usage: %s SHAPEFILE_NAME RESOLUTION [-w]\n\n", string);
  fprintf (stderr, "Where\n");
  fprintf (stderr, "\tSHAPEFILE_NAME = name of shape (.shp) file containing land/water polygons\n");
  fprintf (stderr, "\tRESOLUTION = resolution of output mask in integer meters (1, 2, 3...)\n");
  fprintf (stderr, "\t-w = set this if shape file polygons contain water areas instead of land areas\n\n");
  fprintf (stderr, "Caveats:\n");
  fprintf (stderr, "\tThe shapefile must contain complete polygons for all land (or water) areas needed.\n");
  fprintf (stderr, "\tIn addition to the .shp file there MUST be a generic area file (.are) with the same\n");
  fprintf (stderr, "\tname (e.g. fred.shp, fred.are) that defines the entire area to be covered by the mask.\n\n");
  exit (-1);
}



int32_t main (int32_t argc, char **argv)
{
  int32_t           type, numShapes, resolution = 0, num_threads = 4;
  double            minBounds[4], maxBounds[4];
  char              shpname[512], arename[512], mskname[512], c;
  FILE              *ofp;
  maskThread        mask_thread[4];
  extern int        optind;


  printf ("\n\n%s\n\n", VERSION);


  uint8_t water = NVFalse;
  while ((c = getopt (argc, argv, "w")) != EOF)
    {
      switch (c)
        {
	case 'w':
	  water = NVTrue;
	  break;

	default:
	  usage (argv[0]);
	  break;
        }
    }


  /*  Make sure we got the mandatory arguments.  */

  if (optind >= argc) usage (argv[0]);


  strcpy (shpname, argv[optind]);
  sscanf (argv[optind + 1], "%d", &resolution);


  strcpy (arename, shpname);
  strcpy (strstr (arename, ".shp"), ".are");


  int32_t pc;
  double px[5], py[5];
  NV_F64_XYMBR mbr;
  if (!get_area_mbr (arename, &pc, px, py, &mbr)) usage (argv[0]);


  double **poly_x = NULL;
  double **poly_y = NULL;
  int32_t *poly_count = NULL;
  uint8_t *block = NULL;
  SHPHandle shpHandle;
  SHPObject *shape = NULL;


  //  Initialize variables

  int32_t num_poly = -1;


  //  Open shape file

  shpHandle = SHPOpen (shpname, "rb");

  if (shpHandle == NULL)
    {
      perror (shpname);
      exit (-1);
    }


  fprintf (stderr, "Reading %s                        \n", shpname);
  fflush (stderr);


  //  Get shape file header info

  SHPGetInfo (shpHandle, &numShapes, &type, minBounds, maxBounds);


  //  Convert the resolution (approximately) to decimal degrees at the center of the MBR.

  NV_F64_COORD2 center, xy;
  center.x = mbr.min_x + (mbr.max_x - mbr.min_x) / 2.0;
  center.y = mbr.min_y + (mbr.max_y - mbr.min_y) / 2.0;


  newgp (center.y, center.x, 90.0, (double) resolution, &xy.y, &xy.x);

  double x_resolution = xy.x - center.x;

  newgp (center.y, center.x, 0.0, (double) resolution, &xy.y, &xy.x);

  double y_resolution = xy.y - center.y;


  //  Adjust the MBR to match the computed resolutions.

  int32_t range_x = (int32_t) ((mbr.max_x - mbr.min_x) / x_resolution) + 1;
  int32_t range_y = (int32_t) ((mbr.max_y - mbr.min_y) / y_resolution) + 1;


  int32_t half_range_x = range_x / 2;
  int32_t half_range_y = range_y / 2;

  mbr.min_x = center.x - (double) half_range_x * x_resolution;
  mbr.max_x = center.x + (double) half_range_x * x_resolution;
  mbr.min_y = center.y - (double) half_range_y * y_resolution;
  mbr.max_y = center.y + (double) half_range_y * y_resolution;

  int32_t dim_x = NINT ((mbr.max_x - mbr.min_x) / x_resolution);
  int32_t dim_y = NINT ((mbr.max_y - mbr.min_y) / y_resolution);


  //  Read all shapes

  for (int32_t i = 0 ; i < numShapes ; i++)
    {
      shape = SHPReadObject (shpHandle, i);


      //  Get all vertices

      if (shape->nVertices >= 2)
        {
          for (int32_t j = 0, numParts = 1 ; j < shape->nVertices ; j++)
            {
              uint8_t start_segment = NVFalse;


              //  Check for start of a new segment.

              if (!j && shape->nParts > 0) start_segment = NVTrue;


              //  Check for the start of a new segment inside a larger group of points (this would be a "Ring" point).

              if (numParts < shape->nParts && shape->panPartStart[numParts] == j)
                {
                  start_segment = NVTrue;
                  numParts++;
                }


              //  Start a new segment

              if (start_segment)
                {
                  //  Since num_poly starts at -1 this is perfectly cool.

                  num_poly++;


                  //  Allocate the count array.

                  poly_count = (int32_t *) realloc (poly_count, (num_poly + 1) * sizeof (int32_t));
                  if (poly_count == NULL)
                    {
                      perror ("Allocating poly_count memory");
                      exit (-1);
                    }


                  //  Set the count for the new arrays to zero.

                  poly_count[num_poly] = 0;


                  //  Allocate the polygon arrays.

                  poly_x = (double **) realloc (poly_x, (num_poly + 1) * sizeof (double *));
                  if (poly_x == NULL)
                    {
                      perror ("Allocating poly_x memory");
                      exit (-1);
                    }
                  poly_x[num_poly] = NULL;


                  poly_y = (double **) realloc (poly_y, (num_poly + 1) * sizeof (double *));
                  if (poly_y == NULL)
                    {
                      perror ("Allocating poly_y memory");
                      exit (-1);
                    }
                  poly_y[num_poly] = NULL;
                }


              //  Allocate memory for the new point.

              poly_x[num_poly] = (double *) realloc (poly_x[num_poly], (poly_count[num_poly] + 1) * sizeof (double));
              if (poly_x[num_poly] == NULL)
                {
                  perror ("Allocating poly_x[num_poly] memory");
                  exit (-1);
                }

              poly_y[num_poly] = (double *) realloc (poly_y[num_poly], (poly_count[num_poly] + 1) * sizeof (double));
              if (poly_y[num_poly] == NULL)
                {
                  perror ("Allocating poly_y[num_poly] memory");
                  exit (-1);
                }


              //  Add point to current segment

              poly_x[num_poly][poly_count[num_poly]] = shape->padfX[j];
              poly_y[num_poly][poly_count[num_poly]] = shape->padfY[j];


              //  Increment the point counter.

              poly_count[num_poly]++;
            }
        }


      //  Destroy the shape object.

      SHPDestroyObject (shape);
    }


  //  Increment num_poly to account for the last polygon.

  num_poly++;


  //  Close the input file.

  SHPClose (shpHandle);


  //  Allocate the uint8_t block to put the land/water flags into.

  block = (uint8_t *) calloc (dim_x * dim_y, sizeof (uint8_t));

  if (block == NULL)
    {
      perror ("Allocating block memory");
      exit (-1);
    }


  //  Start all "num_threads" threads to compute the mask.

  uint8_t *complete = (uint8_t *) calloc (num_threads, sizeof (uint8_t));
  if (complete == NULL)
    {
      perror ("Allocating complete memory");
      exit (-1);
    }

  mask_thread[0].mask (block, num_poly, poly_count, poly_y, poly_x, mbr.min_y, mbr.min_x, x_resolution, y_resolution, half_range_x, half_range_y, 0, 0,
                       dim_x, water, complete, 0);
  mask_thread[1].mask (block, num_poly, poly_count, poly_y, poly_x, mbr.min_y, mbr.min_x, x_resolution, y_resolution, half_range_x, dim_y - half_range_y, 0,
                       half_range_y, dim_x, water, complete, 1);
  mask_thread[2].mask (block, num_poly, poly_count, poly_y, poly_x, mbr.min_y, mbr.min_x, x_resolution, y_resolution, dim_x - half_range_x, half_range_y,
                       half_range_x, 0, dim_x, water, complete, 2);
  mask_thread[3].mask (block, num_poly, poly_count, poly_y, poly_x, mbr.min_y, mbr.min_x, x_resolution, y_resolution, dim_x - half_range_x,
                       dim_y - half_range_y, half_range_x, half_range_y, dim_x, water, complete, 3);


  //  We can't move on until all of the threads are complete.

  for (int32_t i = 0 ; i < num_threads ; i++)
    {
      mask_thread[i].wait ();
    }


  //  Free all of the polygon memory.

  for (int32_t i = 0 ; i < num_poly ; i++)
    {
      if (poly_x[i] != NULL) free (poly_x[i]);
      if (poly_y[i] != NULL) free (poly_y[i]);
    }


  if (poly_x != NULL) free (poly_x);
  if (poly_y != NULL) free (poly_y);
  if (poly_count != NULL) free (poly_count);
  poly_x = NULL;
  poly_y = NULL;
  poly_count = NULL;


  strcpy (mskname, shpname);
  strcpy (strstr (mskname, ".shp"), ".msk");


  //  Open the output file.

  if ((ofp = fopen (mskname, "wb")) == NULL)
    {
      perror (mskname);
      exit (-1);
    }


  //  Just for fun, compute the difference in longitudinal size at the north and south of the area.

  double dist_n, dist_s, az, diff;
  invgp (NV_A0, NV_B0, mbr.max_y, mbr.min_x, mbr.max_y, mbr.min_x + x_resolution, &dist_n, &az);
  invgp (NV_A0, NV_B0, mbr.min_y, mbr.min_x, mbr.min_y, mbr.min_x + x_resolution, &dist_s, &az);
  diff = dist_n - dist_s;
  //fprintf (stderr, "%s %s %d %.5f %.5f %.8f\n",NVFFL,dist_n,dist_s,diff);



  //  Write the (minimalist) ASCII header.

  time_t t = time (&t);
  struct tm *cur_tm = gmtime (&t);

  fprintf (ofp, "[HEADER SIZE] = %d\n", 16384);
  fprintf (ofp, "[VERSION] = %s\n", VERSION);
  fprintf (ofp, "[CREATION DATE] = %s", asctime (cur_tm));
  fprintf (ofp, "[START LAT] = %.11f\n", mbr.min_y);
  fprintf (ofp, "[START LON] = %.11f\n", mbr.min_x);
  fprintf (ofp, "[LAT RESOLUTION] = %.11f\n", y_resolution);
  fprintf (ofp, "[LON RESOLUTION] = %.11f\n", x_resolution);
  fprintf (ofp, "[HEIGHT] = %d\n", dim_y);
  fprintf (ofp, "[WIDTH] = %d\n", dim_x);
  fprintf (ofp, "[NOMINAL BIN SIZE IN METERS] = %d\n", resolution);
  fprintf (ofp, "[NORTH SOUTH LON BIN SIZE DIFFERENCE IN METERS] = %.08f\n", diff);
  fprintf (ofp, "[END OF HEADER]\n");


  //  Zero out the remainder of the header.

  uint8_t zero = 0;
  int32_t j = ftell (ofp);
  for (int32_t i = j ; i < 16384 ; i++) fwrite (&zero, 1, 1, ofp);


  //FILE *fp = fopen ("fred.yxz", "w");

  for (int32_t i = 0 ; i < dim_y ; i++)
    {
      //double slat = mbr.min_y + ((double) i + 0.5) * y_resolution;

      for (int32_t j = 0 ; j < dim_x ; j++)
        {
          //double slon = mbr.min_x + ((double) j + 0.5) * x_resolution;

          //if (block[i * dim_x + j]) fprintf (fp, "%.11f,%.11f,1\n", slat, slon); 

          fwrite (&block[i * dim_x + j], 1, 1, ofp);
        }
    }

  //fclose (fp);

  fclose (ofp);

  free (block);


  fprintf (stderr, "100%% processed                         \n\n");
  fflush (stderr);

  return (0);
}
