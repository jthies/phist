#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <png.h>

#include "scamac_collection.h" // because of scamac_err_desc ... remove this include
#include "scamac_aux.h"
#include "scamac_plot.h"

// buffer size = (3+1)*width*height, including alpha channel
// storage is "rows first"!
static int write_png(const char *filename, int width, int height, unsigned char * buffer, char* title) {
  int code = 0;
  FILE *fp = NULL;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  // Open file for writing (binary mode)
  fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Could not open file %s for writing\n", filename);
    code = 1;
    goto finalise;
  }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fprintf(stderr, "Could not allocate write struct\n");
    code = 1;
    goto finalise;
  }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fprintf(stderr, "Could not allocate info struct\n");
    code = 1;
    goto finalise;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "Error during png creation\n");
    code = 1;
    goto finalise;
  }

  png_init_io(png_ptr, fp);

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, info_ptr, width, height,
               8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Set title
  if (title != NULL) {
    png_text title_text;
    title_text.compression = PNG_TEXT_COMPRESSION_NONE;
    title_text.key = "Title";
    title_text.text = title;
    png_set_text(png_ptr, info_ptr, &title_text, 1);
  }

  png_write_info(png_ptr, info_ptr);

  // Write image data
  int y;
  for (y=0; y<height; y++) {
    png_write_row(png_ptr, (png_bytep) &buffer[4*width*y]);
  }

  // End write
  png_write_end(png_ptr, NULL);

finalise:
  if (fp != NULL) fclose(fp);
  if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

  return code;
}


ScamacErrorCode scamac_plot_pattern_onetoone(const scamac_matrix_pattern_st * pat, const char * filename) {
  if (!pat || !filename) {
    return SCAMAC_ENULL;
  }

  int width =pat->px;
  int height=pat->py;
  unsigned char * buffer = malloc(4*width*height * sizeof *buffer);
  if (buffer == NULL) {
    return SCAMAC_EMALLOCFAIL;
  }

  // transparent "black" background
  memset(buffer, 0, 4*width*height * sizeof *buffer);


  int x,y;
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      if (pat->pat[width*y+x]) {
        unsigned char * pixel = &buffer[4*(width*y+x)];
        pixel[0]=255;
        pixel[1]=0;
        pixel[2]=0;
        pixel[3]=255;
      }
    }
  }

  int result = write_png(filename, width, height, buffer, "Sparsity pattern");
  if (result) {
    return SCAMAC_EFAIL;
  }

  free(buffer);

  return SCAMAC_EOK;
}


ScamacErrorCode scamac_plot_pattern(const scamac_matrix_pattern_st * pat, int downscale, const char * filename) {
  if (!pat || !filename) {
    return SCAMAC_ENULL;
  }

  int bufnx = pat->px;
  int bufny = pat->py;

  // convert buffer to image
  int width = 1500;
  int height = 1500;
  unsigned char * image = malloc( 4*width*height * sizeof *image);

  const int rectsize = 15;


  // transparent "black" background
  //  memset(image, 0, 4*width*height * sizeof * image);

  int ix,iy;

  // transparent "white" background
  for (iy=0; iy<width; iy++) {
    for (ix=0; ix<height; ix++) {
      unsigned char * pixel = &image[4*(width*iy + ix)];
      pixel[0]=255;
      pixel[1]=255;
      pixel[2]=255;
      pixel[3]=0;
    }
  }


  for (iy=0; iy<bufny; iy++) {
    int y = floor (  (double) iy * (double) (height-rectsize) / (double) (bufny) );
    for (ix=0; ix<bufnx; ix++) {
      if (pat->pat[iy*bufnx+ix]) {
        int x = floor (  (double) ix * (double) (height-rectsize) / (double) (bufnx) );
        int irx,iry;
        for (iry=0; iry<rectsize; iry++) {
          for (irx=0; irx<rectsize; irx++) {
            unsigned char * pixel = &image[4*(width*(y+iry)+(x+irx))];
            //pixel[0]=255;
            //pixel[1]=0;
            pixel[0]=floor(255.0*(1.0-fabs((double) ix/bufnx- (double) iy/bufny)));
            pixel[1]=floor(255.0* 0.5*fabs((double) ix/bufnx+ (double) iy/bufny));
            pixel[2]=0;
            pixel[3]=255;
          }
        }
      }
    }
  }

  // downscale image
  int width_small;
  int height_small;
  unsigned char * image_small;

  if (downscale == 1) {
    width_small = width;
    height_small = height;
    image_small = image;
  } else {
    width_small = floor((width+1)/downscale -1 );
    height_small = floor((height+1)/downscale -1 );
    image_small = malloc( 4*width_small*height_small * sizeof *image_small);
    double * sclfct = malloc( (2*downscale-1) * sizeof *sclfct);
    sclfct[downscale-1]=1.0;
    for (ix=1; ix<downscale; ix++) {
      sclfct[ downscale -1 -ix ] = 1.0 - ( (double) ix / (double) (downscale) );
      sclfct[ downscale -1 +ix ] = 1.0 - ( (double) ix / (double) (downscale) );
    }
    double sum=0.0;
    for (ix=0; ix<2*downscale-1; ix++) {
      sum = sum + sclfct[ix];
    }
    for (ix=0; ix<2*downscale-1; ix++) {
      sclfct[ix] = sclfct[ix]/sum;
    }
    for (iy=0; iy<height_small; iy++) {
      for (ix=0; ix<width_small; ix++) {
        double r=0.0,g=0.0,b=0.0, al=0.0;
        int irx,iry;
        for (iry=0; iry<2*downscale-1; iry++) {
          for (irx=0; irx<2*downscale-1; irx++) {
            unsigned char * pixel = &image[4*(width*(downscale*iy+iry)+(downscale*ix+irx))];

            r = r + sclfct[irx]*sclfct[iry] * (double) pixel[0];
            g = g + sclfct[irx]*sclfct[iry] * (double) pixel[1];
            b = b + sclfct[irx]*sclfct[iry] * (double) pixel[2];
            al = al + sclfct[irx]*sclfct[iry] * (double) pixel[3];
          }
        }
        unsigned char * pixel = &image_small[4*(width_small*iy+ix)];
        pixel[0]=floor(r);
        pixel[1]=floor(g);
        pixel[2]=floor(b);
        pixel[3]=al;
      }
    }
    free(sclfct);
    free(image);
    image=NULL;
  }




  int result = write_png(filename, width_small, height_small, image_small, "Sparsity pattern");
  if (result) {
    return SCAMAC_EFAIL;
  }

  free(image_small);

  return SCAMAC_EOK;
}
