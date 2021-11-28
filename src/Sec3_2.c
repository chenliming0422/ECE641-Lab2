#include <math.h>
#include "tiff.h"
#include "allocate.h"
#include "randlib.h"
#include "typeutil.h"
#include "qGGMRF.h"

int main (int argc, char **argv) 
{
    FILE *fp;
    struct TIFF_img input_img, output_img;
    double **img1, **img2, **y, **x, **e;
    int32_t pixel;
    double sigma_x, sigma_w;
    double cost[20] = {0};
    int32_t i,j, iter;
    double v, theta_1, theta_2;
    int M,N;
    double filter[5][5] = {{1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {3/81.0, 6/81.0, 9/81.0, 6/81.0, 3/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0}};
    double prediction[3][3] = {{1/12.0, 1/6.0, 1/12.0}, {1/6.0, 0.0, 1/6.0}, {1/12.0, 1/6.0, 1/12.0}};

    if ( argc != 2 ) error( argv[0] );

      /* open image file */
    if ( ( fp = fopen ( argv[1], "rb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file %s\n", argv[1] );
        exit(1);
    }

    /* read image */
    if ( read_TIFF ( fp, &input_img ) ) {
        fprintf ( stderr, "error reading file %s\n", argv[1] );
        exit(1);
    }

    /* close image file */
     fclose ( fp );

    /* check the type of image data */
    if ( input_img.TIFF_type != 'g' ) {
        fprintf ( stderr, "error:  image must be gray scale\n" );
        exit ( 1 );
    }

    /* Allocate image of double precision floats */
    img1 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    img2 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    y = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    x = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    e = (double **)get_img(input_img.width,input_img.height,sizeof(double));


    M = input_img.height;
    N = input_img.width;

    /* copy image pixels to double array */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img1[i][j] = input_img.mono[i][j];
    }

    /* blurring fitler */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img2[i][j] = circ_conv2d(img1, M, N, i, j, filter);
    }

    /* Add noise */
    srandom2(1);
    
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        y[i][j] = img2[i][j] + 4 * normal();
    }

    get_TIFF ( &output_img, input_img.height, input_img.width, 'g' );
    
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)y[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;        
        output_img.mono[i][j] = pixel;
        y[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "Y.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file Y.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* MAP estimation using ICD*/
   sigma_w = 4.0;
   sigma_x = 16.95;

    /* Initialize x <- y */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* Initialize e <- y - Hx*/
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        e[i][j] = y[i][j] - circ_conv2d(x, M, N, i, j, filter);
    }

    theta_2 = pow((9/81.0), 2) + 8 * pow((2/81.0), 2) + 4 * (pow((6/81.0), 2) + pow((4/81.0), 2) + pow((3/81.0), 2) + pow((1/81.0), 2));
    theta_2 = theta_2 / (sigma_w * sigma_w);
}