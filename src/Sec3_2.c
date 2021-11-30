#include <math.h>
#include "lib\tiff.h"
#include "lib\allocate.h"
#include "lib\randlib.h"
#include "lib\typeutil.h"
#include "lib\qGGMRF.h"

void error(char *name);
double prior_model(double **img, int M, int N, int u, int v, double g[3][3]);
double circ_conv2d(double **img, int M, int N, int u, int v, double h[5][5]);
void update_error(double **e, int M, int N, int u, int v, double num, double h[5][5]);
double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double p, double q, double h[5][5], double g[3][3]);
double sum_bsr_xsxr(double **img, int M, int N, int u, int v, double sigma_x, double p, double q, double g[3][3]);
double sum_bsr(double **img, int M, int N, int u, int v, double sigma_x, double p, double q, double g[3][3]);

int main (int argc, char **argv) 
{
    FILE *fp;
    struct TIFF_img input_img, output_img;
    double **img1, **img2, **y, **x, **e;
    int32_t pixel;
    double sigma_x, sigma_w;
    double cost[20] = {0};
    int32_t i,j, iter;
    double alpha_star, theta_1, theta_2, h_sum;
    int M,N;
    double h[5][5] = {{1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {3/81.0, 6/81.0, 9/81.0, 6/81.0, 3/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0}};
    double g[3][3] = {{1/12.0, 1/6.0, 1/12.0}, 
                               {1/6.0, 0.0, 1/6.0}, 
                               {1/12.0, 1/6.0, 1/12.0}};

    double p = 1.2;
    double q = 2.0;

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
        img2[i][j] = circ_conv2d(img1, M, N, i, j, h);
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
    if ( ( fp = fopen ( "../output/Y_subsec2.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file Y.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }


/********************************************************************************************************************************************************/
    /* MAP estimation using ICD*/
   sigma_w = 4.0;
   sigma_x = 6.259;

    /* Initialize x <- y */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* Initialize e <- y - Hx*/
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        e[i][j] = y[i][j] - circ_conv2d(x, M, N, i, j, h);
    }

    h_sum = pow((9/81.0), 2) + 8 * pow((2/81.0), 2) + 4 * (pow((6/81.0), 2) + pow((4/81.0), 2) + pow((3/81.0), 2) + pow((1/81.0), 2));
    h_sum = h_sum / (sigma_w * sigma_w);

    for (iter=0; iter<20; iter++){
        printf("start iteration %02d...", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            theta_1 = (-1.0) * circ_conv2d(e, M, N, i, j, h) / (sigma_w * sigma_w) + sum_bsr_xsxr(x, M, N, i, j, sigma_x, p, q, g);
            theta_2 = h_sum + sum_bsr(x, M, N, i, j, sigma_x, p, q, g);

            alpha_star = (-1.0) * theta_1 / theta_2;
            if (alpha_star < (-1.0) * x[i][j]){
                alpha_star = (-1.0) * x[i][j];
            }

            x[i][j] += alpha_star;

            update_error(e, M, N, i, j, alpha_star, h);
        }

        cost[iter] = cost_function(y, x, M, N, sigma_w, sigma_x, p, q, h, g);
        printf("%.15f\n", cost[iter]);
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = (int32_t)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate1_subsec2.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file X.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

/********************************************************************************************************************************************************/
    sigma_x = 6.259*(1/5.0);  
    printf("change sigma_x to %.3f, start...\n", sigma_x);
   
    /* Initialize x <- y */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* Initialize e <- y - Hx*/
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        e[i][j] = y[i][j] - circ_conv2d(x, M, N, i, j, h);
    }

    h_sum = pow((9/81.0), 2) + 8 * pow((2/81.0), 2) + 4 * (pow((6/81.0), 2) + pow((4/81.0), 2) + pow((3/81.0), 2) + pow((1/81.0), 2));
    h_sum = h_sum / (sigma_w * sigma_w);

    for (iter=0; iter<20; iter++){
        printf("start iteration %02d...\n", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            theta_1 = (-1.0) * circ_conv2d(e, M, N, i, j, h) / (sigma_w * sigma_w) + sum_bsr_xsxr(x, M, N, i, j, sigma_x, p, q, g);
            theta_2 = h_sum + sum_bsr(x, M, N, i, j, sigma_x, p, q, g);

            alpha_star = (-1.0) * theta_1 / theta_2;
            if (alpha_star < (-1.0) * x[i][j]){
                alpha_star = (-1.0) * x[i][j];
            }

            x[i][j] += alpha_star;

            update_error(e, M, N, i, j, alpha_star, h);
        }
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = (int32_t)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate2_subsec2.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file X.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

/***********************************************************************************************************************************************************************/
    sigma_x = 6.259*(5.0);  
    printf("change sigma_x to %.3f, start...\n", sigma_x);
   
    /* Initialize x <- y */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* Initialize e <- y - Hx*/
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        e[i][j] = y[i][j] - circ_conv2d(x, M, N, i, j, h);
    }

    h_sum = pow((9/81.0), 2) + 8 * pow((2/81.0), 2) + 4 * (pow((6/81.0), 2) + pow((4/81.0), 2) + pow((3/81.0), 2) + pow((1/81.0), 2));
    h_sum = h_sum / (sigma_w * sigma_w);

    for (iter=0; iter<20; iter++){
        printf("start iteration %02d...\n", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            theta_1 = (-1.0) * circ_conv2d(e, M, N, i, j, h) / (sigma_w * sigma_w) + sum_bsr_xsxr(x, M, N, i, j, sigma_x, p, q, g);
            theta_2 = h_sum + sum_bsr(x, M, N, i, j, sigma_x, p, q, g);

            alpha_star = (-1.0) * theta_1 / theta_2;
            if (alpha_star < (-1.0) * x[i][j]){
                alpha_star = (-1.0) * x[i][j];
            }

            x[i][j] += alpha_star;

            update_error(e, M, N, i, j, alpha_star, h);
        }
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = (int32_t)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate3_subsec2.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file X.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

/*******************************************************************************************/
    /* de-allocate space which was used for the images */   
    free_TIFF ( &(input_img) );
    free_TIFF ( &(output_img) );

    free_img( (void**)img1 );
    free_img( (void**)img2 );
    free_img( (void**)y );
    free_img( (void**)x );
    free_img( (void**)e );

    return(0);

}

void error(char *name)
{
    printf("usage:  %s  image.tiff \n\n",name);
    printf("this program reads in a 24-bit color TIFF image.\n");
    printf("It then horizontally filters the green component, adds noise,\n");
    printf("and writes out the result as an 8-bit image\n");
    printf("with the name 'green.tiff'.\n");
    printf("It also generates an 8-bit color image,\n");
    printf("that swaps red and green components from the input image");
    exit(1);
}

double sum_bsr_xsxr(double **img, int M, int N, int u, int v, double sigma_x, double p, double q, double g[3][3])
{
    int i, j;
    double delta;
    double T = 1.0;
    double sum = 0.0;
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
        delta = img[u][v] - img[(u+i-1+M)%M][(v+j-1+N)%N];
        sum += get_btilde(delta, g[i][j], sigma_x, p, q, T) * delta;
    }
    return sum;

}

double sum_bsr(double **img, int M, int N, int u, int v, double sigma_x, double p, double q, double g[3][3])
{
    int i, j;
    double delta;
    double T = 1.0;
    double sum = 0.0;
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
        delta = img[u][v] - img[(u+i-1+M)%M][(v+j-1+N)%N];
        sum += get_btilde(delta, g[i][j], sigma_x, p, q, T);
    }
    return sum;
}


double prior_model(double **img, int M, int N, int u, int v, double g[3][3])
{
    int i, j;
    double val = 0.0;
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
        val += (img[(u+i-1+M)%M][(v+j-1+N)%N]) * g[i][j];
    }
    return val;
}

double circ_conv2d(double **img, int M, int N, int u, int v, double h[5][5])
{
    int i, j;
    double val = 0.0;
    for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++){
        val += (img[(u+i-2+M)%M][(v+j-2+N)%N]) * h[i][j];
    }
    return val;
}

void update_error(double **e, int M, int N, int u, int v, double num, double h[5][5])
{
    int i, j;
    for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++) {
        e[(u+i-2+M)%M][(v+j-2+N)%N] -= h[i][j] * num;
    }
}

double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double p, double q, double h[5][5], double g[3][3])
{
    int i, j, m, n;
    double tmp_err;
    double delta;
    double tmp_val1 = 0.0;
    double tmp_val2 = 0.0;
    double T = 1.0;
    double cost;
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++){
        tmp_err = y[i][j] - circ_conv2d(x, M, N, i, j, h);
        tmp_val1 += pow(tmp_err, 2.0) ;

        for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++) {
            delta = x[i][j] - x[(i+m-1+M)%M][(j+n-1+N)%N];
            tmp_val2 += get_rho(delta, g[m][n], sigma_x, p, q, T) * g[m][n];
        }
    }
    cost = tmp_val1 / (2.0 * sigma_w * sigma_w) + tmp_val2 / 2.0;
    return cost;
}