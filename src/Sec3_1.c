#include <math.h>
#include "lib\tiff.h"
#include "lib\allocate.h"
#include "lib\randlib.h"
#include "lib\typeutil.h"
#include "lib\solve.h"

static double theta_1, theta_2;
static int32_t pixel_index_u, pixel_index_v;
static double sigma_x, sigma_w;
static double p = 1.2;
static int M,N;
static double filter[5][5] = {{1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0},
                              {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                              {3/81.0, 6/81.0, 9/81.0, 6/81.0, 3/81.0},
                              {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                              {1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0}};
static double prediction[3][3] = {{1/12.0, 1/6.0, 1/12.0}, 
                                  {1/6.0, 0.0, 1/6.0}, 
                                  {1/12.0, 1/6.0, 1/12.0}};

void error(char *name);
double prior_model(double **img, int M, int N, int u, int v, double g[3][3]);
double circ_conv2d(double **img, int M, int N, int u, int v, double h[5][5]);
double root_function(double u, void * pblock);
void update_error(double **e, int M, int N, int u, int v, double num, double h[5][5]);
double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double h[5][5], double g[3][3]);
double find_clique_max(double **img, int M, int N, int u, int v);
double find_clique_min(double **img, int M, int N, int u, int v);

int main (int argc, char **argv) 
{
    FILE *fp;
    struct TIFF_img input_img, output_img;
    double **img1, **img2, **y, **x, **e;
    int32_t pixel;
    double cost[20] = {0};
    int32_t i,j, iter;
    double tmp_v, tmp_solve;
    double low, high;

    int err_code;

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

    M = input_img.height;
    N = input_img.width;
    sigma_x = 6.259;
    sigma_w = 4.0;

    /* Allocate image of double precision floats */
    img1 = (double **)get_img(N, M, sizeof(double));
    img2 = (double **)get_img(N, M, sizeof(double));
    y = (double **)get_img(N, M, sizeof(double));
    x = (double **)get_img(N, M, sizeof(double));
    e = (double **)get_img(N, M, sizeof(double));

    /* copy image pixels to double array */
    for ( i = 0; i < M; i++ )
    for ( j = 0; j < N; j++ ) {
        img1[i][j] = input_img.mono[i][j];
    }

    /* blurring fitler */
    for ( i = 0; i < M; i++ )
    for ( j = 0; j < N; j++ ) {
        img2[i][j] = circ_conv2d(img1, M, N, i, j, filter);
    }

    /* Add noise */
    srandom2(1);
    
    for ( i = 0; i < M; i++ )
    for ( j = 0; j < N; j++ ) {
        y[i][j] = img2[i][j] + 4 * normal();
    }

    get_TIFF ( &output_img, M, N, 'g' );
    
    for ( i = 0; i < M; i++ )
    for ( j = 0; j < N; j++ ) {
        pixel = (int32_t)y[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;        
        output_img.mono[i][j] = pixel;
        y[i][j] = (double)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/Y.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file Y.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

/*******************************************************************************************/
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

    for (iter = 0; iter < 20; iter++)
    {
        printf("start iteration %02d...", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            pixel_index_u = i;
            pixel_index_v = j;
            tmp_v = x[i][j];
            theta_1 = circ_conv2d(e, M, N, i, j, filter);
            theta_1 = (-1.0) * theta_1 / (sigma_w * sigma_w);

            low = find_clique_min(x, M, N, i, j);
            high = find_clique_max(x, M, N, i, j);
            if(high < tmp_v - theta_1 / theta_2) high = tmp_v - theta_1 / theta_2;
            if(low > tmp_v - theta_1 / theta_2) low = tmp_v - theta_1 / theta_2;

            tmp_solve = solve(root_function, x, low, high, 1e-7, &err_code);
            if(err_code == 0){
                x[i][j] = tmp_solve;
            }
            else {
                printf("err_code = %0d\n", err_code);
            }
            update_error(e, M, N, i, j, x[i][j] - tmp_v, filter);
        }
        cost[iter] = cost_function(y, x, M, N, sigma_w, sigma_x, filter, prediction);
        printf("cost = %.15f\n", cost[iter]);
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = (int32_t)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate1.tif", "wb" ) ) == NULL ) {
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
    sigma_x = 6.259*(1/5.0);  
  
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


    for (iter = 0; iter < 20; iter++)
    {
        printf("start iteration %02d...", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            pixel_index_u = i;
            pixel_index_v = j;
            tmp_v = x[i][j];
            theta_1 = circ_conv2d(e, M, N, i, j, filter);
            theta_1 = (-1.0) * theta_1 / (sigma_w * sigma_w);

            low = find_clique_min(x, M, N, i, j);
            high = find_clique_max(x, M, N, i, j);
            if(high < tmp_v - theta_1 / theta_2) high = tmp_v - theta_1 / theta_2;
            if(low > tmp_v - theta_1 / theta_2) low = tmp_v - theta_1 / theta_2;

            tmp_solve = solve(root_function, x, low, high, 1e-7, &err_code);
            if(err_code == 0){
                x[i][j] = tmp_solve;
            }
            else {
                printf("err_code = %0d\n", err_code);
            }
            update_error(e, M, N, i, j, x[i][j] - tmp_v, filter);
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
    if ( ( fp = fopen ( "../output/MAPestimate2.tif", "wb" ) ) == NULL ) {
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
    sigma_x = 6.259*(5.0);
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

    for (iter = 0; iter < 20; iter++)
    {
        printf("start iteration %02d...", iter+1);
        for ( i = 0; i < M; i++)
        for ( j = 0; j < N; j++){
            pixel_index_u = i;
            pixel_index_v = j;
            tmp_v = x[i][j];
            theta_1 = circ_conv2d(e, M, N, i, j, filter);
            theta_1 = (-1.0) * theta_1 / (sigma_w * sigma_w);

            low = find_clique_min(x, M, N, i, j);
            high = find_clique_max(x, M, N, i, j);
            if(high < tmp_v - theta_1 / theta_2) high = tmp_v - theta_1 / theta_2;
            if(low > tmp_v - theta_1 / theta_2) low = tmp_v - theta_1 / theta_2;

            tmp_solve = solve(root_function, x, low, high, 1e-7, &err_code);
            if(err_code == 0){
                x[i][j] = tmp_solve;
            }
            else {
                printf("err_code = %0d\n", err_code);
            }
            update_error(e, M, N, i, j, x[i][j] - tmp_v, filter);
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
    if ( ( fp = fopen ( "../output/MAPestimate3.tif", "wb" ) ) == NULL ) {
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

double root_function(double u, void * pblock)
{
    int i, j;
    double val = 0.0;
    double sign;
    double ** img;
    img = (double **)pblock;

    for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++) {
        if (u - img[(pixel_index_u+i-1+M)%M][(pixel_index_v+j-1+N)%N] >= 0) sign = 1.0;
        else sign = -1.0;
        val += pow(fabs(u - img[(pixel_index_u+i-1+M)%M][(pixel_index_v+j-1+N)%N]), p-1) * prediction[i][j] * sign;
    }

    val = val / pow(sigma_x, p);

    val += theta_1 + theta_2 * (u - img[pixel_index_u][pixel_index_v]);

    return val;
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

double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double h[5][5], double g[3][3])
{
    int i, j, m, n;
    double tmp_err;
    double tmp_val1 = 0.0;
    double tmp_val2 = 0.0;
    double cost;
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++){
        tmp_err = y[i][j] - circ_conv2d(x, M, N, i, j, h);
        tmp_val1 += pow(tmp_err, 2.0) ;

        for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++) {
            tmp_val2 += pow(fabs(x[i][j] - x[(i+m-1+M)%M][(j+n-1+N)%N]), p) * g[m][n];
        }
    }
    cost = tmp_val1 / (2.0 * sigma_w * sigma_w) + tmp_val2 / (p * pow(sigma_x, p));
    return cost;
}


double find_clique_max(double **img, int M, int N, int u, int v)
{
    int i, j;
    double max = img[u][v];
    for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++){
        if (img[(u+i-1+M)%M][(v+j-1+N)%N] > max) {
            max = img[(u+i-1+M)%M][(v+j-1+N)%N];
        }
    }
    return max;
}

double find_clique_min(double **img, int M, int N, int u, int v)
{
    int i, j;
    double min = img[u][v];
    for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++){
        if (img[(u+i-1+M)%M][(v+j-1+N)%N] < min) {
            min = img[(u+i-1+M)%M][(v+j-1+N)%N];
        }
    }
    return min;
}
