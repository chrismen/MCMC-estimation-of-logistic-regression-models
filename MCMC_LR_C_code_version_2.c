/* -----------------------------------------------------------------------

        This c file do slice samplig for ASV mdoels

     the slice sampler within Gibbs and acceptance-rejection method
       with normal nosies

                            January, 2009
Aug 5, 2009
--------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <time.h>
#include <math_recipes.h>

double min(double x, double y);
double max(double x, double y);

//----------------------------------------------
int main ()
{
    FILE * outputfile;
    FILE * input;
    FILE *para_output;
    FILE *para_result;
    FILE *y_output;
    const int G=1;
    const int N =5000;
    const int rows=10000;
    const int cols=3;
// the burn in number of the sample
    const int n=N/5;
// some temp variables
    float temp1, temp_current,temp_new, u,likelihood;
    float sigma2;
    float x[rows][cols]={};  //
    float y[rows+1]={};  // the response vector
    float beta_current[cols+1]={};
    float beta_new[cols+1]={};
    float yy[cols]={};
    int j,k,h,g;
    long a_time, b_time;
    time_t t1;
    (void) time(&t1);
    a_time=(long) t1;
    b_time=-(a_time+1);
    sigma2=0.05;
    para_output = fopen ( "C-MCMC-time-seris-normal-version-2.txt", "w");
//-------------------------------------------------------------
    input = fopen ("C-simulated-logist-data-x.txt", "r");
    for ( h=0; h<rows;  h++)
    {
        for ( j=0; j<(cols-1);  j++)
        {
            fscanf (input, "%f ", &temp1);
            x[h][j]=temp1;
        }
        fscanf (input, "%f \n", &temp1);
        x[h][(cols-1)]=temp1;
    }
    fclose (input);
    input = fopen ("C-simulated-logist-data-y.txt", "r");
    for ( h=0; h<rows;  h++)
    {
        fscanf (input, "%f \n", &temp1);
        y[h]=temp1;
    }
    fclose (input);

    for ( g=1;g<=N; g++)
    {
        cout<< g<<"\n ";
        for ( k=0;k<(cols+1); k++)
        {
            temp1=beta_current[k] +sigma2*gasdev(&b_time);
            beta_new[k]=temp1;
            likelihood=0.0;
            for (h=0; h<rows; h++)
            {
                temp_new=beta_new[0];
                for (j=1; j<cols+1; j++)
                {
                    temp_new=temp_new+beta_new[j]*x[h][j-1];
                }
                temp_new=exp(temp_new);

                temp_current=beta_current[0];
                for (j=1; j<cols+1; j++)
                {
                    temp_current = temp_current + beta_current[j]*x[h][j-1];
                }
                temp_current=exp(temp_current);

                likelihood=likelihood + log( temp_new /(1.0 + temp_new))*y[h];
                likelihood=likelihood + log(1.0- temp_new /(1.0 + temp_new))*(1.0 - y[h]);
                likelihood=likelihood - log(temp_current /(1.0 + temp_current))*y[h];
                likelihood=likelihood - log(1.0- temp_current /(1.0 + temp_current))*(1.0-y[h]);
            }
            u=ran1(&b_time);
            if ( u<min(1.0, exp(likelihood)))
            {
                beta_current[k]=beta_new[k];
            }
        } // end of k loop
        if ( g>50)
        {
            fprintf (para_output, "%d  %f ",g, likelihood);
            for ( j=0; j<cols+1;  j++)
            {
                fprintf (para_output, " %f", beta_current[j]);
            }
            fprintf (para_output, " \n");
        }
    } /// the end of the loop k with the length N for gererating samples
    fclose (para_output);
    return 0;
}

///------------------------------------------------------------------------------
double min(double x, double y)
{
    if (x<y) return x;
    else
        return y;

}
///----------------------------------------
double max(double x, double y)
{
    if (x<y) return y;
    else
        return x;

}
