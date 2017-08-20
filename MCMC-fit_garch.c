
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
    FILE * input;
    FILE *para_output;
    const int G=1;
    const int N =50000;
    const int T =10000;
    const int n=N/5;

    const int k1=50000;
    float y [T] = {};
    float x1[T]={};
    float x_left[T]={};
    float x_right[T]={};
    float x11,x2;
    float phi  ;
// store the variances for the latent disturbances
    float sigma;
    float tau;
    float psi;
    float rho;
    float mu;
// some temp variables
    float a ,b ,a1, b1,a2,b2,a3,b3,a4,b4,a11,a1_left,a1_right,b11,b1_left,c,c1,c11,c1_left,c1_right,w1,w2,rnd1, latent, bounder,left,right;
    float u1, u2,u3,d;
    float aa,bb,cc,dd;
    float z[T]={};
// remember the states for each loop of Giblls
//    float z1[T]={};
// some temp staff
    float alpha,alpha_t, phi_t, alpha_star, u, rnd , phi_temp, chi_temp,temp1,temp;
    int j,k,h,i,g,burnin,m,count_phi,sign;

    float norm_temp, norm1_temp;
    float unif_temp;
    float phi_para, sigma_para,rho_para,mu_para;
    float mean;
    float mu1, sig,p;
    float mm, alpha1, chris_temp;

    float alpha_mu, sig_mu, alpha_phi, sig_phi, alpha_psi, sig_psi, sig_tau;
    int alpha_tau;
    float alpha0,beta1;


// initialized the random generation seed
    long a_time, b_time;
    time_t t1;

    (void) time(&t1);
    a_time=(long) t1;
    b_time=-(a_time+1);

    /*--------------the code of the body----------------*/

// the initial values for the start of the Gibbs
    phi=0.7;
    sigma=1.5; //rooted
    rho=-0.2;
    mu=0.8;
    mm=0.5;
    count_phi=0;

// the known para and to be estimated, we also use stand variance 1
    phi_para=0.5;
    rho_para=-0.5;
    sigma_para=1.0;//rooted
    mu_para=0.5;
    for (g=1; g<=G; g++) // for analyzing data we set G=1
    {
        cout<<"g="<< g<<"\n";
        para_output = fopen ( "C-mcmc-time-series.txt", "w");

  //-------------------------------------------------------------
        input = fopen ("C-simulated-returns.txt", "r");
        h=0;// the time of the repeating sampling//
        while (h<T)
        {
            fscanf (input, "%f", y+h);
            if (y[j+h]==0.0)
            {
                y[j+h]=0.00001;
            }
            h ++;
        }
        fclose (input);

// the initial value for the state x_0
        alpha0=0.01;
        alpha1=0.02;
        beta1=0.02;
        z[0]=y[0]*y[0];
        for ( j=1; j<T;  j++)
        {
            z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1];
        }


        for ( k=1;k<=N; k++)
        {
            left=0.0;
            right=30.0;
            for (j=1;j< T; j++)
            {
                u1=ran1(&b_time);
                b1=u1/sqrt( alpha0 + alpha1*y[j-1]*y[j-1] +beta1*z[j-1]);
                c1_right =1.0/b1/b1 -alpha1*y[j-1]*y[j-1] -beta1*z[j-1];
                u1=ran1(&b_time);
                b1=u1*exp(  -y[j]*y[j]/2.0/( alpha0 + alpha1*y[j-1]*y[j-1] + beta1*z[j-1]));
                c1_left=-y[j]*y[j]/2.0/log(b1) -alpha1*y[j-1]*y[j-1]- beta1*z[j-1];
                left=max(left,c1_left);
                right=min(c1_right,right);
            }
            u=ran1(&b_time);
            alpha0=left + (right-left)*u;
            left=0.0;
            right=30.0;
            for (j=1;j< T; j++)
            {
                u1=ran1(&b_time);
                b1=u1/sqrt( alpha0 + alpha1*y[j-1]*y[j-1]  +beta1*z[j-1]  );
                c1_right =(  1.0/b1/b1 -alpha0 -beta1*z[j-1])/y[j-1]/y[j-1];
                u1=ran1(&b_time);
                b1=u1*exp(  -y[j]*y[j]/2.0/(alpha0+alpha1*y[j-1]*y[j-1]+beta1*z[j-1] )  );
                c1_left=  ( -y[j]*y[j]/2.0/log(b1) -alpha0 -beta1*z[j-1])/y[j-1]/y[j-1];
                left=max(left,c1_left);
                right=min(c1_right,right);
            }
            right=min(1.0-beta1,right);
            u=ran1(&b_time);
            alpha1=left + (right-left)*u;

            left=0.0;
            right=30.0;
            for (j=1; j< T; j++)
            {
                u1=ran1(&b_time);
                b1=u1/sqrt( alpha0+alpha1*y[j-1]*y[j-1]  +beta1*z[j-1] );
                c1_right =( 1.0/b1/b1 -alpha0 -alpha1*y[j-1]*y[j-1])/z[j-1];
                u1=ran1(&b_time);
                b1=u1*exp(  -y[j]*y[j]/2.0/(alpha0+alpha1*y[j-1]*y[j-1]+beta1*z[j-1] )  );
                c1_left=  (  -y[j]*y[j]/2.0/log(b1) -alpha0 -alpha1*y[j-1]*y[j-1])/z[j-1];
                left=max(left,c1_left);
                right=min(c1_right,right);
            }
            right=min(1.0-alpha1,right);
            u=ran1(&b_time);
            beta1=left + (right-left)*u;
            z[0]=y[0]*y[0];
            for ( j=1; j< T;  j++)
            {
                z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1];
            }
            //------------------------------------------------------------------
            fprintf (para_output, " %d %f %f  %f\n",k, alpha0,alpha1, beta1);
            cout <<k<<"\n";
        }
        fclose (para_output);
    }  /// the end of loop g
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
