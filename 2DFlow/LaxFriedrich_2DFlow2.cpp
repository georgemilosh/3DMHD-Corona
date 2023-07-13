//This mex file simulates inviscid fluid using Lax-Friedrich's method
//Simulation involves 2 dimensional grid. 
//There are 4 hydrodynamic quantities - mass density, momentum density (x,y)
// and energy density (thermal + kinetic)
// Try this : [u t] = LaxFriedrich_2DFlow2(1000,10,1/10,5/3,u0);
//
// use  k = 11; figure; subplot(2,2,1), image(u(:,:,1,k),...
//'CDataMapping','scaled'), axis xy, colormap(jet(256)); colorbar;
//xlabel('x'); ylabel('y'); title(['Density at ' num2str(t(k))]); 
//subplot(2,2,2), image(u(:,:,2,k),'CDataMapping','scaled'); 
//axis xy, colormap(jet(256)); colorbar; xlabel('x'); ylabel('y'); 
//title(['x-momentum at ' num2str(t(k))]); subplot(2,2,3); 
//image(u(:,:,3,k),'CDataMapping','scaled'), axis xy, colormap(jet(256)); 
//colorbar; xlabel('x'); ylabel('y'); 
//title(['y-momentum at ' num2str(t(k))]); 
//subplot(2,2,4), image(u(:,:,4,k),'CDataMapping','scaled'), axis xy;
//colormap(jet(256)); colorbar; xlabel('x'); ylabel('y'); 
//title(['Energy at ' num2str(t(k))]);
// 
#include <string.h>
#include "math.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[]) 
{
    
    double ***U; //This is the discretization for U = U(x,y,i) variable
    double ***U_old; //This is the discretization for U = U(x,y,i) variable
    double **u, **v;
    
    
    double *u_out; //This is output 1 dimensional version of U(x,y,i,t)
    double *u0; //This is the discretization for U = U(x,y,i,0) input
    double *t_out; //This is time output
    int T, T_num;       //Number of t iterations. Number of outputs in time
    int X,Y;     //Size of X and Y dimensions
    double Delta;   //Parameter: dt/dx = dt/dy
    double gamma;   //gamma - the ratio of specific heat
    int i,j,k,n,a;   //loop variables
    
    //Input Data:
    T  = (mwSignedIndex)*mxGetPr(prhs[0]);
    T_num = (mwSignedIndex)*mxGetPr(prhs[1]);
    Delta = *mxGetPr(prhs[2]);
    gamma = *mxGetPr(prhs[3]);
    X  = (mwSignedIndex)(mxGetDimensions(prhs[4]))[1];
    Y  = (mwSignedIndex)(mxGetDimensions(prhs[4]))[0];
    u0 =  mxGetPr(prhs[4]);
    
    //printf("X = %d, Y = %d\n", X, Y);
    
    
    //allocate space for U = U(x,y,i,t) and u0
    U = (double***) mxMalloc(X*sizeof(double**));
    U_old = (double***) mxMalloc(X*sizeof(double**));
    u = (double**) mxMalloc(X*sizeof(double*));
    v = (double**) mxMalloc(X*sizeof(double*));
    for (j = 0; j < X; j++)
    {
        u[j] = (double*) mxMalloc(Y*sizeof(double));
        v[j] = (double*) mxMalloc(Y*sizeof(double));
        U[j] = (double**) mxMalloc(Y*sizeof(double*));
        U_old[j] = (double**) mxMalloc(Y*sizeof(double*));
        for (k = 0; k < Y; k++)
        {
            U[j][k] = (double*) mxMalloc(4*sizeof(double));
            U_old[j][k] = (double*) mxMalloc(4*sizeof(double));
        }
    }
    
    u_out = (double*) mxMalloc((T_num+1)*X*Y*4*sizeof(double));
    t_out = (double*) mxMalloc((T_num+1)*sizeof(double));
   
    int u_dim[4] = {Y, X, 4, (T_num+1)};
    plhs[0] = mxCreateNumericArray(4, u_dim, mxDOUBLE_CLASS, mxREAL);
    u_out = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, (T_num+1), mxREAL);
    t_out = mxGetPr(plhs[1]);
    
    //initialize perturbations:
    
    for (k = 0; k < X; k++)
        for (j = 0; j < Y; j++)
            for (i = 0; i < 4; i++)
            {
                //U_old[k][j][i] = u0[i*Y*X + k*Y + j];
                U[k][j][i] = u0[i*Y*X + k*Y + j];
                u_out[i*X*Y + k*Y + j] = U[k][j][i];
            }
     t_out[0] = 0;
     a = 1;
     
     for (n = 1; n <= T; n++)
     {
        //memcpy(U_old, U, sizeof(U_old));
        for (j = 0; j < X; j++)
            for (k = 0; k < Y; k++)
                for (i = 0; i < 4; i++)
                    U_old[j][k][i] = U[j][k][i];
        //Calculating Velocity and pressure
        for (j = 0; j < X; j++)
            for (k = 0; k < Y; k++)
            {
                u[j][k] = U_old[j][k][1]/U_old[j][k][0];
                v[j][k] = U_old[j][k][2]/U_old[j][k][0];
            }
        //Calculating Density, Momentum and Energy
        for (j = 1; j < (X-1); j++)
            for (k = 1; k < (Y-1); k++)
            {
                //Mass Density Equation
                U[j][k][0] = (U_old[j+1][k][0] + U_old[j][k+1][0] +
                        U_old[j-1][k][0] + U_old[j][k-1][0])/4 -
                        Delta*( U_old[j+1][k][1] - U_old[j-1][k][1] +
                        U_old[j][k+1][2] - U_old[j][k-1][2])/2;
                //Momentum Density Equation - x
                U[j][k][1] = (U_old[j+1][k][1] + U_old[j][k+1][1] +
                        U_old[j-1][k][1] + U_old[j][k-1][1])/4 -
                        Delta*( U_old[j+1][k][1]*u[j+1][k] - 
                        U_old[j-1][k][1]*u[j-1][k] + 
                        U_old[j+1][k][3] - U_old[j-1][k][3] +
                        U_old[j][k+1][1]*v[j][k+1] -  
                        U_old[j][k-1][1]*v[j][k-1])/2;
                //Momentum Density Equation - y
                U[j][k][2] = (U_old[j+1][k][2] + U_old[j][k+1][2] +
                        U_old[j-1][k][2] + U_old[j][k-1][2])/4 -
                        Delta*( U_old[j+1][k][2]*u[j+1][k] - 
                        U_old[j-1][k][2]*u[j-1][k] + 
                        U_old[j][k+1][3] - U_old[j][k-1][3] +
                        U_old[j][k+1][2]*v[j][k+1] -  
                        U_old[j][k-1][2]*v[j][k-1])/2;
                //Energy Density Equation
                U[j][k][3] = (U_old[j+1][k][3] + U_old[j][k+1][3] +
                        U_old[j-1][k][3] + U_old[j][k-1][3])/4 -
                        Delta*( U_old[j+1][k][3]*u[j+1][k] - 
                        U_old[j-1][k][3]*u[j-1][k] +
                        U_old[j][k+1][3]*v[j][k+1] -  
                        U_old[j][k-1][3]*v[j][k-1] + 
                        (gamma-1)*U_old[j][k][3]*(u[j+1][k]-u[j-1][k] +
                        v[j][k+1]-v[j][k-1]))/2;
            }
        //Generating Output
        if (n == ((int)(a*T/T_num + 0.5)))
        {
            for (k = 0; k < X; k++)
                for (j = 0; j < Y; j++)
                    for (i = 0; i < 4; i++)
                        u_out[a*X*Y*4 + i*X*Y + k*Y + j] = U[k][j][i];
            t_out[a] = Delta*n;
            a++;
        }
     }
    
    
}