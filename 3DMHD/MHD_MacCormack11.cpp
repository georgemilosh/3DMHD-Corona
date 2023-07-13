//This mex file simulates inviscid fluid using MacCormack's method
//with MacCormack artificial viscosity
//Simulation involves 2 dimensional grid. 
//There are 6 MHD quantities - mass density, momentum density (x,y)
// ,pressure and mangetic field (x,y). 

//This script can be used for relatively low beta, produces some spurious
// accelerations for large magnetic fields.
// Try this : [u t] = MHD_MacCormack8(1000,10,1/10,5/3,u0);
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

void Forw_Pred(double, int, int, double, double, double***, double***, 
        double**, double**, double**, double**);
void Back_Corr(double, int, int, double, double, double***, double***, 
        double**, double**, double**, double**);
void Prim_Calc(double, int, int, double***, double**, double**);
void Bound_Cond(int, int, double**, int);
void Bound_Cond_w(int, int, double**, int);

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[]) 
{
    
    double ***U; //This is the discretization for U = U(x,y,i) variable
    double ***U_int; //This is the discretization for U = U(x,y,i) variable
    double **u, **v, **artx, **arty;
    double **Ub, **Ub_in;
    
    double *u_out; //This is output 1 dimensional version of U(x,y,i,t)
    double *u0, *ub, *s; //This is the discretization for U = U(x,y,i,0) input
    double *t_out; //This is time output 
    int T, T_num;       //Number of t iterations. Number of outputs in time
    int X,Y;     //Size of X and Y dimensions
    double Delta;   //Parameter: dt/dx = dt/dy
    double gamma;   //gamma - the ratio of specific heat
    int i,j,k,n,a;   //loop variables
    double eps = 20.0; //artificial viscosity constant ~ 1
    int brk, bdim;
    
    //Input Data: (inputs to the mex function in MATLAB)
    
    T_num = (mwSignedIndex)*mxGetPr(prhs[0]);
    Delta = *mxGetPr(prhs[1]);
    gamma = *mxGetPr(prhs[2]);
    X  = (mwSignedIndex)(mxGetDimensions(prhs[3]))[1];
    Y  = (mwSignedIndex)(mxGetDimensions(prhs[3]))[0];
    u0 = mxGetPr(prhs[3]);
    ub = mxGetPr(prhs[4]);
    bdim = (mwSignedIndex)(mxGetDimensions(prhs[4]))[1];
    s  = mxGetPr(prhs[5]);
    T  = (mwSignedIndex)(mxGetDimensions(prhs[5]))[1];
    //printf("X = %d, Y = %d\n", X, Y);
    
    
    //allocate space for U = U(x,y,i,t) and u0
    U = (double***) mxMalloc(6*sizeof(double**));
    U_int = (double***) mxMalloc(6*sizeof(double**));
    Ub = (double**) mxMalloc(6*sizeof(double**));
    for (i = 0; i < 6; i++)
    {
        U[i] = (double**) mxMalloc(X*sizeof(double*));
        U_int[i] = (double**) mxMalloc(X*sizeof(double*));
        for (j = 0; j < X; j++)
        {
            U[i][j] = (double*) mxMalloc(Y*sizeof(double));
            U_int[i][j] = (double*) mxMalloc(Y*sizeof(double));
        }
    }
    Ub = (double**) mxMalloc(bdim*sizeof(double*));
    Ub_in = (double**) mxMalloc(bdim*sizeof(double*));
    for (i = 0; i < bdim; i++)
    {
        Ub[i] = (double*) mxMalloc(X*sizeof(double));
        Ub_in[i] = (double*) mxMalloc(X*sizeof(double));
    }
        
    u = (double**) mxMalloc(X*sizeof(double*));
    v = (double**) mxMalloc(X*sizeof(double*));
    artx = (double**) mxMalloc(X*sizeof(double*));
    arty = (double**) mxMalloc(X*sizeof(double*));
    for (j = 0; j < X; j++)
    {
        u[j] = (double*) mxMalloc(Y*sizeof(double));
        v[j] = (double*) mxMalloc(Y*sizeof(double));
        artx[j] = (double*) mxMalloc(Y*sizeof(double));
        arty[j] = (double*) mxMalloc(Y*sizeof(double));
    }
    
    u_out = (double*) mxMalloc((T_num+1)*X*Y*6*sizeof(double));
    t_out = (double*) mxMalloc((T_num+1)*sizeof(double));

    int u_dim[4] = {Y, X, 6, (T_num+1)};
    plhs[0] = mxCreateNumericArray(4, u_dim, mxDOUBLE_CLASS, mxREAL);
    u_out = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, (T_num+1), mxREAL);
    t_out = mxGetPr(plhs[1]);
    
    //initialize perturbations:
    for (k = 0; k < X; k++)
    {
        
        for (i = 0; i < bdim; i++)
        {
            Ub[i][k] = ub[i*X + k];
            Ub_in[i][k] = u0[i*Y*X + k*Y + 0];
        }
        Ub[1][k] = Ub[1][k]*Ub[0][k];
        Ub[2][k] = Ub[2][k]*Ub[0][k];
        Ub_in[1][k] = Ub_in[1][k]*Ub_in[0][k];
        Ub_in[2][k] = Ub_in[2][k]*Ub_in[0][k];
        
        for (j = 0; j < Y; j++)
        {
            U[0][k][j] = u0[0*Y*X + k*Y + j];
            u[k][j] = u0[1*Y*X + k*Y + j];
            U[1][k][j] = u[k][j]*U[0][k][j];
            v[k][j] = u0[2*Y*X + k*Y + j];
            U[2][k][j] = v[k][j]*U[0][k][j];
            U[3][k][j] = u0[3*Y*X + k*Y + j];
            U[4][k][j] = u0[4*Y*X + k*Y + j];
            U[5][k][j] = u0[5*Y*X + k*Y + j];
            
            
            
            for (i = 0; i < 6; i++)
            {
                u_out[i*X*Y + k*Y + j] = u0[i*Y*X + k*Y + j];
                for (a = 1; a <= T_num; a++)
                    u_out[a*X*Y*6 + i*X*Y + k*Y + j] = 0;
            }
            artx[k][j] = 0;
            arty[k][j] = 0;
            
        }
        
        /*for (i = 0; i < bdim; i++)
            for (j = 1; j < (X-1); j++)
                for (k = 0; k < 1; k++)
                    U[i][j][k] = s[0]*Ub[i][j] + (1-s[0])*Ub_in[i][j];
         */
    }
    t_out[0] = 0;
    a = 1;
    n = 1;
    brk = 0;
    while ((n <= T)&&(brk == 0))
    {
        for (i = 0; i < 6; i++)
            for (j = 0; j < X; j++)
                for (k = 0; k < Y; k++)
                {
                    //Checking whether the solution became NaN
                    if (U[i][j][k] != U[i][j][k])
                        brk = 1;
                    else
                        U_int[i][j][k] = U[i][j][k];
                }
         //stops in case of NaN
        if (brk == 1)
            printf("Stopped on %d iteration\n",n);
        //First Step
        Forw_Pred(Delta, X, Y, gamma, eps, U, U_int, u, v, artx, arty);
        //Second Step
        Back_Corr(Delta, X, Y, gamma, eps, U, U_int, u, v, artx, arty);
        //Boundary Condition:
        for (i = 0; i < bdim; i++)
            for (j = 1; j < (X-1); j++)
                for (k = 0; k < 5; k++)
                    U[i][j][k] = s[n-1]*Ub[i][j] + (1-s[n-1])*Ub_in[i][j];
                    //U[i][j][k] = U[i][j][0];
        //Bound_Cond(X, Y, U[0], 1);
       // Bound_Cond(X, Y, U[1], -1);
       // Bound_Cond(X, Y, U[2], -1);
       // Bound_Cond(X, Y, U[3], 1);
       // Bound_Cond(X, Y, U[4], -1);
       // Bound_Cond(X, Y, U[5], -1);
        
        //for (i = 4; i < 6; i++)
         //   for (j = 1; j < (X-1); j++)
         //       for (k = 1; k < 10; k++)
         //          U[i][j][k] = u0[i*Y*X + j*Y + k];
        //Prim_Calc(gamma, X, Y, U, u, v);
        //Generating Output
        if ((n == ((int)(a*T/T_num + 0.5)))||(brk == 1))
        {
            for (k = 0; k < X; k++)
                for (j = 0; j < Y; j++)
                {
                    //u_out[a*X*Y*6 + 0*X*Y + k*Y + j] = U[i][k][j];
                    u_out[a*X*Y*6 + k*Y + j] = U[0][k][j];
                    u_out[a*X*Y*6 + 1*X*Y + k*Y + j] = u[k][j];
                    u_out[a*X*Y*6 + 2*X*Y + k*Y + j] = v[k][j];
                    u_out[a*X*Y*6 + 3*X*Y + k*Y + j] = U[3][k][j];
                    u_out[a*X*Y*6 + 4*X*Y + k*Y + j] = U[4][k][j];
                    u_out[a*X*Y*6 + 5*X*Y + k*Y + j] = U[5][k][j];
                }
            t_out[a] = Delta*n;
            a++;
        }
        n++;
     }
    
    
}
//FORWARD PREDICTOR
void Forw_Pred(double Delta, int X, int Y, double gamma, double eps,
        double ***U, double ***U_int, double **u, double **v, 
        double **artx, double **arty)
{
    //This procedure updates U_int, p, u and v on the basis of all 
    //other inputs. It incorporates FORWARD PREDICTOR
    int j, k;
    //Artificial Viscosity
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++) 
        {
            artx[j][k] = -eps*(abs(u[j][k]) + 
                    sqrt(gamma*U[3][j][k]/U[0][j][k]) + 
                    sqrt((U[4][j][k]*U[4][j][k] + U[5][j][k]*U[5][j][k])
                    /U[0][j][k]))*abs(U[3][j+1][k] - 2*U[3][j][k] + 
                    U[3][j-1][k])/(U[3][j+1][k] + 2*U[3][j][k] + 
                    U[3][j-1][k]);
            arty[j][k] = -eps*(abs(v[j][k]) + 
                    sqrt(gamma*U[3][j][k]/U[0][j][k]) + 
                    sqrt((U[4][j][k]*U[4][j][k] + U[5][j][k]*U[5][j][k])
                    /U[0][j][k]))*abs(U[3][j][k+1] - 2*U[3][j][k] + 
                    U[3][j][k-1])/(U[3][j][k+1] + 2*U[3][j][k] + 
                    U[3][j][k-1]);
        }
    //Step:
    for (j = 1; j < (X-1); j++)
            for (k = 1; k < (Y-1); k++)
            {
                //Mass Density Equation
                U_int[0][j][k] = U[0][j][k] -
                        Delta*( U[1][j+1][k] - U[1][j][k] +
                        U[2][j][k+1] - U[2][j][k] + 
                        artx[j+1][k]*(U[0][j+1][k] - U[0][j][k]) - 
                        artx[j][k]*(U[0][j][k] - U[0][j-1][k]) + 
                        arty[j][k+1]*(U[0][j][k+1] - U[0][j][k]) - 
                        arty[j][k]*(U[0][j][k] - U[0][j][k-1]) );
                //Momentum Density Equation - x
                U_int[1][j][k] = U[1][j][k] -
                        Delta*(U[1][j+1][k]*u[j+1][k] - 
                        U[1][j][k]*u[j][k] + U[3][j+1][k] - U[3][j][k] +
                        (U[5][j+1][k]*U[5][j+1][k] - 
                        U[5][j][k]*U[5][j][k] - U[4][j+1][k]*U[4][j+1][k] +
                        U[4][j][k]*U[4][j][k])/2  +
                        U[1][j][k+1]*v[j][k+1] - U[1][j][k]*v[j][k] -
                        U[4][j][k+1]*U[5][j][k+1] + 
                        U[4][j][k]*U[5][j][k] +
                        artx[j+1][k]*(U[1][j+1][k] - U[1][j][k]) - 
                        artx[j][k]*(U[1][j][k] - U[1][j-1][k]) + 
                        arty[j][k+1]*(U[1][j][k+1] - U[1][j][k]) - 
                        arty[j][k]*(U[1][j][k] - U[1][j][k-1]));
                //Momentum Density Equation - y
                U_int[2][j][k] = U[2][j][k] -
                        Delta*( U[2][j+1][k]*u[j+1][k] - 
                        U[2][j][k]*u[j][k] -
                        U[4][j+1][k]*U[5][j+1][k] + 
                        U[4][j][k]*U[5][j][k]  + U[3][j][k+1] - U[3][j][k] +
                        U[2][j][k+1]*v[j][k+1] - U[2][j][k]*v[j][k] + 
                        (U[4][j][k+1]*U[4][j][k+1] - 
                        U[4][j][k]*U[4][j][k] - U[5][j][k+1]*U[5][j][k+1] +
                        U[5][j][k]*U[5][j][k])/2 +
                        artx[j+1][k]*(U[2][j+1][k] - U[2][j][k]) - 
                        artx[j][k]*(U[2][j][k] - U[2][j-1][k]) + 
                        arty[j][k+1]*(U[2][j][k+1] - U[2][j][k]) - 
                        arty[j][k]*(U[2][j][k] - U[2][j][k-1]));
                //Pressure Density Equation 
                U_int[3][j][k] = U[3][j][k] -
                        Delta*( U[3][j+1][k]*u[j+1][k] - 
                        U[3][j][k]*u[j][k] + U[3][j][k+1]*v[j][k+1] - 
                        U[3][j][k]*v[j][k] + (gamma - 1)*U[3][j][k]*
                        (u[j+1][k] - u[j][k] +v[j][k+1] - v[j][k]) + 
                        artx[j+1][k]*(U[3][j+1][k] - U[3][j][k]) - 
                        artx[j][k]*(U[3][j][k] - U[3][j-1][k]) + 
                        arty[j][k+1]*(U[3][j][k+1] - U[3][j][k]) - 
                        arty[j][k]*(U[3][j][k] - U[3][j][k-1]));
                //Magnetic Field - x
                U_int[4][j][k] = U[4][j][k] -
                        Delta*(  v[j][k+1]*U[4][j][k+1] -
                        v[j][k]*U[4][j][k] - u[j][k+1]*U[5][j][k+1] + 
                        u[j][k]*U[5][j][k]);
                
                //Magnetic Field - y
                U_int[5][j][k] = U[5][j][k] -
                        Delta*(  - v[j+1][k]*U[4][j+1][k] +
                        v[j][k]*U[4][j][k] + u[j+1][k]*U[5][j+1][k] - 
                        u[j][k]*U[5][j][k]);
                
                
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, U_int, u, v);
}
//BACKWARD CORRECTOR
void Back_Corr(double Delta, int X, int Y, double gamma, double eps,
        double ***U, double ***U_int, double **u, double **v, 
        double **artx, double **arty)
{
    //This procedure updates U, p, u and v on the basis of all 
    //other inputs. It incorporates BACKWARD CORRECTOR
    int j, k;
    //Artificial viscosity
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++) 
        {
            artx[j][k] = -eps*(abs(u[j][k]) + 
                    sqrt(gamma*U_int[3][j][k]/U_int[0][j][k]) + 
                    sqrt((U_int[4][j][k]*U_int[4][j][k] + 
                    U_int[5][j][k]*U_int[5][j][k])/U_int[0][j][k]))*
                    abs(U_int[3][j+1][k] - 2*U_int[3][j][k] + 
                    U_int[3][j-1][k])/(U_int[3][j+1][k] + 2*U_int[3][j][k] + 
                    U_int[3][j-1][k]);
            arty[j][k] = -eps*(abs(v[j][k]) + 
                    sqrt(gamma*U_int[3][j][k]/U_int[0][j][k]) + 
                    sqrt((U_int[4][j][k]*U_int[4][j][k] + 
                    U_int[5][j][k]*U_int[5][j][k])/U_int[0][j][k]))*
                    abs(U_int[3][j][k+1] -  2*U_int[3][j][k] + 
                    U_int[3][j][k-1])/(U_int[3][j][k+1] + 
                    2*U_int[3][j][k] + U_int[3][j][k-1]);
        }
    for (j = 1; j < (X-1); j++)
            for (k = 1; k < (Y-1); k++)
            {
                //Mass Density Equation
                U[0][j][k] = (U[0][j][k] + U_int[0][j][k] -
                        Delta*( U_int[1][j][k] - U_int[1][j-1][k] +
                        U_int[2][j][k] - U_int[2][j][k-1] + 
                        artx[j][k]*(U_int[0][j+1][k] -  U_int[0][j][k]) - 
                        artx[j-1][k]*(U_int[0][j][k] - U_int[0][j-1][k]) + 
                        arty[j][k]*(U_int[0][j][k+1] -  U_int[0][j][k]) - 
                        arty[j][k-1]*(U_int[0][j][k] - U_int[0][j][k-1])
                        ))/2;
                //Momentum Density Equation - x
                U[1][j][k] = (U[1][j][k] + U_int[1][j][k] -
                        Delta*( U_int[1][j][k]*u[j][k] - 
                        U_int[1][j-1][k]*u[j-1][k] + U_int[3][j][k] - 
                        U_int[3][j-1][k] +  (U_int[5][j][k]*U_int[5][j][k] - 
                        U_int[5][j-1][k]*U_int[5][j-1][k] - 
                        U_int[4][j][k]*U_int[4][j][k] +
                        U_int[4][j-1][k]*U_int[4][j-1][k])/2 +
                        U_int[1][j][k]*v[j][k] - 
                        U_int[1][j][k-1]*v[j][k-1]  -
                        U_int[4][j][k]*U_int[5][j][k] + 
                        U_int[4][j][k-1]*U_int[5][j][k-1] +
                        artx[j][k]*(U_int[1][j+1][k] -  U_int[1][j][k]) - 
                        artx[j-1][k]*(U_int[1][j][k] - U_int[1][j-1][k]) + 
                        arty[j][k]*(U_int[1][j][k+1] -  U_int[1][j][k]) - 
                        arty[j][k-1]*(U_int[1][j][k] - U_int[1][j][k-1])
                        ))/2;
                //Momentum Density Equation - y
                U[2][j][k] = (U[2][j][k] + U_int[2][j][k] -
                        Delta*( U_int[2][j][k]*u[j][k] - 
                        U_int[2][j-1][k]*u[j-1][k] -
                        U_int[4][j][k]*U_int[5][j][k] + 
                        U_int[4][j-1][k]*U_int[5][j-1][k] + 
                        U_int[3][j][k] - U_int[3][j][k-1] +
                        U_int[2][j][k]*v[j][k] - 
                        U_int[2][j][k-1]*v[j][k-1] + 
                        (U_int[4][j][k]*U_int[4][j][k] - 
                        U_int[4][j][k-1]*U_int[4][j][k-1] - 
                        U_int[5][j][k]*U_int[5][j][k] +
                        U_int[5][j][k-1]*U_int[5][j][k-1])/2 +
                        artx[j][k]*(U_int[2][j+1][k] -  U_int[2][j][k]) - 
                        artx[j-1][k]*(U_int[2][j][k] - U_int[2][j-1][k]) + 
                        arty[j][k]*(U_int[2][j][k+1] -  U_int[2][j][k]) - 
                        arty[j][k-1]*(U_int[2][j][k] - U_int[2][j][k-1])
                        ))/2;
                //Pressure Density Equation
                U[3][j][k] = (U[3][j][k] + U_int[3][j][k] -
                        Delta*( U_int[3][j][k]*u[j][k] -  
                        U_int[3][j-1][k]*u[j-1][k] + 
                        U_int[3][j][k]*v[j][k] - 
                        U_int[3][j][k-1]*v[j][k-1] + 
                        (gamma - 1)*U_int[3][j][k]*(u[j][k] - u[j-1][k] + 
                        v[j][k] - v[j][k-1]) + 
                        artx[j][k]*(U_int[3][j+1][k] -  U_int[3][j][k]) - 
                        artx[j-1][k]*(U_int[3][j][k] - U_int[3][j-1][k]) + 
                        arty[j][k]*(U_int[3][j][k+1] -  U_int[3][j][k]) - 
                        arty[j][k-1]*(U_int[3][j][k] - U_int[3][j][k-1])
                        ))/2;
                //Magnetic Field - x
                U[4][j][k] = (U[4][j][k] + U_int[4][j][k] -
                        Delta*(  v[j][k]*U_int[4][j][k] -
                        v[j][k-1]*U_int[4][j][k-1] - 
                        u[j][k]*U_int[5][j][k] + 
                        u[j][k-1]*U_int[5][j][k-1] ))/2;
                
                //Magnetic Field - y
                U[5][j][k] = (U[5][j][k] + U_int[5][j][k] -
                        Delta*(  - v[j][k]*U_int[4][j][k] +
                        v[j-1][k]*U_int[4][j-1][k] + 
                        u[j][k]*U_int[5][j][k] - 
                        u[j-1][k]*U_int[5][j-1][k] ))/2;
                
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, U, u, v);
}

void Prim_Calc(double gamma, int X, int Y, double ***U, double **u, 
        double **v)
{
    int j, k;
    //Calculating Velocity and pressure
        for (j = 0; j < X; j++)
            for (k = 0; k < Y; k++)
            {
                if (U[0][j][k] <= 0)
                    printf("Density Negative or Zero!\n");
                u[j][k] = U[1][j][k]/U[0][j][k];
                v[j][k] = U[2][j][k]/U[0][j][k];
                
            }
}

void Bound_Cond(int X, int Y, double **U, int sign)
{
    int j, k;
    
    for (j = 0; j < X; j++)
    {
        //U[j][0] = sign*U[j][1];
        U[j][Y-1] = sign*U[j][Y-2];
    }
    for (k = 0; k < Y; k++)
    {
        U[0][k] = sign*U[1][k];
        U[X-1][k] = sign*U[X-2][k];
    }
            
}

void Bound_Cond_w(int X, int Y, double **U, int sign)
{
    int j, k;
    
    for (j = 0; j < X; j++)
    {
        U[j][0] = sign*U[j][1];
        U[j][Y-1] = sign*U[j][Y-2];
    }
    for (k = 0; k < Y; k++)
    {
        U[0][k] = sign*U[1][k];
        U[X-1][k] = sign*U[X-2][k];
    }
            
}