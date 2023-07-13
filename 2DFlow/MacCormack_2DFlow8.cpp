//This mex file simulates inviscid fluid using MacCormack's method
//with MacCormack artificial viscosity
//Simulation involves 2 dimensional grid. 
//There are 4 hydrodynamic quantities - mass density, momentum density (x,y)
// and energy density (thermal + kinetic)
// Try this : [u t] = MacCormack_2DFlow2(1000,10,1/10,5/3,u0);
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
        double**, double**, double**, double**, double**);
void Back_Corr(double, int, int, double, double, double***, double***, 
        double**, double**, double**, double**, double**);
void Prim_Calc(double, int, int, double***, double**, double**, double**);

void Boundary(int, int, double**);

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[]) 
{
    
    double ***U; //This is the discretization for U = U(x,y,i) variable
    double ***U_int; //This is the discretization for U = U(x,y,i) variable
    double **u, **v, **p, **artx, **arty;
    
    
    double *u_out; //This is output 1 dimensional version of U(x,y,i,t)
    double *u0; //This is the discretization for U = U(x,y,i,0) input
    double *t_out; //This is time output 
    int T, T_num;       //Number of t iterations. Number of outputs in time
    int X,Y;     //Size of X and Y dimensions
    double Delta;   //Parameter: dt/dx = dt/dy
    double gamma;   //gamma - the ratio of specific heat
    int i,j,k,n,a, distx, disty;   //loop variables
    double eps = 1.0; //artificial viscosity constant ~ 1
    
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
    U = (double***) mxMalloc(4*sizeof(double**));
    U_int = (double***) mxMalloc(4*sizeof(double**));
    for (i = 0; i < 4; i++)
    {
        U[i] = (double**) mxMalloc(X*sizeof(double*));
        U_int[i] = (double**) mxMalloc(X*sizeof(double*));
        for (j = 0; j < X; j++)
        {
            U[i][j] = (double*) mxMalloc(Y*sizeof(double));
            U_int[i][j] = (double*) mxMalloc(Y*sizeof(double));
        }
    }
    
    u = (double**) mxMalloc(X*sizeof(double*));
    v = (double**) mxMalloc(X*sizeof(double*));
    p = (double**) mxMalloc(X*sizeof(double*));
    artx = (double**) mxMalloc(X*sizeof(double*));
    arty = (double**) mxMalloc(X*sizeof(double*));
    for (j = 0; j < X; j++)
    {
        u[j] = (double*) mxMalloc(Y*sizeof(double));
        v[j] = (double*) mxMalloc(Y*sizeof(double));
        p[j] = (double*) mxMalloc(Y*sizeof(double));
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
        for (j = 0; j < Y; j++)
        {
            U[0][k][j] = u0[0*Y*X + k*Y + j];
            u[k][j] = u0[1*Y*X + k*Y + j];
            U[1][k][j] = u[k][j]*U[0][k][j];
            v[k][j] = u0[2*Y*X + k*Y + j];
            U[2][k][j] = v[k][j]*U[0][k][j];
            p[k][j] = u0[3*Y*X + k*Y + j];
            U[3][k][j] = p[k][j]/(gamma-1) +  U[0][k][j]*(u[k][j]*u[k][j] +
                    v[k][j]*v[k][j])/2;
            
            for (i = 0; i < 4; i++)
            {
                u_out[i*X*Y + k*Y + j] = u0[i*Y*X + k*Y + j];
            }
            artx[k][j] = 0;
            arty[k][j] = 0;
            u_out[4*X*Y + k*Y + j] = artx[k][j];
            u_out[5*X*Y + k*Y + j] = arty[k][j];
        }
     t_out[0] = 0;
     a = 1;
     n = 1;
     while (n <= T)
     {
        
        for (i = 0; i < 4; i++)
            for (j = 0; j < X; j++)
                for (k = 0; k < Y; k++)
                    U_int[i][j][k] = U[i][j][k];
        //First Step
        Forw_Pred(Delta, X, Y, gamma, eps, U, U_int, p, u, v, artx, arty);
        //Boundary Conditions:
        Boundary(X, Y, u);
        Boundary(X, Y, v);
        Boundary(X, Y, p);
        for (i = 0; i < 4; i++)
            Boundary(X, Y , U_int[i]);
        //Second Step
        Back_Corr(Delta, X, Y, gamma, eps, U, U_int, p, u, v, artx, arty);
        //Boundary Conditions:
        Boundary(X, Y, u);
        Boundary(X, Y, v);
        Boundary(X, Y, p);
        for (i = 0; i < 4; i++)
            Boundary(X, Y , U[i]);
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, U, u, v, p);
        //Generating Output
        if (n == ((int)(a*T/T_num + 0.5)))
        {
            for (k = 0; k < X; k++)
                for (j = 0; j < Y; j++)
                {
                    //u_out[a*X*Y*6 + 0*X*Y + k*Y + j] = U[i][k][j];
                    u_out[a*X*Y*6 + k*Y + j] = U[0][k][j];
                    u_out[a*X*Y*6 + 1*X*Y + k*Y + j] = u[k][j];
                    u_out[a*X*Y*6 + 2*X*Y + k*Y + j] = v[k][j];
                    u_out[a*X*Y*6 + 3*X*Y + k*Y + j] = p[k][j];
                    u_out[a*X*Y*6 + 4*X*Y + k*Y + j] = artx[k][j];
                    u_out[a*X*Y*6 + 5*X*Y + k*Y + j] = arty[k][j];
                }
            t_out[a] = Delta*n;
            a++;
        }
        n++;
     }
    
    
}
//FORWARD PREDICTOR
void Forw_Pred(double Delta, int X, int Y, double gamma, double eps,
        double ***U, double ***U_int, double **p, double **u, double **v, 
        double **artx, double **arty)
{
    //This procedure updates U_int, p, u and v on the basis of all 
    //other inputs. It incorporates FORWARD PREDICTOR
    int j, k;
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++) 
        {
            artx[j][k] = -eps*(abs(u[j][k]) + 
                    sqrt(gamma*p[j][k]/U[0][j][k]))*abs(p[j+1][k] - 
                    2*p[j][k] + p[j-1][k])/(p[j+1][k] + 2*p[j][k] + 
                    p[j-1][k]);
            arty[j][k] = -eps*(abs(v[j][k]) + 
                    sqrt(gamma*p[j][k]/U[0][j][k]))*abs(p[j][k+1] - 
                    2*p[j][k] + p[j][k-1])/(p[j][k+1] + 2*p[j][k] + 
                    p[j][k-1]);
            artx[j][k] = 10;
            arty[j][k] = 10;
                
        }
    //Step:
    for (j = 0; j < (X-1); j++)
            for (k = 0; k < (Y-1); k++)
            {
                //Mass Density Equation
                U_int[0][j][k] = U[0][j][k] -
                        Delta*( U[1][j+1][k] - U[1][j][k] +
                        U[2][j][k+1] - U[2][j][k] /*+ 
                        artx[j+1][k]*(U[0][j+1][k] - U[0][j][k]) - 
                        artx[j][k]*(U[0][j][k] - U[0][j-1][k]) + 
                        arty[j][k+1]*(U[0][j][k+1] - U[0][j][k]) - 
                        arty[j][k]*(U[0][j][k] - U[0][j][k-1]) */);
                //Momentum Density Equation - x
                U_int[1][j][k] = U[1][j][k] -
                        Delta*( U[1][j+1][k]*u[j+1][k] - 
                        U[1][j][k]*u[j][k] + p[j+1][k] - p[j][k] +
                        U[1][j][k+1]*v[j][k+1] - U[1][j][k]*v[j][k] /*+ 
                        artx[j+1][k]*(U[1][j+1][k] - U[1][j][k]) - 
                        artx[j][k]*(U[1][j][k] - U[1][j-1][k]) + 
                        arty[j][k+1]*(U[1][j][k+1] - U[1][j][k]) - 
                        arty[j][k]*(U[1][j][k] - U[1][j][k-1])*/);
                //Momentum Density Equation - y
                U_int[2][j][k] = U[2][j][k] -
                        Delta*( U[2][j+1][k]*u[j+1][k] - 
                        U[2][j][k]*u[j][k] + p[j][k+1] - p[j][k] +
                        U[2][j][k+1]*v[j][k+1] - U[2][j][k]*v[j][k] /*+ 
                        artx[j+1][k]*(U[2][j+1][k] - U[2][j][k]) - 
                        artx[j][k]*(U[2][j][k] - U[2][j-1][k]) + 
                        arty[j][k+1]*(U[2][j][k+1] - U[2][j][k]) - 
                        arty[j][k]*(U[2][j][k] - U[2][j][k-1])*/);
                //Energy Density Equation
                U_int[3][j][k] = U[3][j][k] -
                        Delta*( (U[3][j+1][k] + p[j+1][k])*u[j+1][k] - 
                        (U[3][j][k] + p[j][k])*u[j][k]   +
                        (U[3][j][k+1] + p[j][k+1])*v[j][k+1] -  
                        (U[3][j][k] + p[j][k])*v[j][k] /*+ 
                        artx[j+1][k]*(U[3][j+1][k] - U[3][j][k]) - 
                        artx[j][k]*(U[3][j][k] - U[3][j-1][k]) + 
                        arty[j][k+1]*(U[3][j][k+1] - U[3][j][k]) - 
                        arty[j][k]*(U[3][j][k] - U[3][j][k-1])*/);
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, U_int, u, v, p);
}
//BACKWARD CORRECTOR
void Back_Corr(double Delta, int X, int Y, double gamma, double eps,
        double ***U, double ***U_int, double **p, double **u, double **v, 
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
                    sqrt(gamma*p[j][k]/U_int[0][j][k]))*abs(p[j+1][k] - 
                    2*p[j][k] + p[j-1][k])/(p[j+1][k] + 2*p[j][k] + 
                    p[j-1][k]);
            arty[j][k] = -eps*(abs(v[j][k]) + 
                    sqrt(gamma*p[j][k]/U_int[0][j][k]))*abs(p[j][k+1] - 
                    2*p[j][k] + p[j][k-1])/(p[j][k+1] + 
                    2*p[j][k] + p[j][k-1]);
        }
    for (j = 1; j < X; j++)
            for (k = 1; k < Y; k++)
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
                        U_int[1][j-1][k]*u[j-1][k] + p[j][k] - p[j-1][k] +
                        U_int[1][j][k]*v[j][k] -  
                        U_int[1][j][k-1]*v[j][k-1] + 
                        artx[j][k]*(U_int[1][j+1][k] -  U_int[1][j][k]) - 
                        artx[j-1][k]*(U_int[1][j][k] - U_int[1][j-1][k]) + 
                        arty[j][k]*(U_int[1][j][k+1] -  U_int[1][j][k]) - 
                        arty[j][k-1]*(U_int[1][j][k] - U_int[1][j][k-1])
                        ))/2;
                //Momentum Density Equation - y
                U[2][j][k] = (U[2][j][k] + U_int[2][j][k] -
                        Delta*( U_int[2][j][k]*u[j][k] - 
                        U_int[2][j-1][k]*u[j-1][k] + p[j][k] - p[j][k-1] +
                        U_int[2][j][k]*v[j][k] -  
                        U_int[2][j][k-1]*v[j][k-1] + 
                        artx[j][k]*(U_int[2][j+1][k] -  U_int[2][j][k]) - 
                        artx[j-1][k]*(U_int[2][j][k] - U_int[2][j-1][k]) + 
                        arty[j][k]*(U_int[2][j][k+1] -  U_int[2][j][k]) - 
                        arty[j][k-1]*(U_int[2][j][k] - U_int[2][j][k-1])
                        ))/2;
                //Energy Density Equation
                U[3][j][k] = (U[3][j][k] + U_int[3][j][k] -
                        Delta*( (U_int[3][j][k] + p[j][k])*u[j][k] - 
                        (U_int[3][j-1][k] + p[j-1][k])*u[j-1][k]   +
                        (U_int[3][j][k] + p[j][k])*v[j][k] -  
                        (U_int[3][j][k-1] + p[j][k-1])*v[j][k-1] + 
                        artx[j][k]*(U_int[3][j+1][k] -  U_int[3][j][k]) - 
                        artx[j-1][k]*(U_int[3][j][k] - U_int[3][j-1][k]) + 
                        arty[j][k]*(U_int[3][j][k+1] -  U_int[3][j][k]) - 
                        arty[j][k-1]*(U_int[3][j][k] - U_int[3][j][k-1])
                        ))/2;
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, U, u, v, p);
}

void Prim_Calc(double gamma, int X, int Y, double ***U, double **u, 
        double **v, double **p)
{
    int j, k;
    //Calculating Velocity and pressure
        for (j = 0; j < X; j++)
            for (k = 0; k < Y; k++)
            {
                u[j][k] = U[1][j][k]/U[0][j][k];
                v[j][k] = U[2][j][k]/U[0][j][k];
                p[j][k] = (gamma-1)*(U[3][j][k] - 
                        U[0][j][k]*(u[j][k]*u[j][k] + 
                        v[j][k]*v[j][k])/2);
            }
}
//Imposes Boundary Conditions:
void Boundary(int X, int Y , double **var)
{
    int j, k;
    for (j = 1; j < (X-1); j++)
        var[j][0] = 2*var[j][1] - var[j][2];
        //var[j][0] = (5*var[j][1] - 4*var[j][2] + var[j][3])/2;
    for (k = 0; k < (Y-1); k++)
        var[X-1][k] = 2*var[X-2][k] - var[X-3][k];
        //var[X-1][k] = (5*var[X-2][k] - 4*var[X-3][k] + var[X-4][k])/2;
    for (j = (X-1); j > 0; j--)
        var[j][Y-1] = 2*var[j][Y-2] - var[j][Y-3];
        //var[j][Y-1] = (5*var[j][Y-2] - 4*var[j][Y-3] + var[j][Y-4])/2;
    for (k = 0; k < Y; k++)
        var[0][k] = 2*var[1][k] - var[2][k];
        //var[0][k] = (5*var[1][k] - 4*var[2][k] + var[3][k])/2;
}