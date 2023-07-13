//This mex file simulates ideal MHD using MacCormack's method
//with MacCormack artificial viscosity
//Simulation involves 3 dimensional grid. The simulation is fully 3D
//There are 8 quantities - mass density, momentum density (x,y,z), pressure,
// magnetic field (x,y,z)
// Try this : [u t] = MC3DMHD1(16000,10,1/320,5/3,u0);


#include <string.h>
#include "math.h"
#include "mex.h"

void Forw_Pred(double ,int ,int ,int,double ,double , double ****, 
        double ****, double ***, double ***, double ***, double ***, 
        double ***, double ***);
void Back_Corr(double ,int ,int , int,  double, double, double ****, 
        double ****, double ***, double ***, double ***, double ***, 
        double ***, double ***);
void Prim_Calc(double, int , int , int , double ****, double ***, 
        double ***, double ***);
void Bound_Cond(int , int , int , double ***, int);

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[]) 
{
    
    double ****U; //This is the discretization for U = U(x,y,z,i) variable
    double ****U_int; //This is the discretization for U = U(x,y,z,i) variable
    double ***u, ***v, ***w, ***artx, ***arty, ***artz;
    
    
    double *u_out; //This is output 1 dimensional version of U(i,x,y,t)
    double *u0; //This is the discretization for U = U(i,x,y,0) input
    double *t_out; //This is time output 
    int T, T_num;       //Number of t iterations. Number of outputs in time
    int X,Y,Z;     //Size of X and Y dimensions
    double Delta;   //Parameter: dt/dx = dt/dy
    double gamma;   //gamma - the ratio of specific heat
    int i,j,k,r,n,a;   //loop variables
    double eps = 0.0; //artificial viscosity constant ~ 1
    int brk;
    
    //Input Data:
    T  = (mwSignedIndex)*mxGetPr(prhs[0]);
    T_num = (mwSignedIndex)*mxGetPr(prhs[1]);
    Delta = *mxGetPr(prhs[2]);
    gamma = *mxGetPr(prhs[3]);
    X  = (mwSignedIndex)(mxGetDimensions(prhs[4]))[1];
    Y  = (mwSignedIndex)(mxGetDimensions(prhs[4]))[0];
    Z  = (mwSignedIndex)(mxGetDimensions(prhs[4]))[2];
    u0 =  mxGetPr(prhs[4]);
    
    //printf("X = %d, Y = %d\n", X, Y);
    
    
    //allocate space for U = U(i,x,y,z) and u0
    U = (double****) mxMalloc(8*sizeof(double***));
    U_int = (double****) mxMalloc(8*sizeof(double***));
    for (i = 0; i < 8; i++)
    {
        U[i] = (double***) mxMalloc(X*sizeof(double**));
        U_int[i] = (double***) mxMalloc(X*sizeof(double**));
        for (j = 0; j < X; j++)
        {
            U[i][j] = (double**) mxMalloc(Y*sizeof(double*));
            U_int[i][j] = (double**) mxMalloc(Y*sizeof(double*));
            for (k = 0; k < Y; k++)
            {
                U[i][j][k] = (double*) mxMalloc(Z*sizeof(double));
                U_int[i][j][k] = (double*) mxMalloc(Z*sizeof(double));
            }
        }
    }
    u = (double***) mxMalloc(X*sizeof(double**));
    v = (double***) mxMalloc(X*sizeof(double**));
    w = (double***) mxMalloc(X*sizeof(double**));
    artx = (double***) mxMalloc(X*sizeof(double**));
    arty = (double***) mxMalloc(X*sizeof(double**));
    artz = (double***) mxMalloc(X*sizeof(double**));
    for (j = 0; j < X; j++)
    {
        u[j] = (double**) mxMalloc(Y*sizeof(double*));
        v[j] = (double**) mxMalloc(Y*sizeof(double*));
        w[j] = (double**) mxMalloc(Y*sizeof(double*));
        artx[j] = (double**) mxMalloc(Y*sizeof(double*));
        arty[j] = (double**) mxMalloc(Y*sizeof(double*));
        artz[j] = (double**) mxMalloc(Y*sizeof(double*));
        for (k = 0; k < Y; k++)
        {
            u[j][k] = (double*) mxMalloc(Z*sizeof(double));
            v[j][k] = (double*) mxMalloc(Z*sizeof(double));
            w[j][k] = (double*) mxMalloc(Z*sizeof(double));
            artx[j][k] = (double*) mxMalloc(Z*sizeof(double));
            arty[j][k] = (double*) mxMalloc(Z*sizeof(double));
            artz[j][k] = (double*) mxMalloc(Z*sizeof(double));
        }
    }
    
    u_out = (double*) mxMalloc(T_num*X*Y*8*sizeof(double));
    t_out = (double*) mxMalloc(T_num*sizeof(double));
  
    int u_dim[5] = {Y, X, Z, 8, T_num};
    plhs[0] = mxCreateNumericArray(5, u_dim, mxDOUBLE_CLASS, mxREAL);
    u_out = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, T_num, mxREAL);
    t_out = mxGetPr(plhs[1]);
    
    //initialize fields:
    for (k = 0; k < X; k++)
        for (j = 0; j < Y; j++)
            for (r = 0; r < Z; r++)
            {
                U[0][k][j][r] = u0[0*X*Y*Z + r*Y*X + k*Y + j];
                u[k][j][r] = u0[1*X*Y*Z + r*Y*X + k*Y + j];
                U[1][k][j][r] = u[k][j][r]*U[0][k][j][r];
                v[k][j][r] = u0[2*X*Y*Z + r*Y*X + k*Y + j];
                U[2][k][j][r] = v[k][j][r]*U[0][k][j][r];
                w[k][j][r] = u0[3*X*Y*Z + r*Y*X + k*Y + j];
                U[3][k][j][r] = w[k][j][r]*U[0][k][j][r];
                U[4][k][j][r]= u0[4*X*Y*Z + r*Y*X + k*Y + j];
                U[5][k][j][r] = u0[5*X*Y*Z + r*Y*X + k*Y + j];
                U[6][k][j][r] = u0[6*X*Y*Z + r*Y*X + k*Y + j];
                U[7][k][j][r] = u0[7*X*Y*Z + r*Y*X + k*Y + j];
            
                for (i = 0; i < 8; i++)
                    for (a = 0; a < T_num; a++)
                        u_out[a*X*Y*Z*8 + i*X*Y*Z + r*Y*X + k*Y + j] = 0;
                
                artx[k][j][r] = 0;
                arty[k][j][r] = 0;
                artz[k][j][r] = 0;
            }
     //Simulation Starts Here
     t_out[0] = 0;
     a = 0;
     n = 1;
     brk = 0;
     while ((n <= T)&&(brk == 0))
     {
         for (i = 0; i < 8; i++)
            for (j = 0; j < X; j++)
                for (k = 0; k < Y; k++)
                    for (r = 0; r < Z; r++)
                    {
                        //Checking whether the solution became NaN
                        if (U[i][j][k][r] != U[i][j][k][r])
                            brk = 1;
                        else
                            U_int[i][j][k][r] = U[i][j][k][r];
                    }
         //stops in case of NaN
         if (brk == 1)
             printf("Stopped on %d iteration\n",n);
        //First Step
        Forw_Pred(Delta, X, Y, Z, gamma, eps, U, U_int, u, v, w, artx, arty,
                artz);
        //Second Step
        //Back_Corr(Delta, X, Y, Z, gamma, eps, U, U_int, u, v, w, artx, arty,
        //        artz);
        //Boundary Condition:
        
        
       
        //for (r = 0; r < 5; r++)
        //    for (i = 0; i < X; i++)
        //        for (j = 0; j < Y; j++)
        //            for (k = 1; k < 5; k++)
        //                U[r][i][j][k] = U[r][i][j][0]; //Hydrodynamical Quantities should be kept constant near the boundary
													   //where Magnetic Field is strong. Otherwise spurious acceleration is produced.
        
        
        //Bound_Cond(X, Y, Z, U[0], 1);
        //Bound_Cond(X, Y, Z, U[1], -1);
        //Bound_Cond(X, Y, Z, U[2], -1);
        //Bound_Cond(X, Y, Z, U[3], -1);
        //Bound_Cond(X, Y, Z, U[4], 1);
        //Bound_Cond(X, Y, Z, U[5], -1);
        //Bound_Cond(X, Y, Z, U[6], -1);
        //Bound_Cond(X, Y, Z, U[7], -1);
        
        //Prim_Calc(gamma, X, Y, Z, U, u, v, w);
        //Generating Output
        if ((n == ((int)((a+1)*T/T_num + 0.5)))||(brk == 1))
        {
            for (k = 0; k < X; k++)
                for (j = 0; j < Y; j++)
                    for (r = 0; r < Z; r++)
                    {
                        u_out[a*X*Y*Z*8 + 0*X*Y*Z + r*Y*X + k*Y + j] = 
                                U_int[0][k][j][r];
                        u_out[a*X*Y*Z*8 + 1*X*Y*Z + r*Y*X + k*Y + j] = 
                                u[k][j][r];
                        u_out[a*X*Y*Z*8 + 2*X*Y*Z + r*Y*X + k*Y + j] = 
                                v[k][j][r];
                        u_out[a*X*Y*Z*8 + 3*X*Y*Z + r*Y*X + k*Y + j] = 
                                w[k][j][r];
                        u_out[a*X*Y*Z*8 + 4*X*Y*Z + r*Y*X + k*Y + j] = 
                                U_int[4][k][j][r];
                        u_out[a*X*Y*Z*8 + 5*X*Y*Z + r*Y*X + k*Y + j] = 
                                U_int[5][k][j][r];
                        u_out[a*X*Y*Z*8 + 6*X*Y*Z + r*Y*X + k*Y + j] = 
                                U_int[6][k][j][r];
                        u_out[a*X*Y*Z*8 + 7*X*Y*Z + r*Y*X + k*Y + j] = 
                                U_int[7][k][j][r];
                    }
            t_out[a] = Delta*n;
            a++;
        }
        n++;
     }
    
    
}
//FORWARD PREDICTOR
void Forw_Pred(double Delta, int X, int Y, int Z, double gamma, double eps,
        double ****U, double ****U_int, double ***u, double ***v, 
        double ***w, double ***artx, double ***arty, double ***artz)
{
    //This procedure updates U_int, p, u and v on the basis of all 
    //other inputs. It incorporates FORWARD PREDICTOR
    int j, k, r;
    //Artificial Viscosity
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++) 
            for (r = 1; r < (Z-1); r++) 
            {
                artx[j][k][r] = -eps*(abs(u[j][k][r]) + 
                    sqrt(gamma*U[4][j][k][r]/U[0][j][k][r]) + 
                    sqrt((U[5][j][k][r]*U[5][j][k][r] + U[6][j][k][r]*
                    U[6][j][k][r] + U[7][j][k][r]*U[7][j][k][r])
                    /U[0][j][k][r]))*abs(U[4][j+1][k][r] - 2*U[4][j][k][r] + 
                    U[4][j-1][k][r])/(U[4][j+1][k][r] + 2*U[4][j][k][r] + 
                    U[4][j-1][k][r]);
                arty[j][k][r] = -eps*(abs(v[j][k][r]) + 
                    sqrt(gamma*U[4][j][k][r]/U[0][j][k][r]) + 
                    sqrt((U[5][j][k][r]*U[5][j][k][r] + U[6][j][k][r]*
                    U[6][j][k][r] + U[7][j][k][r]*U[7][j][k][r])
                    /U[0][j][k][r]))*abs(U[4][j][k+1][r] - 2*U[4][j][k][r] + 
                    U[4][j][k-1][r])/(U[4][j][k+1][r] + 2*U[4][j][k][r] + 
                    U[4][j][k-1][r]);
                artz[j][k][r] = -eps*(abs(w[j][k][r]) + 
                    sqrt(gamma*U[4][j][k][r]/U[0][j][k][r]) + 
                    sqrt((U[5][j][k][r]*U[5][j][k][r] + U[6][j][k][r]*
                    U[6][j][k][r] + U[7][j][k][r]*U[7][j][k][r])
                    /U[0][j][k][r]))*abs(U[4][j][k][r+1] - 2*U[4][j][k][r] + 
                    U[4][j][k][r-1])/(U[4][j][k][r+1] + 2*U[4][j][k][r] + 
                    U[4][j][k][r-1]);
            }
    //Step:
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++)
            for (r = 1; r < (Z-1); r++)
            {
                //Mass Density Equation
                U_int[0][j][k][r] = U[0][j][k][r] -
                        Delta*( U[1][j+1][k][r] - U[1][j][k][r] +
                        U[2][j][k+1][r] - U[2][j][k][r] + 
                        U[3][j][k][r+1] - U[3][j][k][r] + 
                        //artificial viscosity:
                        artx[j+1][k][r]*(U[0][j+1][k][r] - U[0][j][k][r]) - 
                        artx[j][k][r]*(U[0][j][k][r] - U[0][j-1][k][r]) + 
                        arty[j][k+1][r]*(U[0][j][k+1][r] - U[0][j][k][r]) - 
                        arty[j][k][r]*(U[0][j][k][r] - U[0][j][k-1][r]) +
                        artz[j][k][r+1]*(U[0][j][k][r+1] - U[0][j][k][r]) - 
                        artz[j][k][r]*(U[0][j][k][r] - U[0][j][k][r-1]) );
                
                //Momentum Density Equation - x
                U_int[1][j][k][r] = U[1][j][k][r] -
                        Delta*( //x - derivative:
                        U[1][j+1][k][r]*u[j+1][k][r] - 
                        U[1][j][k][r]*u[j][k][r] + U[4][j+1][k][r] - 
                        U[4][j][k][r] + (U[6][j+1][k][r]*U[6][j+1][k][r] - 
                        U[6][j][k][r]*U[6][j][k][r] + 
                        U[7][j+1][k][r]*U[7][j+1][k][r] - 
                        U[7][j][k][r]*U[7][j][k][r] - 
                        U[5][j+1][k][r]*U[5][j+1][k][r] +
                        U[5][j][k][r]*U[5][j][k][r])/2  +
                        //y - derivative:
                        U[1][j][k+1][r]*v[j][k+1][r] - 
                        U[1][j][k][r]*v[j][k][r] -
                        U[5][j][k+1][r]*U[6][j][k+1][r] + 
                        U[5][j][k][r]*U[6][j][k][r] +
                        //z - derivative:
                        U[1][j][k][r+1]*w[j][k][r+1] - 
                        U[1][j][k][r]*w[j][k][r] -
                        U[5][j][k][r+1]*U[7][j][k][r+1] + 
                        U[5][j][k][r]*U[7][j][k][r] +
                        //artificial viscosity:
                        artx[j+1][k][r]*(U[1][j+1][k][r] - U[1][j][k][r]) - 
                        artx[j][k][r]*(U[1][j][k][r] - U[1][j-1][k][r]) + 
                        arty[j][k+1][r]*(U[1][j][k+1][r] - U[1][j][k][r]) - 
                        arty[j][k][r]*(U[1][j][k][r] - U[1][j][k-1][r]) + 
                        artz[j][k][r+1]*(U[1][j][k][r+1] - U[1][j][k][r]) - 
                        artz[j][k][r]*(U[1][j][k][r] - U[1][j][k][r-1]));
                //Momentum Density Equation - y
                U_int[2][j][k][r] = U[2][j][k][r] -
                        Delta*( //x - derivative:
                        U[2][j+1][k][r]*u[j+1][k][r] - 
                        U[2][j][k][r]*u[j][k][r] -
                        U[6][j+1][k][r]*U[5][j+1][k][r] + 
                        U[6][j][k][r]*U[5][j][k][r]  + 
                        //y - derivative:
                        U[4][j][k+1][r] - U[4][j][k][r] +
                        U[2][j][k+1][r]*v[j][k+1][r] - 
                        U[2][j][k][r]*v[j][k][r] + 
                        (U[5][j][k+1][r]*U[5][j][k+1][r] - 
                        U[5][j][k][r]*U[5][j][k][r] + 
                        U[7][j][k+1][r]*U[7][j][k+1][r] - 
                        U[7][j][k][r]*U[7][j][k][r] - 
                        U[6][j][k+1][r]*U[6][j][k+1][r] +
                        U[6][j][k][r]*U[6][j][k][r])/2  +
                        //z - derivative:
                        U[2][j][k][r+1]*w[j][k][r+1] - 
                        U[2][j][k][r]*w[j][k][r] -
                        U[6][j][k][r+1]*U[7][j][k][r+1] + 
                        U[6][j][k][r]*U[7][j][k][r]  + 
                        //artificial viscosity:
                        artx[j+1][k][r]*(U[2][j+1][k][r] - U[2][j][k][r]) - 
                        artx[j][k][r]*(U[2][j][k][r] - U[2][j-1][k][r]) + 
                        arty[j][k+1][r]*(U[2][j][k+1][r] - U[2][j][k][r]) - 
                        arty[j][k][r]*(U[2][j][k][r] - U[2][j][k-1][r]) + 
                        artz[j][k][r+1]*(U[2][j][k][r+1] - U[2][j][k][r]) - 
                        artz[j][k][r]*(U[2][j][k][r] - U[2][j][k][r-1]));
                //Momentum Density Equation - z
                U_int[3][j][k][r] = U[3][j][k][r] -
                        Delta*( //x - derivative:
                        U[3][j+1][k][r]*u[j+1][k][r] - 
                        U[3][j][k][r]*u[j][k][r] -
                        U[7][j+1][k][r]*U[5][j+1][k][r] + 
                        U[7][j][k][r]*U[5][j][k][r] +
                        //y - derivative:
                        U[3][j][k+1][r]*v[j][k+1][r] - 
                        U[3][j][k][r]*v[j][k][r] -
                        U[7][j][k+1][r]*U[6][j][k+1][r] + 
                        U[7][j][k][r]*U[6][j][k][r] +
                        //z - derivative:
                        U[3][j][k][r+1]*w[j][k][r+1] - 
                        U[3][j][k][r]*w[j][k][r] +
                        U[4][j][k][r+1] - U[4][j][k][r] +
                        (U[5][j][k][r+1]*U[5][j][k][r+1] - 
                        U[5][j][k][r]*U[5][j][k][r] + 
                        U[6][j][k][r+1]*U[6][j][k][r+1] - 
                        U[6][j][k][r]*U[6][j][k][r] - 
                        U[7][j][k][r+1]*U[7][j][k][r+1] +
                        U[7][j][k][r]*U[7][j][k][r])/2  +
                        //artificial viscosity:
                        artx[j+1][k][r]*(U[3][j+1][k][r] - U[3][j][k][r]) - 
                        artx[j][k][r]*(U[3][j][k][r] - U[3][j-1][k][r]) + 
                        arty[j][k+1][r]*(U[3][j][k+1][r] - U[3][j][k][r]) - 
                        arty[j][k][r]*(U[3][j][k][r] - U[3][j][k-1][r]) + 
                        artz[j][k][r+1]*(U[3][j][k][r+1] - U[3][j][k][r]) - 
                        artz[j][k][r]*(U[3][j][k][r] - U[3][j][k][r-1]));
                //Pressure Density Equation 
                U_int[4][j][k][r] = U[4][j][k][r] -
                        Delta*( U[4][j+1][k][r]*u[j+1][k][r] - 
                        U[4][j][k][r]*u[j][k][r] + 
                        U[4][j][k+1][r]*v[j][k+1][r] - 
                        U[4][j][k][r]*v[j][k][r] + 
                        U[4][j][k][r+1]*w[j][k][r+1] - 
                        U[4][j][k][r]*w[j][k][r] + (gamma - 1)*U[4][j][k][r]*
                        (u[j+1][k][r] - u[j][k][r] + v[j][k+1][r] - 
                        v[j][k][r] + w[j][k][r+1] - w[j][k][r]) + 
                        //artificial viscosity:
                        artx[j+1][k][r]*(U[4][j+1][k][r] - U[4][j][k][r]) - 
                        artx[j][k][r]*(U[4][j][k][r] - U[4][j-1][k][r]) + 
                        arty[j][k+1][r]*(U[4][j][k+1][r] - U[4][j][k][r]) - 
                        arty[j][k][r]*(U[4][j][k][r] - U[4][j][k-1][r]) + 
                        artz[j][k][r+1]*(U[4][j][k][r+1] - U[4][j][k][r]) - 
                        artz[j][k][r]*(U[4][j][k][r] - U[4][j][k][r-1]));
                //Magnetic Field - x
                U_int[5][j][k][r] = U[5][j][k][r] -
                        Delta*(  //y - derivative:
                        v[j][k+1][r]*U[5][j][k+1][r] -
                        v[j][k][r]*U[5][j][k][r] - 
                        u[j][k+1][r]*U[6][j][k+1][r] + 
                        u[j][k][r]*U[6][j][k][r] +
                        //z - derivative:
                        w[j][k][r+1]*U[5][j][k][r+1] -
                        w[j][k][r]*U[5][j][k][r] - 
                        u[j][k][r+1]*U[7][j][k][r+1] + 
                        u[j][k][r]*U[7][j][k][r] );
                //Magnetic Field - y
                U_int[6][j][k][r] = U[6][j][k][r] -
                        Delta*( //x - derivative:
                        u[j+1][k][r]*U[6][j+1][k][r] - 
                        u[j][k][r]*U[6][j][k][r] - 
                        v[j+1][k][r]*U[5][j+1][k][r] +
                        v[j][k][r]*U[5][j][k][r] +
                        //z - derivative:
                        w[j][k][r+1]*U[6][j][k][r+1] - 
                        w[j][k][r]*U[6][j][k][r] - 
                        v[j][k][r+1]*U[7][j][k][r+1] +
                        v[j][k][r]*U[7][j][k][r] );
                //Magnetic Field - z
                U_int[7][j][k][r] = U[7][j][k][r] -
                        Delta*(  //x - derivative:
                        u[j+1][k][r]*U[7][j+1][k][r] - 
                        u[j][k][r]*U[7][j][k][r] - 
                        w[j+1][k][r]*U[5][j+1][k][r] +
                        w[j][k][r]*U[5][j][k][r] +
                        //y - derivative:
                        v[j][k+1][r]*U[7][j][k+1][r] -
                        v[j][k][r]*U[7][j][k][r] - 
                        w[j][k+1][r]*U[6][j][k+1][r] + 
                        w[j][k][r]*U[6][j][k][r] );
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, Z, U_int, u, v, w);
}
//BACKWARD CORRECTOR
void Back_Corr(double Delta, int X, int Y, int Z, double gamma, double eps,
        double ****U, double ****U_int, double ***u, double ***v, 
        double ***w, double ***artx, double ***arty, double ***artz)
{
    //This procedure updates U, p, u and v on the basis of all 
    //other inputs. It incorporates BACKWARD CORRECTOR
    int j, k, r;
    //Artificial viscosity
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++) 
            for (r = 1; r < (Z-1); r++)
            {
                artx[j][k][r] = -eps*(abs(u[j][k][r]) + 
                    sqrt(gamma*U_int[4][j][k][r]/U_int[0][j][k][r]) + 
                    sqrt((U_int[5][j][k][r]*U_int[5][j][k][r] + 
                    U_int[6][j][k][r]*U_int[6][j][k][r] + U_int[7][j][k][r]*
                    U_int[7][j][k][r])/U_int[0][j][k][r]))*
                    abs(U_int[4][j+1][k][r] - 2*U_int[4][j][k][r] + 
                    U_int[4][j-1][k][r])/(U_int[4][j+1][k][r] + 
                    2*U_int[4][j][k][r] + U_int[4][j-1][k][r]);
                arty[j][k][r] = -eps*(abs(v[j][k][r]) + 
                    sqrt(gamma*U_int[4][j][k][r]/U_int[0][j][k][r]) + 
                    sqrt((U_int[5][j][k][r]*U_int[5][j][k][r] + 
                    U_int[6][j][k][r]*U_int[6][j][k][r] + 
                    U_int[7][j][k][r]*U_int[7][j][k][r])/U_int[0][j][k][r]))*
                    abs(U_int[4][j][k+1][r] -  2*U_int[4][j][k][r] + 
                    U_int[4][j][k-1][r])/(U_int[4][j][k+1][r] + 
                    2*U_int[4][j][k][r] + U_int[4][j][k-1][r]);
                artz[j][k][r] = -eps*(abs(w[j][k][r]) + 
                    sqrt(gamma*U_int[4][j][k][r]/U_int[0][j][k][r]) + 
                    sqrt((U_int[5][j][k][r]*U_int[5][j][k][r] + 
                    U_int[6][j][k][r]*U_int[6][j][k][r] + 
                    U_int[7][j][k][r]*U_int[7][j][k][r])/U_int[0][j][k][r]))*
                    abs(U_int[4][j][k][r+1] -  2*U_int[4][j][k][r] + 
                    U_int[4][j][k][r-1])/(U_int[4][j][k][r+1] + 
                    2*U_int[4][j][k][r] + U_int[4][j][k][r-1]);
            }
    for (j = 1; j < (X-1); j++)
        for (k = 1; k < (Y-1); k++)
            for (r = 1; r < (Z-1); r++)
            {
                //Mass Density Equation
                U[0][j][k][r] = (U[0][j][k][r] + U_int[0][j][k][r] -
                        Delta*( U_int[1][j][k][r] - U_int[1][j-1][k][r] +
                        U_int[2][j][k][r] - U_int[2][j][k-1][r] + 
                        U_int[3][j][k][r] - U_int[3][j][k][r-1] + 
                        //artificial viscosity:
                        artx[j][k][r]*(U_int[0][j+1][k][r] -  
                        U_int[0][j][k][r]) - artx[j-1][k][r]*
                        (U_int[0][j][k][r] - U_int[0][j-1][k][r]) + 
                        arty[j][k][r]*(U_int[0][j][k+1][r] -  
                        U_int[0][j][k][r]) - arty[j][k-1][r]*
                        (U_int[0][j][k][r] - U_int[0][j][k-1][r]) + 
                        artz[j][k][r]*(U_int[0][j][k][r+1] -  
                        U_int[0][j][k][r]) - artz[j][k][r-1]*
                        (U_int[0][j][k][r] - U_int[0][j][k][r-1])
                        ))/2;
                //Momentum Density Equation - x
                U[1][j][k][r] = (U[1][j][k][r] + U_int[1][j][k][r] -
                        Delta*( //x - derivative:
                        U_int[1][j][k][r]*u[j][k][r] - 
                        U_int[1][j-1][k][r]*u[j-1][k][r] + 
                        U_int[4][j][k][r] - U_int[4][j-1][k][r] +  
                        (U_int[6][j][k][r]*U_int[6][j][k][r] - 
                        U_int[6][j-1][k][r]*U_int[6][j-1][k][r] +
                        U_int[7][j][k][r]*U_int[7][j][k][r] - 
                        U_int[7][j-1][k][r]*U_int[7][j-1][k][r] - 
                        U_int[5][j][k][r]*U_int[5][j][k][r] +
                        U_int[5][j-1][k][r]*U_int[5][j-1][k][r])/2 +
                        //y - derivative:
                        U_int[1][j][k][r]*v[j][k][r] - 
                        U_int[1][j][k-1][r]*v[j][k-1][r]  -
                        U_int[5][j][k][r]*U_int[6][j][k][r] + 
                        U_int[5][j][k-1][r]*U_int[6][j][k-1][r] +
                        // z -derivative:
                        U_int[1][j][k][r]*w[j][k][r] - 
                        U_int[1][j][k][r-1]*w[j][k][r-1]  -
                        U_int[5][j][k][r]*U_int[7][j][k][r] + 
                        U_int[5][j][k][r-1]*U_int[7][j][k][r-1] +
                        //artificial viscosity:
                        artx[j][k][r]*(U_int[1][j+1][k][r] -  
                        U_int[1][j][k][r]) - artx[j-1][k][r]*
                        (U_int[1][j][k][r] - U_int[1][j-1][k][r]) + 
                        arty[j][k][r]*(U_int[1][j][k+1][r] -  
                        U_int[1][j][k][r]) - arty[j][k-1][r]*
                        (U_int[1][j][k][r] - U_int[1][j][k-1][r]) + 
                        artz[j][k][r]*(U_int[1][j][k][r+1] -  
                        U_int[1][j][k][r]) - artz[j][k][r-1]*
                        (U_int[1][j][k][r] - U_int[1][j][k][r-1])
                        ))/2;
                //Momentum Density Equation - y
                U[2][j][k][r] = (U[2][j][k][r] + U_int[2][j][k][r] -
                        Delta*( //x - derivative:
                        U_int[2][j][k][r]*u[j][k][r] - 
                        U_int[2][j-1][k][r]*u[j-1][k][r] -
                        U_int[6][j][k][r]*U_int[5][j][k][r] + 
                        U_int[6][j-1][k][r]*U_int[5][j-1][k][r] + 
                        //y - derivative:
                        U_int[4][j][k][r] - U_int[4][j][k-1][r] +
                        U_int[2][j][k][r]*v[j][k][r] - 
                        U_int[2][j][k-1][r]*v[j][k-1][r] + 
                        (U_int[5][j][k][r]*U_int[5][j][k][r] - 
                        U_int[5][j][k-1][r]*U_int[5][j][k-1][r] +
                        U_int[7][j][k][r]*U_int[7][j][k][r] - 
                        U_int[7][j][k-1][r]*U_int[7][j][k-1][r] - 
                        U_int[6][j][k][r]*U_int[6][j][k][r] +
                        U_int[6][j][k-1][r]*U_int[6][j][k-1][r])/2 +
                        //z - derivative:
                        U_int[2][j][k][r]*w[j][k][r] - 
                        U_int[2][j][k][r-1]*w[j][k][r-1] -
                        U_int[6][j][k][r]*U_int[7][j][k][r] + 
                        U_int[6][j][k][r-1]*U_int[7][j][k][r-1] + 
                        //artificial viscosity:
                        artx[j][k][r]*(U_int[2][j+1][k][r] -  
                        U_int[2][j][k][r]) - artx[j-1][k][r]*
                        (U_int[2][j][k][r] - U_int[2][j-1][k][r]) + 
                        arty[j][k][r]*(U_int[2][j][k+1][r] -  
                        U_int[2][j][k][r]) - arty[j][k-1][r]*
                        (U_int[2][j][k][r] - U_int[2][j][k-1][r]) + 
                        artz[j][k][r]*(U_int[2][j][k][r+1] -  
                        U_int[2][j][k][r]) - artz[j][k][r-1]*
                        (U_int[2][j][k][r] - U_int[2][j][k][r-1])
                        ))/2;
                //Momentum Density Equation - z
                U[3][j][k][r] = (U[3][j][k][r] + U_int[3][j][k][r] -
                        Delta*( //x - derivative:
                        U_int[3][j][k][r]*u[j][k][r] - 
                        U_int[3][j-1][k][r]*u[j-1][k][r] -
                        U_int[7][j][k][r]*U_int[5][j][k][r] + 
                        U_int[7][j-1][k][r]*U_int[5][j-1][k][r] + 
                        //y - derivative:
                        U_int[3][j][k][r]*v[j][k][r] - 
                        U_int[3][j][k-1][r]*v[j][k-1][r] -
                        U_int[7][j][k][r]*U_int[6][j][k][r] + 
                        U_int[7][j][k-1][r]*U_int[6][j][k-1][r] + 
                        //z - derivative:
                        U_int[4][j][k][r] - U_int[4][j][k][r-1] +
                        U_int[3][j][k][r]*w[j][k][r] - 
                        U_int[3][j][k][r-1]*w[j][k][r-1] + 
                        (U_int[5][j][k][r]*U_int[5][j][k][r] - 
                        U_int[5][j][k][r-1]*U_int[5][j][k][r-1] +
                        U_int[6][j][k][r]*U_int[6][j][k][r] - 
                        U_int[6][j][k][r-1]*U_int[6][j][k][r-1] - 
                        U_int[7][j][k][r]*U_int[7][j][k][r] +
                        U_int[7][j][k][r-1]*U_int[7][j][k][r-1])/2 +
                        //artificial viscosity:
                        artx[j][k][r]*(U_int[3][j+1][k][r] -  
                        U_int[3][j][k][r]) - artx[j-1][k][r]*
                        (U_int[3][j][k][r] - U_int[3][j-1][k][r]) + 
                        arty[j][k][r]*(U_int[3][j][k+1][r] -  
                        U_int[3][j][k][r]) - arty[j][k-1][r]*
                        (U_int[3][j][k][r] - U_int[3][j][k-1][r]) + 
                        artz[j][k][r]*(U_int[3][j][k][r+1] -  
                        U_int[3][j][k][r]) - artz[j][k][r-1]*
                        (U_int[3][j][k][r] - U_int[3][j][k][r-1])
                        ))/2;
                //Pressure Density Equation
                 U[4][j][k][r] = (U[4][j][k][r] + U_int[4][j][k][r] -
                        Delta*( U_int[4][j][k][r]*u[j][k][r] -  
                        U_int[4][j-1][k][r]*u[j-1][k][r] + 
                        U_int[4][j][k][r]*v[j][k][r] - 
                        U_int[4][j][k-1][r]*v[j][k-1][r] + 
                        U_int[4][j][k][r]*w[j][k][r] -  
                        U_int[4][j][k][r-1]*w[j][k][r-1] +
                        (gamma - 1)*U_int[4][j][k][r]*(u[j][k][r] - 
                        u[j-1][k][r] + v[j][k][r] - v[j][k-1][r] + 
                        w[j][k][r] - w[j][k][r-1]) + 
                        //artificial viscosity:
                        artx[j][k][r]*(U_int[4][j+1][k][r] -  
                        U_int[4][j][k][r]) - artx[j-1][k][r]*
                        (U_int[4][j][k][r] - U_int[4][j-1][k][r]) + 
                        arty[j][k][r]*(U_int[4][j][k+1][r] -  
                        U_int[4][j][k][r]) - arty[j][k-1][r]*
                        (U_int[4][j][k][r] - U_int[4][j][k-1][r]) + 
                        artz[j][k][r]*(U_int[4][j][k][r+1] -  
                        U_int[4][j][k][r]) - artz[j][k][r-1]*
                        (U_int[4][j][k][r] - U_int[4][j][k][r-1])
                        ))/2;
                //Magnetic Field - x
                U[5][j][k][r] = (U[5][j][k][r] + U_int[5][j][k][r] -
                        Delta*(  //y - derivative:
                        v[j][k][r]*U_int[5][j][k][r] -
                        v[j][k-1][r]*U_int[5][j][k-1][r] - 
                        u[j][k][r]*U_int[6][j][k][r] + 
                        u[j][k-1][r]*U_int[6][j][k-1][r] +
                        //z - derivative:
                        w[j][k][r]*U_int[5][j][k][r] -
                        w[j][k][r-1]*U_int[5][j][k][r-1] - 
                        u[j][k][r]*U_int[7][j][k][r] + 
                        u[j][k][r-1]*U_int[7][j][k][r-1]
                        ))/2;
                //Magnetic Field - y
                U[6][j][k][r] = (U[6][j][k][r] + U_int[6][j][k][r] -
                        Delta*( //x - derivative:
                        u[j][k][r]*U_int[6][j][k][r] - 
                        u[j-1][k][r]*U_int[6][j-1][k][r] -
                        v[j][k][r]*U_int[5][j][k][r] +
                        v[j-1][k][r]*U_int[5][j-1][k][r] + 
                        //z - derivative:
                        w[j][k][r]*U_int[6][j][k][r] - 
                        w[j][k][r-1]*U_int[6][j][k][r-1] -
                        v[j][k][r]*U_int[7][j][k][r] +
                        v[j][k][r-1]*U_int[7][j][k][r-1]
                        ))/2;
                //Magnetic Field - z
                U[7][j][k][r] = (U[7][j][k][r] + U_int[7][j][k][r] -
                        Delta*( //x - derivative:
                        u[j][k][r]*U_int[7][j][k][r] - 
                        u[j-1][k][r]*U_int[7][j-1][k][r] -
                        w[j][k][r]*U_int[5][j][k][r] +
                        w[j-1][k][r]*U_int[5][j-1][k][r] + 
                        //y - derivative:
                        v[j][k][r]*U_int[7][j][k][r] -
                        v[j][k-1][r]*U_int[7][j][k-1][r] - 
                        w[j][k][r]*U_int[6][j][k][r] + 
                        w[j][k-1][r]*U_int[6][j][k-1][r]
                        ))/2;
            }
        //Calculating Velocity and pressure
        Prim_Calc(gamma, X, Y, Z, U, u, v, w);
}

void Prim_Calc(double gamma, int X, int Y, int Z, double ****U, double ***u, 
        double ***v, double ***w)
{
    int j, k, r;
    //Calculating Velocity and pressure
        for (j = 0; j < X; j++)
            for (k = 0; k < Y; k++)
                for (r = 0; r < Z; r++)
                {
                    if (U[0][j][k][r] <= 0)
                    {
                        printf("Density Negative or Zero!\n");
                        printf("x = %d; y = %d; z = %d!\n",j,k,r);
                    }
                    u[j][k][r] = U[1][j][k][r]/U[0][j][k][r];
                    v[j][k][r] = U[2][j][k][r]/U[0][j][k][r];
                    w[j][k][r] = U[3][j][k][r]/U[0][j][k][r];
                }
}

void Bound_Cond(int X, int Y, int Z, double ***U, int sign)
{
    int j, k, r;
    
    for (j = 0; j < X; j++)
        for (k = 0; k < Y; k++)
        {
            U[j][k][Z-1] = sign*U[j][k][Z-2];
        }
    for (k = 0; k < Y; k++)
        for (r = 0; r < Z; r++)
        {
            U[0][k][r] = sign*U[1][k][r];
            U[X-1][k][r] = sign*U[X-2][k][r];
        }
    for (j = 0; j < X; j++)
        for (r = 0; r < Z; r++)
        {
            U[j][0][r] = sign*U[j][1][r];
            U[j][Y-1][r] = sign*U[j][Y-2][r];
        }           
}