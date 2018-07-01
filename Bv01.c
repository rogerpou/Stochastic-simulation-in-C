/*
------------------------------------------------------
Modificaciones respecto versión anterior:
1-La función box_muller se llama dentro del bucle esto
 ahorra tiempo.
2-Cambio de la función potencial
------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tiempo (float)(10)//Tiempo final
#define dt (float)(0.01)//Paso en tiempo
#define Temperatura (float)(0.3)// Constante Difusión
#define dT (float)(0.01)//Paso en tiempo
#define trayec (int)(10000)
#define M   (float)(1.0)
#define Eta (float)(1.0)
#define K   (float)(2*3.141592)
#define A   (float)(0.5)
#define omega   (float)(0.01)
#define PasosT (int) (0)  //pasos termalización



void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float box_muller(float m, float s);

unsigned char ind_ran,ig1,ig2,ig3;
unsigned int Wheel[256],ir1;
FILE *D1;
FILE *D2;
int main()
{
    ini_ran(123456789);
    D1=fopen("T-E_cin-D.txt","w");
    D2=fopen("2.txt","w");
    float x,w,v,D[2],epsilon,**zz,coseno[(int)(tiempo/dt)];
    register float temp,temp1,temp2;temp=K/M;temp1=Eta/M,temp2=0.5*dt;
    float vmedia[trayec];


	zz = (float **)malloc(trayec*sizeof(float*));
	for (int i=0;i<trayec;i++)
        {
            zz[i]=(float*)malloc((int)(tiempo/dt)*sizeof(float));
            vmedia[i] = 0;
        }

    for(float T=0;T<=Temperatura;T+=dT)/*Para varias temperaturas*/
    {
        D[0]=D[1]=0;epsilon=pow((dt*2*T),0.5);
        for(int k=0;k<trayec;k++)
        {
            x=w=v=0;
            for(int t=0;t<(int)(tiempo/dt);t++)
            {
                register float g1,g2;
                if(T==0&&k<trayec&&t<(int)(tiempo/dt))zz[k][t]=box_muller(0,1);
                if(T==0&&k==0&&t<(int)(tiempo/dt))coseno[t]=cos(omega*t*dt);
                g1=cos(K*(x+epsilon*zz[k][t]))+0.5*cos(2*K*(x+epsilon*zz[k][t]))+A*sin(w);  /*Calculo las funciones g11,g12,... para el rk estocástico*/
                g2=cos(K*(x+g1*dt))+0.5*cos(2*K*(x+g1*dt))+A*sin(w+g1*dt);
                x=x+temp2*(g1+g2)+epsilon*zz[k][t];
                w = w+temp2*(g1+g2);
                v = cos(K*(x))+0.5*cos(2*K*(x))+A*sin(w);
                //printf(" zz = %f\n", zz[k][t]);
                if(t>PasosT){
                vmedia[k] = vmedia[k] + v;
                }
            vmedia[k] = vmedia[k]/((tiempo/dt)-PasosT);
            }

                D[0]+=vmedia[k]; //velocidad

        }
        printf("%f\t %f\n",T,D[0]/((float)trayec));
        fprintf(D1,"%f\t %f\n",T,D[0]/((float)trayec));
    }

    fclose(D1);
    fclose(D2);
    return 0;
}
void ini_ran(int SEMILLA)
{
    srand(SEMILLA);
    for(int k=0;k<256;k++)Wheel[k]=(rand()<<16)+rand();
    ind_ran=ig1=ig2=ig3=0;
}
float RandomC(float Max,float Min)/*Modifico rand() con Parisi-Rapuano*/
{
    float r;
    ig1=ind_ran-11;ig2=ind_ran-29;ig3=ind_ran-17;
    Wheel[ind_ran]=Wheel[ig1]+Wheel[ig2];
    ir1=Wheel[ind_ran]^Wheel[ig3];
    ind_ran++;
    r=Min+(Max-Min)*ir1*(2.3283063671E-10F);
    return r;
}
float box_muller(float m, float s)/*Con media M y varianza s*/
{
    register float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last){y1 = y2; use_last = 0;}
    else
    {
        do {
            x1 = RandomC(1,-1);
            x2 = RandomC(1,-1);
            w = x1 * x1 + x2 * x2;
        } while(w>=1.0);

        w = sqrt( (-2.0*log(w))/w);
        y1 = x1*w;
        y2 = x2*w;
        use_last=1;
    }
    return(m + y1*s);
}
