#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define tiempo (float)(7000)//Tiempo final
#define dt (float)(0.01)//Paso en tiempo
#define Afin (float)(7)//Tiempo final
#define dA (float)(0.01)//Paso en amplitud
#define trayec (int)(1)
#define K   (float)(2*3.141592)
#define D   (float)(0)
#define omega   (float)(0.25)
#define PasosTerm (float)(100)


void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float box_muller(float m, float s);
float f(float var1,float var2, float A);

unsigned char ind_ran,ig1,ig2,ig3;
unsigned int Wheel[256],ir1;
FILE *D1;
int main()
{
    ini_ran(123456789);
    D1=fopen("omega_1.txt","w");
    float x,v,v_med,epsilon,**zz;

	zz = (float **)malloc(trayec*sizeof(float*));
	for (int i=0;i<trayec;i++)zz[i]=(float*)malloc((int)(tiempo/dt)*sizeof(float));

    for(float A=0;A<=Afin;A+=dA)/*Para varias temperaturas*/
    {
        v_med=0;epsilon=pow(dt*2*D,0.5);
        for(int k=0;k<trayec;k++)
        {
            x=v=0;
            for(int t=0;t<(int)(tiempo/dt);t++)
            {
                register float g12,g22;
                if(A==0&&k<trayec&&t<(int)(tiempo/dt))zz[k][t]=box_muller(0,1);
                g12=f(omega*t*dt,x+epsilon*zz[k][t],A);   g22=f(omega*(t+1)*dt,x+epsilon*zz[k][t],A);
                x+=0.5*dt*(g22+g12)+epsilon*zz[k][t];
                //v+=x;
                if(t>PasosTerm)
                v+=f(omega*(t+1)*dt,x,A);
            }
            v_med+=v/(float)(tiempo/dt);
        }
        printf("%f\t %f\t\n",A,v_med);

        fprintf(D1,"%f\t %f\n",A,v_med);
    }

    fclose(D1);
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
float f(float var1,float var2, float A){return cos(K*var2)+0.5*cos(2*K*var2)+A*sin(var1);}//var1->Omega  var2->x
