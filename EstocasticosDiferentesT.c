#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tiempo (float)(100)//Tiempo final
#define dt (float)(0.01)//Paso en tiempo
#define trayec (int)(10000)
#define M   (float)(0.1)
#define Eta (float)(1.0)
#define K   (float)(0.0)
//#define T   (float)(1.0)
//#define epsilon   (pow((dt*2*T*Eta/(M*M)),0.5))

void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float box_muller(float m, float s);
void avanza (float *x,float *v, float epsilon);


unsigned char ind_ran,ig1,ig2,ig3;
unsigned int Wheel[256],ir1;
FILE *D1;
FILE *D2;
int main()
{
    ini_ran(123456789);
    D1=fopen("1.txt","w");
    D2=fopen("2.txt","w");
    float x,v;
    float XX_[trayec],E_c[trayec];
    double Ek[(int)(tiempo/dt)];
    double Ek_m;
    float Tini, Tfinal, deltaT;
    float epsilon;
    int k,t, pasosT,i;

    /// Variación Temperatura
    Tini = 0;
    Tfinal = 10;
    deltaT=0.5;
    pasosT = (int)((Tfinal-Tini)/deltaT); /// N de temperaturas
    float T[pasosT], Te, E_cin[pasosT],Difusion[pasosT];
    //------------------------
    //printf("pasosT = %d \n",pasosT);

    for(i=0; i<=pasosT; i++)
        {
            T[i] = Tini + deltaT*i;
            //printf("T[%d] = %f \n",i,T[i]);
            Te = T[i];
            epsilon = pow((dt*2*T[i]*Eta/(M*M)),0.5);

            for(k=0;k<trayec;k++)
                {
                    Ek_m = 0;
                    x=v=0;
                    for(t=0;t<(int)(tiempo/dt);t++)
                        {
                            avanza(&x,&v,epsilon);
                            Ek[t]=0.5*M*v*v;
                            Ek_m = Ek_m+Ek[t]/(tiempo/dt);
                        }
                    XX_[k]=x*x;
                    E_c[k]=Ek_m;
                }

            for(k=0; k<=pasosT; k++)
                {
                    E_cin[k]=0;
                    Difusion[k]=0;
                }

            for(k=200;k<trayec;k++)
                {
                    E_cin[i]+=E_c[k]/(float)(trayec);
                    Difusion[i]+=XX_[k]/(float)(trayec);
                    // fprintf(D1,"%d\t%f\n",k,E_cin);

                }
            printf("<E(T=%f)>=%f\t D=%f\n",Te,E_cin[i],Difusion[i]/(2*tiempo));
            fprintf(D1,"%f\t%f\t%f\n",T[i],E_cin[i],Difusion[i]/(2*tiempo));

        }

    fclose(D1);
    fclose(D2);
    return 0;
}

void ini_ran(int SEMILLA)
{
    int k;
    srand(SEMILLA);
    for(k=0;k<256;k++)Wheel[k]=(rand()<<16)+rand();
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
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last){y1 = y2; use_last = 0;}
    else
    {
        do {
            x1 = RandomC(1,-1);
            x2 = RandomC(1,-1);
           // printf("%f\n",x1);
            w = x1 * x1 + x2 * x2;
        } while(w>=1.0);

        w = sqrt( (-2.0*log(w))/w);
        y1 = x1*w;
        y2 = x2*w;
        use_last=1;
    }
    return(m + y1*s);
}
void avanza (float *x,float *v, float epsilon)
{

    float g11,g12,g21,g22,z;

    z=box_muller(pow(10,-6),1);

    //printf("%f\n",z);
    g11=*v+epsilon*z;  /*Calculo las funciones g11,g12,... para el rk estocástico*/
    g12=-(Eta/M)*g11-(K/M)*(*x);
    g21=*v+g12*dt;
    g22=-(Eta/M)*g21-(K/M)*(*x+g11*dt);
    *x=*x+(dt/2)*(g11+g21);    /*Calculo posiciones y velocidades en cada momento*/
    *v=*v+(dt/2)*(g12+g22)+epsilon*z;
}





