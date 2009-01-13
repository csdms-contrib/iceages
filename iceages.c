#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define maxoutputlength 10000000

#define FREE_ARG char*
#define NR_END 1

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(idum)
int *idum;
{
        static int iset=0;
        static float gset;
        float fac,r,v1,v2;
        float ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

main()
{    FILE *fp1;
     int n,t,t2,idum,duration,nsteps;
     float tau,per,drift,eta,*T,*B,Bh,*v,deltaT,timestep;

     fp1=fopen("iceagemodel","w");
     T=vector(1,maxoutputlength);   
     B=vector(1,maxoutputlength);
     v=vector(1,maxoutputlength);
     eta=0.9;
     tau=75;
     drift=0.2;
     T[1]=0;B[1]=0;v[1]=0;
     idum=-7834;
     timestep=100;      /* yr */
     duration=3000000;  /* yr */
     nsteps=(int)(duration/timestep);
     for (t=1;t<=nsteps;t++) 
      {if (T[t]<-1) deltaT=-T[t]-3.2*sqrt(-1-T[t])+eta*gasdev(&idum)+drift; 
        else deltaT=-T[t]+eta*gasdev(&idum)+drift;
       Bh=0;
       for (t2=t-2*tau;t2<t;t2++)
        if (t2>1) Bh+=-(T[t2]-1)*exp(-fabs(t-t2)/tau)/tau; 
       B[t+1]=Bh;
       if (T[t]<-1) v[t+1]=fabs(B[t+1]-B[t]); else v[t+1]=0; 
       deltaT+=25*v[t+1];
       T[t+1]=T[t]+deltaT;
       fprintf(fp1,"%d %f %f %f\n",t*100.0,T[t],v[t+1],Bh);}
     fclose(fp1);
}     
