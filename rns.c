#include "stdio.h"
#include "time.h"
#include <Windows.h>
#include "unistd.h"
#include "math.h"

void findRNS(int64_t*,int64_t*,int64_t,unsigned int);
int64_t findCRT(int64_t*,int64_t*,unsigned int);
int64_t findGCD(int64_t, int64_t, int64_t*, int64_t*);

//Parameters pointing to moduli, array to store RNS rep and input x
void findRNS(int64_t* m,int64_t* X,int64_t x,unsigned int k){
    int i;
    for(i=0;i<k;i++){
        X[i]=x % m[i];
    }
}

int64_t findCRT(int64_t* n,int64_t* Y,unsigned int k){
    int64_t i,j,m=1,m_i[k],c[k],t,gcd,value=0,c_ij[k][k],r[k][k],prod;
    for(i=0;i<k;i++){
        m = m * n[i];
    }
    for(i=0;i<k;i++){
        m_i[i] = m/n[i];
        gcd=findGCD(m_i[i],n[i],&c[i],&t);
        c[i] = (c[i]+n[i])%n[i];
    }
    /* // This implementation uses  multi-precision arithmetic.. so implementing MRC
    for(i=0;i<k;i++){
        value+= Y[i]*c[i]*m_i[i];
    } */
    for(j=0;j<k;j++) r[j][0]=Y[j];

    for(i=1;i<k;i++){
        for(j=0;j<i;j++){
            gcd=findGCD(n[j],n[i],&c_ij[i][j],&t);
            c_ij[i][j]= (c_ij[i][j]+n[i])%n[i];
            r[i][j+1]=((((r[i][j]-r[j][j])*c_ij[i][j])%n[i])+n[i])%n[i]; //solve this issue later... doesn't work if r[i][j]-r[j][j] goes much negative
            //printf("c_ij[i][j] - %llu , rij - %llu \n",c_ij[i][j],r[i][j+1]);
            }
    }
    for(i=0;i<k;i++){
        prod=1;
        for(j=0;j<i;j++){
            prod*=n[j];
        }
        //printf("rij - %llu , prod- %llu",r[i][i],prod);
        value+=r[i][i]*prod;
    }
    return value%m;
}

// C function for extended Euclidean Algorithm
int64_t findGCD(int64_t a, int64_t b, int64_t *x, int64_t *y){
    // Base Case
    if (a == 0)
    {
        *x = 0;
        *y = 1;
        return b;
    }

    int64_t x1, y1; // To store results of recursive call
    int64_t gcd = findGCD(b%a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b/a) * x1;
    *y = x1;

    return gcd;
}

main() {

    //Procedure:
    //First, find n' and r^-1 using EEA
    //decide on RNS moduli set (how many m_i's, k=?) 1280bit RNS arithmetic? 20m_i's.. store them globally
    //Implement functions to find RNS and CRT representations
    //

    unsigned int k=5; //value used in the slides
    int i;
    //printf("sizeof(long long int) = %d bytes\n", (long long int) sizeof(long long int));
    //Create an array to store moduli - use long long int to get 8byte storage
    int64_t m[]={3,5,7,11,13};
    //Declaring variables to store inputs
    int64_t a=19,b=21,n=29,r,kk,gcd;
    //Declaring arrays to store RNS representation of inputs
    int64_t A[k],B[k],N[k],R[k],N_bar[k],R_inv[k],T[k],Q[k],T_bar[k]; //Note to prof: Mistake in Step5 - Should be T = (A.B+Q.N).R^-1. Defined T_bar[] to store A.B
    findRNS(m,A,a,k); //findRNS checked and working
    findRNS(m,B,b,k);
    findRNS(m,N,n,k);
    /*for(i=0;i<k;i++){
        printf("%llu \n",B[i]);
    } */
    //CRT checked and working?
    int64_t value,nn[]={7,9,11},Z[]={2,2,3,6,4},r_hat,n_hat,n_bar,t,q;

    kk=floor(log(n)/log(2)+1);
    r=pow(2,kk);
    findRNS(m,R,r,k);
    gcd=findGCD(r,n,&r_hat,&n_hat);
    n_bar= -n_hat;
    findRNS(m,N_bar,n_bar,k);
    for(i=0;i<k;i++){
        gcd=findGCD(R[i],m[i],&R_inv[i],&n_hat);
        R_inv[i]=(R_inv[i]+m[i])%m[i];
    }
    double time;
    clock_t start=clock();
    //Step1 of RNS
    for(i=0;i<k;i++){
        T_bar[i]=(A[i]*B[i])%m[i];
        T[i]=(A[i]*B[i]*N_bar[i])%m[i];
        //printf("%lld %lld %lld %lld \n",A[i],B[i],N_bar[i],T[i]);
    }
    //Step2
    t=findCRT(m,T,k);
    //Step3
    q=t%r;
    //Step4
    findRNS(m,Q,q,k);
    //Step5
    for(i=0;i<k;i++){
        T[i]=((T_bar[i]+Q[i]*N[i])*R_inv[i])%m[i];
        //printf("%lld %lld %lld %lld \n",Q[i],N[i],R_inv[i],T[i]);
    }
    //convert back from rns
    t=findCRT(m,T,k);
    printf("%lld \n",t);
    Sleep(100);
    clock_t stop=clock();

    //

    time = (double)(stop - start)/CLOCKS_PER_SEC;
    printf("Time : %f",time);
    return 0;
}
