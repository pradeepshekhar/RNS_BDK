#include "stdio.h"
#include "time.h"
#include <Windows.h>
#include "sys/time.h"
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
        X[i]=(X[i]+m[i])% m[i];
    }
}

int64_t findCRT(int64_t* n,int64_t* Y,unsigned int k){
    int64_t i,j,m=1,m_i[k],c[k],t,gcd,value=0,c_ij[k][k],r[k][k],prod;
    for(i=0;i<k;i++){
        m = m * n[i];
    }
    /*
    for(i=0;i<k;i++){
        m_i[i] = m/n[i];
        gcd=findGCD(m_i[i],n[i],&c[i],&t);
        c[i] = (c[i]+n[i])%n[i];
    }*/
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

    int loop;
    double RNS=0,BDK=0;
    clock_t start,stop;

    for(loop=0;loop<100000;loop++){

        unsigned int k=5;
        int64_t m[]={3,5,7,11,13},m_BDK,M[]={3,7,13,19,29},p_BDK,P[]={5,11,17,23,31},m_P_inv;
        //Declaring variables to store inputs
        //int64_t i,j,prod,a=19,b=21,n=29,r,kk,gcd,c=26386,d=72931,o=14527;
        int64_t i,j,prod,a=26386,b=72931,n=14527,r,kk,gcd,c=26386,d=72931,o=14527;
        //Declaring arrays to store RNS representation of inputs
        int64_t A[k],B[k],N[k],R[k],N_bar[k],R_inv[k],T[k],Q[k],T_bar[k],C_M[k],C_P[k],D_M[k],D_P[k],O_M[k],O_P[k],O_M_inv[k],C_M_neg[k],S[k],r_MRC[k][k],c_ij[k][k],S_MRC[k],M_ji[k-1][k],S_P[k],M_P_inv[k],T_P[k]; //Note to prof: Mistake in Step5 - Should be T = (A.B+Q.N).R^-1. Defined T_bar[] to store A.B
        findRNS(m,A,a,k); //findRNS checked and working
        findRNS(m,B,b,k);
        findRNS(m,N,n,k);
        /*for(i=0;i<k;i++){
            printf("%llu \n",B[i]);
        } */
        //CRT checked and working?
        int64_t value,nn[]={7,9,11},Z[]={2,2,3,6,4},r_hat,n_hat,n_bar,t,q;

        kk=floor(log(n)/log(2)+1);
        //r=pow(2,kk);
        r=150423;
        findRNS(m,R,r,k);
        gcd=findGCD(r,n,&r_hat,&n_hat);
        n_bar= -n_hat;
        if(n_bar<0){
            n_bar+=r;
        }
        //printf("n_bar - %lld \n",n_bar);
        findRNS(m,N_bar,n_bar,k);
        //Finding R^-1
        for(i=0;i<k;i++){
            gcd=findGCD(R[i],m[i],&R_inv[i],&n_hat);
            R_inv[i]=(R_inv[i]+m[i])%m[i];
        }


        //Step1 of RNS
        for(i=0;i<k;i++){
            start=clock();
            T_bar[i]=(A[i]*B[i])%m[i];
            T[i]=(A[i]*B[i]*N_bar[i])%m[i];
            //printf("%lld %lld %lld %lld \n",A[i],B[i],N_bar[i],T[i]);
            stop=clock();
        }
        RNS+=(double)(stop-start)/CLOCKS_PER_SEC;
        start=clock();
        //Step2
        t=findCRT(m,T,k);
        //Step3
        q=t%r;
        //printf("t - %lld q=%llu \n",t,q);
        //Step4
        findRNS(m,Q,q,k);
        stop=clock();
        RNS+=(double)(stop-start)/CLOCKS_PER_SEC;
        //Step5
        for(i=0;i<k;i++){
            start=clock();
            T[i]=((T_bar[i]+Q[i]*N[i])*R_inv[i])%m[i];
            //printf("%lld %lld %lld %lld \n",Q[i],N[i],R_inv[i],T[i]);
            stop=clock();
        }
        RNS+=(double)(stop-start)/CLOCKS_PER_SEC;
        //convert back from rns
        t=findCRT(m,T,k);
        //printf("%lld \n",t);
        //Sleep(0.1);
        //Implementation of BDK
        //Setup steps
        findRNS(M,C_M,c,k);
        findRNS(P,C_P,c,k);
        findRNS(M,D_M,d,k);
        findRNS(P,D_P,d,k);
        findRNS(M,O_M,o,k);
        findRNS(P,O_P,o,k);
        //printf("%lld %lld %lld %lld %lld\n",O_P[0],O_P[1],O_P[2],O_P[3],O_P[4]);
        for(i=0;i<k;i++){
            gcd=findGCD(O_M[i],M[i],&O_M_inv[i],&n_hat);
            O_M_inv[i]=(O_M_inv[i]+M[i])%M[i];
            //printf("%lld \n",O_M_inv[i]);
        }
        for(i=0;i<k;i++)  C_M_neg[i]=M[i]-C_M[i];
        for(i=0;i<k;i++) {
            start=clock();
            S[i]=(C_M_neg[i]*D_M[i]*O_M_inv[i])%M[i];
            stop=clock();
            //printf("%lld \n",S[i]);
        }
        BDK+=(double)(stop-start)/CLOCKS_PER_SEC;
        //Step2 - Basis Conversion from M to P
        start=clock();
        for(j=0;j<k;j++) r_MRC[j][0]=S[j];

        for(i=1;i<k;i++){
            for(j=0;j<i;j++){
                gcd=findGCD(M[j],M[i],&c_ij[i][j],&n_hat);
                c_ij[i][j]= (c_ij[i][j]+M[i])%M[i];
                r_MRC[i][j+1]=((((r_MRC[i][j]-r_MRC[j][j])*c_ij[i][j])%M[i])+M[i])%M[i];
                //printf("c_ij[i][j] - %llu , rij - %llu \n",c_ij[i][j],r_MRC[i][j+1]);
                }
        }

        for(i=0;i<k;i++){

            S_MRC[i]=r_MRC[i][i];
        }
        stop=clock();
        BDK+=(double)(stop-start)/CLOCKS_PER_SEC;
        //precomputing M_ji
        for(i=1;i<k;i++){
            prod=1;
            for(j=0;j<i;j++){
                prod*=M[j];
            }
            findRNS(P,&M_ji[i-1],prod,k);
            //printf("%lld %lld %lld %lld %lld\n",M_ji[i-1][0],M_ji[i-1][1],M_ji[i-1][2],M_ji[i-1][3],M_ji[i-1][4]);
            //value+=r[i][i]*prod;
        }
        //Comverting S from M to P
        for(i=0;i<k;i++){
            start=clock();
            S_P[i]=S_MRC[0];
            for(j=1;j<k;j++){
                S_P[i]+=M_ji[j-1][i]*S_MRC[j];
                //printf("%lld %lld \n",M_ji[j-1][i],S_MRC[j]);
            }
            S_P[i]=S_P[i]%P[i];
            //printf("%lld \n",S_P[i]);
            stop=clock();
        }

        BDK+=(double)(stop-start)/CLOCKS_PER_SEC;
        //Computation of T
        m_BDK=1;p_BDK=1;
        for(i=0;i<k;i++){
            m_BDK*=M[i];
            p_BDK*=P[i];
        }
        //printf("%lld \n",m_BDK);
        gcd=findGCD(m_BDK,p_BDK,&m_P_inv,&n_hat);
        findRNS(P,M_P_inv,m_P_inv,k);
        for(i=0;i<k;i++){
            start=clock();
            T_P[i]=((C_P[i]*D_P[i]+S_P[i]*O_P[i])*M_P_inv[i])%P[i];
            //printf("%lld \n",T_P[i]);
            stop=clock();
        }
        BDK+=(double)(stop-start)/CLOCKS_PER_SEC;

    }

    printf("Times : %f %f",RNS,BDK);
}
