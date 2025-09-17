#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double sigma = 0.01;
double epsilon = 1;
#define M 1;
static double dt = 0.1;
static double L = 20;
static int N = 1000;
static double T = 1;
static double t_max = 1;

int main () {
    srand(time(0));
    struct Part {
        double x;
        double y;
        double vx;
        double vy;
        double ax;
        double ay;
        double m;
    };

    int constructeur(struct Part *P,double x, double y) {
        P->x = x;
        P->y = y;
        P->ax = 0;
        P->ay = 0;
        P->vx = 0;
        P->vy = 0;
        P->m = M;
        return 0;
    }

    int afficher(struct Part *P) {
        printf("x = %f : y = %f\n",P->x,P->y);
        printf("vx = %f : vy = %f\n",P->vx,P->vy);
        printf("ax = %f : ay = %f\n\n",P->ax,P->ay);
        return 0;
    }

    double rng() {
        double rd = rand();
        rd = (rd/RAND_MAX - 0.5)*L;
        return rd;
    }

    double modulo(double x){
        double res = x;
        if (res>L/2) {
            long div = (long) (res+L/2)/L;
            res = res-L*(double) div;
            // printf("%f : %li : %f\n",res,div,x);
        } if (res<-L/2) {
            long div = (long) (-res+L/2)/L;
            res = res+L*(double) div;
            // printf("%f : %li : %f\n",res,div,x);
        }
        return res;
    }

    int update(struct Part *P) {
        P->vx = P->vx + P->ax * dt;
        P->vy = P->vy + P->ay * dt;
        P->x = modulo(P->x + P->vx * dt);
        P->y = modulo(P->y + P->vy * dt);
        P->ax = 0;
        P->ay = 0;
        return 0;
    }

    int Force(struct Part *P1, struct Part *P2) {
        double dx = P2->x-P1->x;;
        double dy = P2->y-P1->y;
        double r = sqrt(dx*dx + dy*dy);
    
            int Force_r(double* Fr, double r) {
                *Fr = 24*epsilon*(2*pow(sigma/r,12)/r - pow(sigma/r,6)/r);
                return 0;
            }
    
        double Fr = 0;
        Force_r(&Fr,r);
        // printf("F_r : %f\n",Fr);
        
        double Fx = Fr * dx/r;
        double Fy = Fr * dy/r;
    
        P1->ax = P1->ax - Fx;
        P1->ay = P1->ay - Fy;
        P2->ax = P2->ax + Fx;
        P2->ay = P2->ay + Fy;
        return 0;
    
        }

    

    int Force_liste(struct Part Liste[]) {
        for (int i = 0;i<N-1;i++) {
            for (int j = i+1;j<N;j++) {
                Force(&Liste[i],&Liste[j]);
            }
        }
        return 0;
    }

    int update_Liste(struct Part Liste[]) {
        for (int i = 0;i<N;i++) {
            update(&Liste[i]);
        }
    }

    struct Part Liste[N];
    for (int i = 0; i < N;i++) {
        constructeur(&Liste[i],rng(),rng());
    }

    for (double t = 0; t < t_max;t=t+dt) {
        Force_liste(Liste);
        for (int i = 0; i < N;i++) {
            afficher(&Liste[i]);
        }
        update_Liste(Liste);
        printf("Temps : %f\n",t);
    }



    for (int i = 0; i < N;i++) {
        afficher(&Liste[i]);
    }
    

    // constructeur(&P2,1.5,0);
    // double t;
    // static double t_max = 10;
    // for (t = 0; t<T; t = t + 0.1) {
    //     Force(&P1,&P2);
    //     afficher(&P1);
    //     afficher(&P2);
    //     update(&P1);
    //     update(&P2);
    //     printf("%f : %f\n",rng(),rng());
    // }
    // afficher(&P1);
    // afficher(&P2);
    printf("%f : %f : %f",modulo(0.4),modulo(234),modulo(-654.3));

    return 0;
}