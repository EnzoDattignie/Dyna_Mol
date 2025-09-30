#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


static char nom_fichier[] = "./temp.txt";

const double t_star = 1;
const double n_pas = 1e6  ;
static double dt = t_star/n_pas;

// static double dt = 0.001;

double sigma = 1;
double epsilon = 1;
static double M = 1;
static double L = 3;
const int n = 3;
static int N = 15;
static double dl = 0.70 ;
static double T = 1;
static double t_max = 10;

int main () {
    FILE *fichier;
    fichier = fopen(nom_fichier,"w");
    fprintf(fichier,"dt = %f\nTemps;Energie\n",dt);
    fclose(fichier);
    srand(time(0));
    struct Part {
        double x;
        double y;
        double vx;
        double vy;
        double ax;
        double ay;
        double u;
        double m;
    };

    int constructeur(struct Part *P,double x, double y) {
        P->x = x;
        P->y = y;
        P->ax = 0;
        P->ay = 0;
        P->vx = 0;
        P->vy = 0;
        P->u = 0;
        P->m = M;
        return 0;
    }

    int afficher(struct Part *P) {
        printf("x = %f : y = %f\n",P->x,P->y);
        printf("vx = %f : vy = %f\n",P->vx,P->vy);
        printf("ax = %f : ay = %f\n\n",P->ax,P->ay);
        printf("u = %f\n", P->u);
        return 0;
    }

    double rng() {
        double rd = rand();
        rd = (rd/RAND_MAX);
        return rd;
    }

    double modulo(double x,double y){
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

    int Euler(struct Part *P) {
        P->vx = P->vx + P->ax * dt;
        P->vy = P->vy + P->ay * dt;
        P->x = P->x + P->vx * dt;
        P->y = P->y + P->vy * dt;
        P->ax = 0;
        P->ay = 0;
        P->u = 0;
        return 0;
    }

    int fct_u(struct Part *P1, struct Part *P2) {
        double dx = P2->x-P1->x;
        double dy = P2->y-P1->y;
        double r2 = sigma*sigma/(dx*dx + dy*dy);
        double r6 = r2*r2*r2;
        double r12 = r6*r6;
        double val_u = 4*(r12-r6);
        P1->u = P1->u + val_u;
        P2->u = P2->u + val_u;
        return 0;
    }

    int update_u(struct Part Liste[]) {
        for (int i = 0; i < N; i++) {
            Liste[i].u = 0;
        }
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {
                fct_u(&Liste[i],&Liste[j]); 
            }
        }
        return 0;
    }

    int somme_E(struct Part Liste[], double *E_cin,double *E_pot) {
        *E_cin = 0;
        *E_pot = 0;
        for (int i = 0; i < N; i++) {
            *E_cin = *E_cin+0.5*Liste[i].m*(Liste[i].vx*Liste[i].vx + Liste[i].vy*Liste[i].vy);
            *E_pot = *E_pot+0.5*Liste[i].u;
        }
    }



    int Force(struct Part *P1, struct Part *P2) {
        double dx = P2->x-P1->x;
        double dy = P2->y-P1->y;
        double r2i = 1/(dx*dx + dy*dy);
        double r6i = r2i*r2i*r2i;
        double r12i = r6i*r6i;
        double Fr = 24*epsilon*(2*r12i - r6i)*r2i;
        // printf("F_r : %f\n",Fr);
        
        double Fx = Fr * dx;
        double Fy = Fr * dy;
    
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

    int Euler_Liste(struct Part Liste[]) {
        for (int i = 0;i<N;i++) {
            Euler(&Liste[i]);
        }
        return 0;
    }

    int Verlet_Liste(struct Part Liste[],double *E_cin, double *E_pot) {
        double old_ax[N];
        double old_ay[N];
        update_u(Liste);
        somme_E(Liste,E_cin,E_pot);
        for (int i = 0; i<N;i++) {
            Liste[i].x = Liste[i].x + Liste[i].vx*dt + 0.5*Liste[i].ax*dt*dt;
            Liste[i].y = Liste[i].y + Liste[i].vy*dt + 0.5*Liste[i].ay*dt*dt;
            old_ax[i] = Liste[i].ax;
            old_ay[i] = Liste[i].ay;
            Liste[i].ax = 0;
            Liste[i].ay = 0;
        }
        Force_liste(Liste);
        for (int i = 0; i<N;i++) {
            Liste[i].vx = Liste[i].vx + 0.5*(Liste[i].ax+old_ax[i])*dt;
            Liste[i].vy = Liste[i].vy + 0.5*(Liste[i].ay+old_ay[i])*dt;
        }
        return 0;
    }

    int config_crist(struct Part Liste[]) {
        double l = L/(n);
        int compteur = 0;
        double posx = (l-L)/2;
        double posy = (l-L)/2;
        for (int i = 0;i<N;i++) {
            compteur ++;
            constructeur(&Liste[i],posx,posy);
            posx = posx + dl;
            if (compteur == n) {
                posy = posy + dl;
                posx = (dl-L)/2;;
                compteur = 0;
            }
        }
        return 0;
    }

    int config_crist2(struct Part Liste[]) {
        int compteur = 0;
        double posx = 0;
        double posy = 0;
        for (int i = 0;i<N;i++) {
            compteur ++;
            constructeur(&Liste[i],posx,posy);
            posx = posx + dl;
            if (compteur == n) {
                posy = posy + dl;
                posx = 0;
                compteur = 0;
            }
        }
        return 0;
    }

    double modulo2(double x, double y) {
        int div = x/y;
        return x-div*y;
    }


    double E_cin = 0;
    double E_pot = 0;
    

    struct Part Liste[N];
    config_crist2(Liste);
    // constructeur(&Liste[0],0,0);
    // constructeur(&Liste[0],1,1);

    // for (int i = 0; i < N;i++) {
    //     constructeur(&Liste[i],rng(),rng());
    // }

    update_u(Liste);
    somme_E(Liste,&E_cin,&E_pot);

    // double dE = (E-(E_cin+E_pot))/N;
    // double dv = sqrt(dE/M);
    // for (int i = 0; i < N;i++) {
    //     double theta = rng()*2*M_PI;
    //     Liste[i].vx = Liste[i].vx+dv*cos(theta);
    //     Liste[i].vy = Liste[i].vx+dv*sin(theta);
    // }


    printf("Energie initiale %f\nE_cin = %f : E_pot = %f\n",E_cin+E_pot,E_cin, E_pot);
    fichier = fopen(nom_fichier,"a");
    fprintf(fichier,"%f;%f\n",0.0,E_cin+E_pot);
    double compteur = 0;
    Force_liste(Liste);
    for (double t = 0; t < t_max;t=t+dt) {
        // update_u(Liste);
        // somme_E(Liste,&E_cin,&E_pot);
        // Euler_Liste(Liste);
        // Force_liste(Liste);
        Verlet_Liste(Liste,&E_cin,&E_pot);
        if (compteur >= t_max/20-0.5*dt) {
            fprintf(fichier,"%f;%f\n",t,E_cin+E_pot);
            printf("Temps : %f; Energie : %f\n",t,E_cin+E_pot);
            compteur = 0;
        }
        compteur = compteur + dt;
        
    }
    fclose(fichier);
    // for (int i = 0; i < N;i++) {
    //     afficher(&Liste[i]);
    // }
    

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
    // printf("%f : %f : %f",modulo(0.4),modulo(234),modulo(-654.3));

    return 0;
}