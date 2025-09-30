#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double sigma = 1;
double epsilon = 1;
#define M 1;
static double dt = 0.01;
static double L = 4;
const int n = 3;
static int N = n*n;
// static int N = 2;
static double T = 1;
static double t_max = 10;
static double Rc = 2.5;

int main () {
    // FILE *fichier;
    // fichier = fopen("Mesure.txt","w");
    // fprintf(fichier,"Temps;Energie\n");
    // fclose(fichier);
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
        rd = (rd/RAND_MAX - 0.5)*L;
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

    int update(struct Part *P) {
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
        double x1 = P1->x;
        double x2 = P2->x;
        double y1 = P1->y;
        double y2 = P2->y;
        double dx = x2-x1;
        double dy = y2-y1;
        double r2 = (dx*dx + dy*dy);
        printf("r initial : %f\n",sqrt(r2));
        if ( r2 > Rc*Rc ) {
            if (x1 > 0 && x2 < 0) {
                x2 += L;
            } else {
                if (x1 < 0 && x2 > 0) {
                    x2 += -L;
                }
            }
            if (y1 > 0 && y2 < 0) {
                y2 += L;
            } else {
                if (y1 < 0 && y2 > 0) {
                    y2 += -L;
                }
            }
            dx = x2-x1;
            dy = y2-y1;
            r2 = (dx*dx + dy*dy);
            printf("Copie en x = %f, y = %f\nValeur de r2 : %f, Rc2 : %f\n",x2,y2,r2,Rc*Rc);
            
        }
        if ( r2 < Rc*Rc ) {
            double r2i = sigma*sigma/r2;
            double r6 = r2i*r2i*r2i;
            double r12 = r6*r6;
            double val_u = 4*(r12-r6);
            P1->u = P1->u + val_u;
            P2->u = P2->u + val_u;
        } else {
            printf("r toujours trop grand : %f",sqrt(r2));
        }

        return 0;
    }

    int update_u(struct Part Liste[]) {
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
            *E_cin = *E_cin+Liste[i].m*(Liste[i].vx*Liste[i].vx + Liste[i].vy*Liste[i].vy);
            *E_pot = *E_pot+Liste[i].u;
        }
    }



    int Force(struct Part *P1, struct Part *P2) {
        double dx = P2->x-P1->x;
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
        return 0;
    }

    int config_crist(struct Part Liste[]) {
        double dl = L/(n);
        int compteur = 0;
        double posx = (dl-L)/2;
        double posy = (dl-L)/2;
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

    double modulo2(double x, double y) {
        int div = x/y;
        return x-div*y;
    }


    double E_cin = 0;
    double E_pot = 0;
    

    struct Part Liste[2];
    constructeur(&Liste[0],-1,0);
    constructeur(&Liste[1],2,0);

    // config_crist(Liste);
    // constructeur(&Liste[0],0,0);
    // constructeur(&Liste[0],1,1);

    // for (int i = 0; i < N;i++) {
    //     constructeur(&Liste[i],rng(),rng());
    // }

    // update_u(Liste);
    // somme_E(Liste,&E_cin,&E_pot);
    // printf("Energie initiale %f\n",E_cin+E_pot);
    // fichier = fopen("Mesure.txt","a");
    // for (double t = 0; t < t_max;t=t+dt) {
    //     Force_liste(Liste);
    //     // for (int i = 0; i < N;i++) {
    //     //     afficher(&Liste[i]);
    //     // }
    //     update_Liste(Liste);
    //     update_u(Liste);
    //     somme_E(Liste,&E_cin,&E_pot);
    //     fprintf(fichier,"%f;%f\n",t,E_cin+E_pot);

    //     if (fabs(modulo2(t,1)-0.1)<0.1*dt) {
    //     printf("Temps : %f; Energie : %f\n",t,E_cin+E_pot);
    //     }
        
    // }
    // fclose(fichier);
    for (int i = 0; i < 2;i++) {
        afficher(&Liste[i]);
    }
    fct_u(&Liste[0],&Liste[1]);
    for (int i = 0; i < 2;i++) {
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
    // printf("%f : %f : %f",modulo(0.4),modulo(234),modulo(-654.3));

    return 0;
}