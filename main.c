#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// ==== Initialisation des constantes ==== //


int const seed = 220;

#define N (4) //Nombre total de particule
double sigma = 1;
double epsilon = 1;
double seuil = 0.8;

char nom_fichier[] = "./temp.txt"; //Nom du fichier d'enregistrement

double t_star_defaut = 1; 
double n_pas_defaut = 1e5;
static double t_max = 0.5;

static double M = 1;
static const double n = sqrt(N);


static const double L = 4;

// static double dl = 0.4 ;
static const double dl = L/n; //Ne marche que si N est un carré
static const double Rc = 2.5;

int main (int argc, char *argv[]) {
    srand(seed);
    double val;
    int mode;
    double t_val;
    if (argc == 4) {
        sscanf(argv[1],"%lf",&val);
        sscanf(argv[2],"%d",&mode);
        sscanf(argv[3],"%lf",&t_val);
    } else {
        val = n_pas_defaut;
        mode = 1;
        t_val = t_star_defaut;
    }
    const double n_pas = val;
    const double t_star = t_val;
    const double dt = t_star/n_pas;
    // ==== Initialisation ==== //
    FILE *fichier;
    fichier = fopen(nom_fichier,"w");
    fprintf(fichier,"dt = %.12f\nTemps;Energie\n",dt);
    fclose(fichier);

    // ==== Fonctions ==== //
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

  
double modulo(double x){
    double res = x;
    if (res > L/2.0) {
        long div = (long) ((res + L/2.0) / L);
        res = res - L * (double) div;
    } 
    if (res < -L/2.0) {
        long div = (long) ((-res + L/2.0) / L);
        res = res + L * (double) div;
    }
    return res;
}

    int fct_u(struct Part *P1, struct Part *P2) {
        double x1 = P1->x;
        double x2 = P2->x;
        double y1 = P1->y;
        double y2 = P2->y;
        double dx = x2-x1;
        double dy = y2-y1;
        double r2 = (dx*dx + dy*dy);
        if (r2 > Rc*Rc) {
            if (x2-x1 > L/2) {
                x2 = x2 - L;
            } else {
                if (x2-x1 < -L/2) {
                    x2 = x2 + L;
                }
            }
            if (y2-y1 > L/2) {
                y2 = y2 - L;
            } else {
                if (y2-y1 < -L/2) {
                    y2 = y2 + L;
                }
            }
            dx = x2-x1;
            dy = y2-y1;
            r2 = (dx*dx + dy*dy);
        }
        if (r2 < Rc*Rc) {
        double r2i = 1/r2;
        double r6 = r2i*r2i*r2i;
        double r12 = r6*r6;
        double val_u = 4*(r12-r6);
        P1->u = P1->u + val_u;
        P2->u = P2->u + val_u;
        }
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
        double x1 = P1->x;
        double x2 = P2->x;
        double y1 = P1->y;
        double y2 = P2->y;
        double dx = x2-x1;
        double dy = y2-y1;
        double r2 = (dx*dx + dy*dy);
        if (r2 > Rc*Rc) {
            if (x2-x1 > L/2) {
                x2 = x2 - L;
            } else {
                if (x2-x1 < -L/2) {
                    x2 = x2 + L;
                }
            }
            if (y2-y1 > L/2) {
                y2 = y2 - L;
            } else {
                if (y2-y1 < -L/2) {
                    y2 = y2 + L;
                }
            }
            dx = x2-x1;
            dy = y2-y1;
            r2 = (dx*dx + dy*dy);
        }
        if (r2 < Rc*Rc) {
            double r2i = 1/r2;
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
        }
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

    int Euler_step(struct Part Liste[],double *E_cin, double *E_pot) {
        // printf("AAAA");
        update_u(Liste);
        somme_E(Liste,E_cin,E_pot);
        for (int i = 0;i<N;i++) {
            Liste[i].vx = Liste[i].vx + Liste[i].ax * dt;
            Liste[i].vy = Liste[i].vy + Liste[i].ay * dt;
            Liste[i].x = Liste[i].x + Liste[i].vx * dt;
            Liste[i].y = Liste[i].y + Liste[i].vy * dt;
            if (Liste[i].x*Liste[i].x > L*L*0.25) {
                Liste[i].x = modulo(Liste[i].x);
            }
            if (Liste[i].y*Liste[i].y > L*L*0.25) {
                Liste[i].y = modulo(Liste[i].y);
            }
            Liste[i].ax = 0;
            Liste[i].ay = 0;
        }
        Force_liste(Liste);
        return 0;
    }

    int Verlet_step(struct Part Liste[],double *E_cin, double *E_pot) {
        double old_ax[N];
        double old_ay[N];
        update_u(Liste);
        somme_E(Liste,E_cin,E_pot);
        for (int i = 0; i<N;i++) {
            Liste[i].x = Liste[i].x + Liste[i].vx*dt + 0.5*Liste[i].ax*dt*dt;
            if (Liste[i].x*Liste[i].x > L*L*0.25) {
                Liste[i].x = modulo(Liste[i].x);
            }
            Liste[i].y = Liste[i].y + Liste[i].vy*dt + 0.5*Liste[i].ay*dt*dt;
            if (Liste[i].y*Liste[i].y > L*L*0.25) {
                Liste[i].y = modulo(Liste[i].y);
            }
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
        int compteur = 0;
        double origin = (-L+dl)/2;
        double posx = origin;
        double posy = origin;
        for (int i = 0;i<N;i++) {
            compteur ++;
            constructeur(&Liste[i],posx,posy);
            posx = posx + dl;
            if (compteur == (int)n) {
                posy = posy + dl;
                posx = origin;
                compteur = 0;
            }
        }
        return 0;
    }

    int config_rdm(struct Part Liste[]) {
        double posx;
        double posy;
        double rdm1;
        double rdm2;
        int flag;
        for (int i = 0;i<N;i++) {
            flag = 0;
            rdm1 = (double)rand();
            rdm2 = (double)rand();
            posx = (rdm1/RAND_MAX - 0.5)*L;
            posy = (rdm2/RAND_MAX - 0.5)*L;
            for (int j = 0;j<i;i++) {
                
            }
            constructeur(&Liste[i],posx,posy);          
    }}


    // ==== Initialisation du système ==== //

    double E_cin = 0;
    double E_pot = 0;
    

    struct Part Liste[N];
    config_rdm(Liste);
    // constructeur(&Liste[0],-0.879013,0.495180);
    // constructeur(&Liste[1],1.789280,-0.639103);
    // constructeur(&Liste[2],-1.984904,-1.381889);
    // constructeur(&Liste[3],1.565183,-0.247625);

    for (int i = 0; i < N;i++) { 
        afficher(&Liste[i]);
    }
    Force_liste(Liste);
    update_u(Liste);
    somme_E(Liste,&E_cin,&E_pot);
    printf("Dt : %f, Energie initiale %f\nE_cin = %f : E_pot = %f\n",dt,E_cin+E_pot,E_cin, E_pot);

    fichier = fopen(nom_fichier,"a");
    fprintf(fichier,"%f;%f\n",0.0,E_cin+E_pot);
    
    // ==== Boucle d'itération principale ==== //
    double compteur = 0;
    double compteur2 = 0;
    for (double t = 0; t < t_max;t=t+dt) {
        // printf("Mode : %d",mode);
        if (mode == 0) {
            Euler_step(Liste,&E_cin,&E_pot);
            // printf("Euler");
        } else {
            // printf("Verlet");
            Verlet_step(Liste,&E_cin,&E_pot);
        }
        // printf("Temps : %f; Energie : %f\n",t,E_cin+E_pot);
        if (compteur >= t_max/20-0.5*dt) { //Permet d'afficher 20 valeurs
            printf("Temps : %f; Energie : %f\n",t,E_cin+E_pot);
    //          for (int i = 0; i < N;i++) {
    //     afficher(&Liste[i]);
    // }
            compteur = 0;
        }
        if (compteur2 >= t_max/10000-0.5*dt) { //Permet d'enregistrer 10000 valeurs
            fprintf(fichier,"%f;%f\n",t,E_cin+E_pot);
            compteur2 = 0;
        }
        compteur = compteur + dt;
        compteur2 = compteur2 + dt;
    }

    // ==== Fin du programme ==== //
    fclose(fichier);
    return 0;
}