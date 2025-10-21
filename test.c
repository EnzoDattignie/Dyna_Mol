#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
    const double L = 2;
    double x = -63713144576421233.2;
    double modulo(double x){
        double res = x;
        if (res>L/2) {
            long div = (long) (res+L/2)/L;
            res = res-L*(double) div;
            printf("%f : %li : %f\n",res,div,x);
        } if (res<-L/2) {
            long div = (long) (-res+L/2)/L;
            res = res+L*(double) div;
            printf("%f : %li : %f\n",res,div,x);
        }
        return res;
    }
    modulo(x);
}