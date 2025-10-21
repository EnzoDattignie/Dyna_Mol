#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

    static double L = 4;
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

    double x = -2.3;
    printf("%f",modulo(x));


}