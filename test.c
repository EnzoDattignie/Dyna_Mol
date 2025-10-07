#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
    printf("argc : %d\n",argc);
    printf("argv : %s\n",argv[1]);
    double val;
    sscanf(argv[1],"%lf",&val);
    printf("val = %lf",val);

}