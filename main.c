#include <stdio.h>
#include "lapacke.h"
#include "polytope_optimized.h"

int main(int argc, char* argv[]) {
    printf("--------------- PolyVest, but faster ----------------\n");
    // pritnf("If you have any questions or if you have found some bugs," <<
    // endl << "please contact me at <gecj@ios.ac.cn>." << endl;
    printf("=========================================\n");

    if (argc < 2 || argc > 4) {
        printf("ERROR: Invalid arguments.\n");
        printf("USAGE: %s <input-file>  <step-size-coef> [output-file]\n",
               argv[0]);
        return 1;
    }

    FILE* input = fopen(argv[1], "rt");

    if (!input) {
        printf("ERROR: Can't open target file\n");
        return 1;
    }

    polytope* p = NULL;

    read_polytope(&p, input);

    // print_polytope(&p);

    preprocess(p);

    // Default coef
    int coef = 1600;
    double vol = estimate_volume(p, coef);

    printf("volume: %lf\n", vol);

    return 0;
}