#include <stdio.h>
#include <stdlib.h>
#include "flow.h"

void pairingProfile(const char *filePath, const char *filePathPtn, const char *filePathAlpha, long matching_algorithm) {
    
    //read alpha as cap_add cap_orig (two longs)
    //all variable declarations to the top (int, long, pointer, double, FILE)
     
    int route_flag = 1;

    long nedges;
    long fflow;
    long size_cut;
    long matching_index = matching_algorithm;

    long N, n, M, m, temp;
    long cap_add, cap_orig;

    int *output_set;
    long *mheads = NULL;
    long *mtails = NULL;
    long *mweights = NULL;
    long *tails;
    long *heads;
    long *weights;
    int *mask;

    FILE *eg2File;
    FILE *ptnFile;
    FILE *alphaFile;
 
    alphaFile = fopen(filePathAlpha, "r");
    if (!alphaFile) {
        fprintf(stderr, "Unable to open file: %s\n", filePathAlpha);
    }

    fscanf(alphaFile, "%ld %ld", &cap_add, &cap_orig);

    eg2File = fopen(filePath, "r");
    if (!eg2File) {
        fprintf(stderr, "Unable to open file: %s\n", filePath);
        return;
    }

    fscanf(eg2File, "%ld %ld %ld", &N, &temp, &M);
 
    tails = calloc(sizeof(*tails), M + N);
    heads = calloc(sizeof(*heads), M + N);
    weights = calloc(sizeof(*weights), M + N);
    mask = calloc(sizeof(*mask), N + 1);

    for (int i = 0; i < M; i++) {
        fscanf(eg2File, "%ld %ld %ld", &tails[i], &heads[i], &weights[i]);
        tails[i]++;
        heads[i]++;
        weights[i] *= cap_orig;
    }
    
    fclose(eg2File);

    ptnFile = fopen(filePathPtn, "r");
    if (!ptnFile) {
        fprintf(stderr, "Unable to open file: %s\n", filePathPtn);
    }

    for (int h = 0; h < N; h++) {
        fscanf(ptnFile, "%d", &mask[h+1]);
        if (mask[h+1] == 1) {
            tails[M+h] = N+1;
            heads[M+h] = h+1;
        } else {
            tails[M+h] = h+1;
            heads[M+h] = N+2;
        }
        weights[M+h] = cap_add;
    }

    /* CALL HI_PR - modified to output flow - would prefer for hipr to allocate this memory*/

    n = N + 2;
    m = M + N;


    hipr(n, m, tails, heads, weights, N + 1, N + 2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow,
         route_flag);

    free(heads);
    free(tails);
    free(weights);
    if (mheads) free(mheads);
    if (mtails) free(mtails);
    if (mweights) free(mweights);
    free(output_set);
}
