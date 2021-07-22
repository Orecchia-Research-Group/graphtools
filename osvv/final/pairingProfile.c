#include <stdio.h>
#include <stdlib.h>
#include "flow.h"
#include "timer.h"

void pairingProfile(const char *filePath, const char *filePathPtn, const char *filePathAlpha, int runNumber, long matching_algorithm) {
     
    int route_flag = 1;

    long nedges;
    long fflow;
    
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

    float time_init;
    float time_S1;
    float time_S2;
    float time_match;

    float time_run;
    float time_sum = 0;

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

    printf("Allocated graph variables\n");
    for (int i = 0; i < M; i++) {
        fscanf(eg2File, "%ld %ld %ld", &tails[i], &heads[i], &weights[i]);
        tails[i]++;
        heads[i]++;
        weights[i] *= cap_orig;
    }
    fclose(eg2File);
    printf("Read file\n");

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
    printf("Read mask and created additional edges\n");

    /* CALL HI_PR - modified to output flow - would prefer for hipr to allocate this memory*/

    n = N + 2;
    m = M + N;

    printf("Graph %s has %ld nodes and %ld edges\n", filePath, N, M);

    for(int i = 0; i < runNumber; i++){
        time_run = timer();
        #ifdef DEBUG
            hipr(n, m, tails, heads, weights, N + 1, N + 2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow, route_flag, matching_algorithm, &time_init, &time_S1, &time_S2, &time_match);
        #else
            hipr(n, m, tails, heads, weights, N + 1, N + 2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow, route_flag, matching_algorithm);
        #endif
        time_run = timer() - time_run;
        time_sum += time_run;
    }

    time_sum /= runNumber;
    fprintf(stderr, "Average hipr runtime: %.5f\n", time_sum);
    fflush(stderr);

    printf("Finished with hipr call\n");
    free(heads);
    free(tails);
    free(weights);
    if (mheads) free(mheads);
    if (mtails) free(mtails);
    if (mweights) free(mweights);
    free(output_set);
}
