#include <stdio.h>
#include <stdlib.h>
#include "pairingProfile.h"

int main (int argc, char *argv[]) {
    if((argc < 4) || (argc > 6)){
        fprintf(stderr, "Wrong number of arguments\n");
        fprintf(stderr, "Example usage: %s filePath filePathPtn filePathAlpha [runNumber [matchingAlgorithm]] ", argv[0]);
        exit(-1);
    }

    int runNumber = 5;
    long matching_algorithm = 1;

    if(argc > 4){
        sscanf(argv[4], "%d", &runNumber);
    }

    if(argc > 5){
        sscanf(argv[5], "%ld", &matching_algorithm);
    }

    pairingProfile(argv[1], argv[2], argv[3], runNumber, matching_algorithm);
}