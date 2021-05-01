#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>
#include "pairingProfile.h"

#define ARRAYSIZE 100

void profileFolder(const char *folderPath, int runNumber) {
    unsigned int eg2Len, currLen;
    int fileCount = 0;
//instead of realloc iterate through the folder multiple times by calling opendir
//use strdup
    char *eg2Arr[ARRAYSIZE], *ptnArr[ARRAYSIZE], *alphaArr[ARRAYSIZE];
    char *eg2Path, *ptnPath, *alphaPath, *currPath, *currentDirectory;
    char graphExt[] = ".eg2";
    char ptnExt[] = ".ptn";
    char alphaExt[] = ".alpha";
    unsigned int graphExtLen = strlen(graphExt);
    unsigned int ptnExtLen = strlen(ptnExt);
    unsigned int alphaExtLen = strlen(alphaExt);

    DIR *folder;
    DIR *iterFolder;

    struct dirent *currFile;
    struct dirent *iterFile;

    folder = opendir(folderPath);

    if (folder == NULL) {
        fprintf(stderr, "Could not open graph directory %s\n", folderPath);
        return;
    }

    // Sort files into eg2 or ptn arrays based on extension (last 3 characters)
    while ((currFile = readdir(folder)) != NULL) {
        eg2Path = currFile->d_name;
        eg2Len = strlen(eg2Path);
        if (strcmp(eg2Path + eg2Len - graphExtLen, graphExt)) {
            continue;
        }

        iterFolder = opendir(folderPath);
        ptnPath = alphaPath = NULL;

        while((iterFile = readdir(iterFolder)) != NULL) {
            currPath = iterFile->d_name;
            currLen = strlen(currPath);

            if ((ptnPath == NULL) && !strncmp(currPath, eg2Path, eg2Len-graphExtLen) && !strcmp(currPath + currLen - ptnExtLen, ptnExt)) {
                ptnPath = currPath;
            }

            if ((alphaPath == NULL) && !strncmp(currPath, eg2Path, eg2Len-strlen(graphExt)) && !strcmp(currPath + currLen - alphaExtLen, alphaExt)) {
                alphaPath = currPath; 
            }

            if (ptnPath && alphaPath) {
                break;
            }
        }
        closedir(iterFolder);
        if (ptnPath && alphaPath) {
            printf("Graph %s\n", eg2Path);
            eg2Arr[fileCount] = strdup(eg2Path);
            ptnArr[fileCount] = strdup(ptnPath);
            alphaArr[fileCount] = strdup(alphaPath);
            currentDirectory = getcwd(NULL, 0);
            chdir(folderPath);
            pairingProfile(eg2Arr[fileCount], ptnArr[fileCount], alphaArr[fileCount], runNumber, 0);
            pairingProfile(eg2Arr[fileCount], ptnArr[fileCount], alphaArr[fileCount], runNumber, 1);
            chdir(currentDirectory);
            free(currentDirectory);
            fileCount++;
        }
    }

    closedir(folder);
}

int main(int argc, char *argv[]) {
    if ((argc < 2) || (argc > 3)) {
        fprintf(stderr, "Usage: %s folderpath\n [runNumber]", argv[0]);
        exit(-1);
    }

    int runNumber = 5;
    if (argc > 2) {
        sscanf(argv[2], "%d", &runNumber);
    }

    profileFolder(argv[1], runNumber);
}
