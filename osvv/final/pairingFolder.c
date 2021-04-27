#include "pairingProfile.c"
#include "flow.h"
#include <dirent.h>
#include <string.h>

#define ARRAYSIZE 100

void main(const char *folderPath, const char *alphaPath, int runNumber){
    int eg2_len=0;
    int ptn_len=0;
    int currLen, iterLen;
    int fileCount = 0;
//instead of realloc iteratre through the folder multiple times by calling opendir
//use strdup
    char *eg2Arr[ARRAYSIZE], *ptnArr[ARRAYSIZE];
    int eg2Len[ARRAYSIZE], ptnLen[ARRAYSIZE];
    char currName[1000], iterName[1000];
    char graphExt[] = "eg2";
    char ptnExt[] = "ptn";

    DIR *folder;
    DIR *iterFolder;

    struct dirent *currFile;
    struct dirent *iterFile;

    folder = opendir(folderPath);

    if (folder == NULL){
        fprintf(stderr, "Could not open this directory");
        return;
    }

    // Sort files into eg2 or ptn arrays based on extension (last 3 characters)
    while ((currFile = readdir(folder)) != NULL){
        strcpy(currName, currFile -> d_name);
        currLen = strlen(currName);
    
        iterFolder = opendir(folderPath);

        if(!strcmp(currName+currLen-strlen(graphExt), graphExt)){
            while((iterFile = readdir(iterFolder)) != NULL){
                strcpy(iterName, iterFile -> d_name);
                iterLen = strlen(iterName);

                if(!strcmp(iterName+currLen-strlen(ptnExt), ptnExt) && !strncmp(currName, iterName, currLen-strlen(graphExt))){
                    eg2Arr[fileCount] = currName;
                    eg2Len[fileCount] = currLen;

                    ptnArr[fileCount] = iterName;
                    ptnLen[fileCount] = iterLen;

                    fileCount++;

                    pairingProfile(currName, iterName, alphaPath, (long) 0);
                    pairingProfile(currName, iterName, alphaPath, (long) 1);
                }
            }
            closedir(iterFolder);
        }
    }
    closedir(folder);
}
