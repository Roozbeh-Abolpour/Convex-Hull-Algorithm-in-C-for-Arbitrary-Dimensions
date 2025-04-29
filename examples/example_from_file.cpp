#include <stdio.h>
#include "src/MatLib.h"
#include "src/convexhull.h"

int main(int argc,char* argv[]){
    Mat* V=readmat("file1.txt");
    showmat(V);
    convexhull* cvh=convh(V);
    writemat("file2.txt",cvh->facets);
    writemat("file3.txt",cvh->neighbors_indices);
    writemat("file4.txt",cvh->outpoints_indices);
    freemat(V);
    freemat(cvh->facets);
    freemat(cvh->neighbors_indices);
    freemat(cvh->outpoints_indices);
    freemat(cvh->center);
    free(cvh);
    return 0;
}
