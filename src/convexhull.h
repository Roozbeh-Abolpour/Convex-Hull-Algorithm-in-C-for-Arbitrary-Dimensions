#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include "MatLib.h"

struct convexhull{
    Mat* facets;
    Mat* neighbors_indices;
    Mat* outpoints_indices;
    Mat* points;
    Mat* center;
    int dim;
};

#ifdef __cplusplus
extern "C" {
#endif

convexhull* convh(Mat* points);
void init(convexhull* cvh);
Mat* normal(convexhull* cvh,int facet_index);
Mat* facetcenter(convexhull* cvh,int facet_index);
bool isoutside(convexhull* cvh,int facet_index,int point_index);
int furthestpoint(convexhull* cvh,int facet_index);
void vectoradd(Mat* A,double x);
void localsearch(convexhull* cvh,Mat* local_facets,Mat* boundary_facets,int facet_index,int prev_facet_index,int point_index);
void deleterows(Mat* A,Mat* rows);
void deletefacets(convexhull* cvh,Mat* facets);
void updatenewfacets(convexhull* cvh,Mat* local_facets,Mat* boundary_facets,int Nfb);

#ifdef __cplusplus
}
#endif

#endif
