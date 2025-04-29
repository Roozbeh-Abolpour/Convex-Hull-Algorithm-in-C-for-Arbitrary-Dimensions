#include "convexhull.h"
#include "MatLib.h"

struct convexhull{
	Mat* facets;
	Mat* neighbors_indices;
	Mat* outpoints_indices;
	Mat* points;	
	Mat* center;
	int dim;
};

Mat* normal(convexhull* cvh,int facet_index){
	int d=cvh->dim;	
	Mat* T=newmat(d-1,d,0);		
	for(int j=2;j<=d;j++){			
		for(int k=1;k<=d;k++){
			int p1_index=(int)cvh->facets->entries[(facet_index-1)*d+j-1];
			double x1=cvh->points->entries[(p1_index-1)*d+k-1];
			int p2_index=(int)(cvh->facets->entries[(facet_index-1)*d]);
			double x2=cvh->points->entries[(p2_index-1)*d+k-1];
			T->entries[(j-2)*d+k-1]=x1-x2;
		}
	}			
	Mat* N=null(T);
	Mat* center=cvh->center;
	double de=0;
	for(int j=1;j<=d;j++){
		int p_index=(int)cvh->facets->entries[(facet_index-1)*d];
		double x1=cvh->points->entries[(p_index-1)*d+j-1];
		de+=(center->entries[j-1]-x1)*N->entries[j-1];
	}
	if (de<0){
		scalermultiply2(N,N,-1);
	}
	freemat(T);	
	return N;
}
Mat* facetcenter(convexhull* cvh,int facet_index){
	int dim=cvh->dim;		
	Mat* center=newmat(1,dim,0);
	double c=((double)1)/((double)dim);
	for(int i=1;i<=dim;i++){		
		for(int j=1;j<=dim;j++){
			int p_index=(int)cvh->facets->entries[(facet_index-1)*dim+i-1];
			double x1=cvh->points->entries[(p_index-1)*dim+j-1];
			center->entries[j-1]+=c*x1;
		}
	}
	return center;
}
bool isoutside(convexhull* cvh,int facet_index,int point_index){
	int d=cvh->dim;	
	bool b=false;
	double de=0;	
	Mat* N=normal(cvh,facet_index);	
	for(int k=1;k<=d;k++){
		double x1=cvh->points->entries[(point_index-1)*d+k-1];
		int p_index=(int)cvh->facets->entries[(facet_index-1)*d];
		double x2=cvh->points->entries[(p_index-1)*d+k-1];
		de+=(x1-x2)*N->entries[k-1];
	}	
	if(de<-1e-4){
		b=true;
	}
	freemat(N);
	return b;
}
void init(convexhull* cvh){
	int d=cvh->dim;	
	double c1=((double)1)/((double)(d+1));	
	Mat* temp1=newmat(d+1,d,0);
	Mat* center=newmat(1,d,0);
	for(int i=1;i<=d+1;i++){	
		int k=0;
		for(int j=1;j<=d+1;j++){
			center->entries[j-1]+=c1*cvh->points->entries[(i-1)*d+j-1];
			if(i==j){
				continue;
			}
			temp1->entries[(i-1)*d+k]=j;k++;
		}				
	}	
	cvh->facets=temp1;	
	cvh->neighbors_indices=copyvalue(temp1);
	cvh->center=center;
	int Np=cvh->points->row;
	int Nf=cvh->facets->row;
	Mat* Op=newmat(Nf,Np,0);	
	for(int i=1;i<=Nf;i++){							
		int k=0;
		for(int j=1;j<=Np;j++){										
			if(isoutside(cvh,i,j)){
				Op->entries[(i-1)*Np+k]=j;k++;
			}
		}			
	}	
	cvh->outpoints_indices=Op;
}
int furthestpoint(convexhull* cvh,int facet_index){
	int dim=cvh->dim;
	int Np=cvh->points->row;
	Mat* Ope=newmat(1,Np,0);
	Mat* centerf=facetcenter(cvh,facet_index);
	Mat* N=normal(cvh,facet_index);
	getrow2(cvh->outpoints_indices,Ope,facet_index);
	double d=0;
	int ind=0;
	for(int i=1;i<=Np;i++){
		if(Ope->entries[i-1]==0){
			break;
		}
		int op=(int)Ope->entries[i-1];
		double de=0;
		for(int j=1;j<=dim;j++){
			double x=cvh->points->entries[(op-1)*dim+j-1]-centerf->entries[j-1];
			x=x*N->entries[j-1];
			de+=x;
		}
		de=abs(de);
		if(ind==0||de>d){
			ind=op;
			d=de;
		}
	}
	freemat(centerf);
	freemat(Ope);
	freemat(N);
	return ind;
}
void vectoradd(Mat* A,double x){
	double* temp=(double*)realloc(A->entries,sizeof(double)*(A->col+1));
	if(temp==NULL){exit(1);}
	A->entries=temp;
	A->entries[A->col]=x;
	A->col++;
}
void localsearch(convexhull* cvh,Mat* local_facets,Mat* boundary_facets,int facet_index,int prev_facet_index,int point_index){
	int dim=cvh->dim;
	bool repeated=false;
	for(int i=1;i<=local_facets->col;i++){
		if (local_facets->entries[i-1]==facet_index){
			repeated=true;
			break;
		}
	}	
	if(repeated){
		return;
	}
	if(!isoutside(cvh,facet_index,point_index)){
		vectoradd(boundary_facets,(double)facet_index);
		int Nfb=cvh->facets->row;		
		double* temp=(double*)realloc(cvh->facets->entries,(Nfb+1)*dim*sizeof(double));
		if(temp==NULL){exit(1);}
		cvh->facets->entries=temp;										
		cvh->facets->row++;
		int k=0;
		for(int i=1;i<=dim;i++){			
			int v1=(int)cvh->facets->entries[(facet_index-1)*dim+i-1];
			for(int j=1;j<=dim;j++){
				int v2=(int)cvh->facets->entries[(prev_facet_index-1)*dim+j-1];
				if (v1==v2){
					cvh->facets->entries[Nfb*dim+k]=v1;k++;
					break;
				}
			}			
		}			
		cvh->facets->entries[Nfb*dim+k]=point_index;				
		return;
	}
	vectoradd(local_facets,(double)facet_index);
	for(int i=1;i<=dim;i++){
		int nf=(int)cvh->neighbors_indices->entries[(facet_index-1)*dim+i-1];	
		localsearch(cvh,local_facets,boundary_facets,nf,facet_index,point_index);
	}
}

void deleterows(Mat* A,Mat* rows){
	double* temp=(double*)malloc(sizeof(double)*(A->row-rows->col)*A->col);
	int k=0;
	for(int i=1;i<=A->row;i++){
		bool remove=false;
		for(int j=1;j<=rows->col;j++){
			if(i==rows->entries[j-1]){
				remove=true;
				break;
			}
		}
		if(!remove){
			for(int j=1;j<=A->col;j++){
				temp[k++]=A->entries[(i-1)*A->col+j-1];
			}
		}
	}
	free(A->entries);
	A->entries=temp;
	A->row-=rows->col;
}

void deletefacets(convexhull* cvh,Mat* facets){
	int dim=cvh->dim;
	deleterows(cvh->facets,facets);
	deleterows(cvh->neighbors_indices,facets);
	deleterows(cvh->outpoints_indices,facets);
	for(int i=1;i<=cvh->neighbors_indices->row;i++){
		for(int j=1;j<=cvh->neighbors_indices->col;j++){
			int t=(int)cvh->neighbors_indices->entries[(i-1)*dim+j-1];
			int nd=0;
			for(int k=1;k<=facets->col;k++){
				if(t>facets->entries[k-1]){
					nd++;
				}
			}
			cvh->neighbors_indices->entries[(i-1)*dim+j-1]=t-nd;
		}
	}
}

void updatenewfacets(convexhull* cvh,Mat* local_facets,Mat* boundary_facets,int Nfb){
	int dim=cvh->dim;
	int m=cvh->facets->row-Nfb;
	double* temp=(double*)realloc(cvh->neighbors_indices->entries,sizeof(double)*cvh->facets->row*dim);
	if(temp==NULL){exit(1);}
	cvh->neighbors_indices->entries=temp;
	for(int i=1;i<=m;i++){		
		for(int j=1;j<=dim-1;j++){
			int t=(int)cvh->facets->entries[(Nfb+i-1)*dim+j-1];
			for(int k=1;k<=m;k++){
				if(k==i){
					continue;
				}	
				int nc=0;	
				for(int l1=1;l1<=dim;l1++){
					if(l1==j){
						continue;
					}
					int t1=(int)cvh->facets->entries[(Nfb+i-1)*dim+l1-1];
					for(int l2=1;l2<=dim;l2++){
						int t2=(int)cvh->facets->entries[(Nfb+k-1)*dim+l2-1];
						if(t1==t2){
							nc++;
							break;
						}
					}
				}
				if(nc==(dim-1)){
					cvh->neighbors_indices->entries[(Nfb+i-1)*dim+j-1]=Nfb+k;					
					break;
				}
			}

		}		
		for(int j=1;j<=boundary_facets->col;j++){					
			int nc=0;int mc=0;	
			for(int k=1;k<=dim;k++){				
				int t1=(int)cvh->facets->entries[(int)(boundary_facets->entries[j-1]-1)*dim+k-1];
				bool changed=false;
				for(int l=1;l<=dim-1;l++){
					int t2=(int)cvh->facets->entries[(Nfb+i-1)*dim+l-1];					
					if(t1==t2){
						nc++;
						changed=true;
						break;
					}
				}
				if(!changed){
					mc=k;
				}
			}
			if(nc==(dim-1)){
				cvh->neighbors_indices->entries[(Nfb+i-1)*dim+dim-1]=boundary_facets->entries[j-1];
				cvh->neighbors_indices->entries[(int)(boundary_facets->entries[j-1]-1)*dim+mc-1]=Nfb+i;
			}
		}		
	}	
	cvh->neighbors_indices->row+=m;	
	int Np=cvh->points->row;
	double* temp2=(double*)realloc(cvh->outpoints_indices->entries,sizeof(double)*cvh->facets->row*Np);
	if(temp2==NULL){exit(1);}
	cvh->outpoints_indices->entries=temp2;
	cvh->outpoints_indices->row+=m;	
	for(int i=1;i<=m;i++){
		for(int j=1;j<=Np;j++){
			cvh->outpoints_indices->entries[(Nfb+i-1)*Np+j-1]=0;
		}
	}
	for(int i=1;i<=m;i++){
		int ka=0;
		for(int j=1;j<=local_facets->col;j++){
			for(int k=1;k<=Np;k++){
				int op=(int)cvh->outpoints_indices->entries[(int)(local_facets->entries[j-1]-1)*Np+k-1];
				if(op==0){
					break;
				}
				bool repeated=false;
				for(int l=1;l<=Np;l++){
					if(cvh->outpoints_indices->entries[(Nfb+i-1)*Np+l-1]==0){
						break;
					}
					if(op==cvh->outpoints_indices->entries[(Nfb+i-1)*Np+l-1]){
						repeated=true;
						break;
					}
				}
				if(!repeated&&isoutside(cvh,Nfb+i,op)){
					cvh->outpoints_indices->entries[(Nfb+i-1)*Np+ka]=op;
					ka++;
				}
			}
		}
	}		
	deletefacets(cvh,local_facets);	
}

convexhull* convh(Mat* points){
	int dim=points->col;
	int Np=points->row;
	convexhull* cvh=(convexhull*)malloc(sizeof(convexhull));
	cvh->points=points;cvh->dim=dim;	
	init(cvh);	
	while(1==1){				
		int i=1;
		while(i<=cvh->facets->row&&cvh->outpoints_indices->entries[(i-1)*Np]==0){
			i++;
			continue;
		}		
		if(i>cvh->facets->row){
			break;
		}
		int fp=furthestpoint(cvh,i);				
		Mat* local_facets=newmat(1,0,0);
		Mat* boundary_facets=newmat(1,0,0);
		int Nf1=cvh->facets->row;
		localsearch(cvh,local_facets,boundary_facets,i,0,fp);		
		updatenewfacets(cvh,local_facets,boundary_facets,Nf1);									
		freemat(local_facets);
		freemat(boundary_facets);		
	}
	return cvh;
}
