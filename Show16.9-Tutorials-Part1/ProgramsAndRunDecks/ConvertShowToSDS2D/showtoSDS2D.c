#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <df.h>

/*

	gcc -I/usr/local/hdf/include -O3 -Wall -o showtoSDS2D showtoSDS2D.c -L/usr/local/hdf/lib -ldf -ljpeg -lz

*/

#ifndef FALSE
#define FALSE 0
#endif


char WarningBuff[256];

#define min(a,b)         ((a)<(b)?(a):(b))

#define max(a,b)         ((a)>(b)?(a):(b))

#define uint16 unsigned int
#define DFTAG_SDG   ((uint16)700) /* Scientific Data Group*/

FILE *nout,*in;


int DFerror;

struct SDS2Dout{
    char *path;
    char *name;
    char *pioName;
    long ixmax;
    long iymax;
    double *data;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    double vmin;
    double vmax;
    double dt;
    int n;
};
 
struct SHOW{
    long nodeCount;
    long valueCount;
    long normalCount;
    long dxyzCount;
    long rsCount;
    long localSystem;
    long elementCount;
    long stressCount;
    long nodesMax;
    long offset;
    double *x;
    double *y;
    double *z;
    double *v;
    double *nx;
    double *ny;
    double *nz;
    double *dx;
    double *dy;
    double *dz;
    double *r;
    double *s;
    double *stress;
    long *elements;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    double vmin;
    double vmax;
};

char fin[256];

int writesds2d(struct SDS2Dout *sdsout);

int getFloatList(FILE *in,double *FloatData,unsigned long count);
int getLongList(FILE *in,long *LongData,unsigned long count);
int getString(unsigned char *data,long length,FILE *out);
int cFree(char *ptr);
char *cRealloc(char *p,unsigned long r,int tag);
char *cMalloc(unsigned long length,int tag);
int getLong(long *n,FILE *in);
int zerol(char *s,unsigned long n);
int FreeBlock(char *name,struct SHOW *s);
int LoadBlock(char *name,struct SHOW *s);
int SplitBlock(char *name,struct SHOW *s);
int WriteBlock(char *name,struct SHOW *s,FILE *out,long *elements,long ne);
int showWrite(FILE *out,unsigned char *data,unsigned long length);
int showFloat(FILE *out,float *x,unsigned long nn);
int showLong(FILE *out,long *node,unsigned long nn);

int WriteCellCenteredStresses(char *name,struct SHOW *s,FILE *out,long *elements,long ne);
int WriteNodeCenteredValues(char *name,struct SHOW *s,FILE *out,long *elements,long ne);

double dt,scale;

char *nameFile;

int main(int argc,char *argv[])
{
	struct SHOW s;
	unsigned long start,total;
	time_t ship;
	int n;

	
		scale=1.0;
		
        n=1;
        start=time(&ship);
        for(n=1;n<argc;++n){
            if(!strcmp(argv[n],"-dt")){
                dt=atof(argv[++n]);
                continue;
            }else if(!strcmp(argv[n],"-name")){
                nameFile=argv[++n];
                continue;
            }else if(!strcmp(argv[n],"-scale")){
               	scale=atof(argv[++n]);
               	printf("scale = %g\n",scale);
               	continue;
            }
            zerol((char *)&s,sizeof(struct SHOW));
            if(LoadBlock(argv[n],&s))continue;

	}
        total=time(&ship)-start;

        printf("Time To Write Split Files  %ld Seconds\n",total);
        return 0;
}
int SplitBlock(char *name,struct SHOW *s)
{
	long mat[21],nmat,nskip,n,*elements,m,ne,k;
	char buff[256],buff2[256],*bp;
	FILE *out;

	if(!name || !s)return 1;

	out = NULL;

	for(m=0;m<21;++m){
	    mat[m] = 0;
	}

	nskip=s->nodesMax+1;
	for(n=0;n<s->elementCount;++n){
	   nmat=s->elements[n*nskip];
	   if(nmat >= 0 && nmat <= 20)++mat[nmat];
	}
	nmat=0;
	for(m=0;m<21;++m){
	    if(mat[m]>0){
	       printf("File %s mat %ld num %ld\n",name,m,mat[m]);
	       ++nmat;
	    }
	}

	if(nmat == 0){
	   printf("File %s No Materials\n",name);
	   return 0;
	}

	elements=(long *)cMalloc(s->elementCount*nskip*sizeof(long),3500);
	if(!elements){
	   printf("File %s Error Out Of Memory\n",name);
	   return 0;
	}

	for(m=0;m<21;++m){
	    if(mat[m]>0){
			printf("File %s mat %ld num %ld\n",name,m,mat[m]);
			strcpy(buff2,name);
			bp=buff2;
			for(k=0;k<strlen(name);++k){
			   if(buff2[k] != ' ')*bp++ = buff2[k];
			}
			*bp++ = '\0';
			sprintf(buff,"%s.%04ld.s2d",buff2,m);
			out = NULL;
			ne=0;
			for(n=0;n<s->elementCount;++n){
			nmat=s->elements[n*nskip];
			if(nmat == m){
				for(k=0;k<nskip;++k){
					elements[ne*nskip+k]=s->elements[n*nskip+k];
				}
				++ne;
			}
			}
			if(ne != mat[m]){
				printf("File %s Error ne %ld mat %ld\n",name,ne,mat[m]);
			}
			WriteBlock(buff,s,out,elements,ne);
			out = NULL;

	    }
	}

	if(elements)cFree((char *)elements);

	return 0;
}
int WriteBlock(char *name,struct SHOW *s,FILE *out,long *elements,long ne)
{
	double xmin,xmax,ymin,ymax;
	double x,y;
	long n;
	
	
	if(!name || !s || !elements)return 1;

	if(s->nodeCount){
	    fprintf(stderr,"nodeCount %ld ",s->nodeCount);
	}

	if(s->valueCount){
	    fprintf(stderr,"valueCount %ld ",s->valueCount);
	}

	if(s->normalCount){
	    fprintf(stderr,"normalCount %ld ",s->normalCount);
	}

	if(s->dxyzCount){
	    fprintf(stderr,"dxyzCount %ld \n",s->dxyzCount);

	    xmin = 1e30;
	    xmax =  -1e30;
	    ymin = 1e30;
	    ymax = -1e30;
	    for(n=0;n<s->dxyzCount;++n){
	        x=scale*s->dx[n];
	        y=scale*s->dy[n];
	        if(x < xmin)xmin=x;
	        if(x > xmax)xmax=x;
	        if(y < ymin)ymin=y;
	        if(y > ymax)ymax=y;
	    }
	    fprintf(stderr,"dxmin %f dxmax %f\n",xmin,xmax);
	    fprintf(stderr,"dymin %f dymax %f\n",ymin,ymax);


	}

	if(s->rsCount){
	    fprintf(stderr,"rsCount %ld ",s->rsCount);
	}

	if(ne){
	    fprintf(stderr,"ne %ld ",ne);
	}

	if(s->stressCount){
	    fprintf(stderr,"stressCount %ld ",s->stressCount);
	}

	fprintf(stderr,"\n");

	if((s->nodeCount > 0) && (ne > 0) && (s->stressCount == ne)){
	    fprintf(stderr,"WriteCellCenteredStresses nodesMax %ld\n",(long)s->nodesMax);
	    if(WriteCellCenteredStresses(name,s,out,elements,ne))goto ErrorOut;
	}

	if((s->nodeCount > 0) && (ne > 0) && (s->valueCount == s->nodeCount)){
		fprintf(stderr,"WriteNodeCenteredValues nodesMax %ld\n",(long)s->nodesMax);
	    if(WriteNodeCenteredValues(name,s,out,elements,ne))goto ErrorOut;
	}

	return 0;
ErrorOut:
	return 1;
}
int WriteNodeCenteredValues(char *name,struct SHOW *s,FILE *out,long *elements,long ne)
{
	double xmin,xmax,ymin,ymax,vmin,vmax;
	static char lastname[256];
	static int count;
	double x,y,v;
	double didx,djdy,dxdi,dydj;
	long ixmin,iymin,ixmax,iymax,length;
	double *data;
	long n,np,nn,nc;
	long n1,n2,n3,n4;
	double x1,x2,x3,x4;
	double y1,y2,y3,y4;
	double v1,v2,v3,v4;
	double dx1,dx2,dy1,dy2;
	int i;
	struct SDS2Dout sdsout;
	double stress;
	double d,od,s1,s2;
	
	
	if(!name || !s || !elements)return 1;
	
	
	if(s->nodesMax != 3 && s->nodesMax != 4){
	    fprintf(stderr,"WriteNodeCenteredValues Handles only 3 or 4 Noded Elements \n");
	    return 1;
	}
	
	xmin = 1e30;
	xmax =  -1e30;
	ymin = 1e30;
	ymax = -1e30;
	vmin = 1e30;
	vmax = -1e30;
	for(n=0;n<ne;++n){
		np=n*(s->nodesMax+1)+1;
		for(i=0;i<s->nodesMax;++i){
		    nn=elements[np+i];
		    if(nn < 0 || nn >= s->nodeCount)continue;
			x=s->x[nn];
			y=s->y[nn];
			if(s->dxyzCount){
			  x += scale*s->dx[nn];
			  y += scale*s->dy[nn];
			}
			if(x < xmin)xmin=x;
			if(x > xmax)xmax=x;
			if(y < ymin)ymin=y;
			if(y > ymax)ymax=y;
			v=s->v[nn];
			if(v < vmin)vmin=v;
			if(v > vmax)vmax=v;
		}
	}
	s->xmin=xmin;
	s->ymin=ymin;
	s->vmin=vmin;
	s->zmin=0;
	s->xmax=xmax;
	s->ymax=ymax;
	s->vmax=vmax;
	s->zmax=0;
	fprintf(stderr,"xmin %f xmax %f\n",xmin,xmax);
	fprintf(stderr,"ymin %f ymax %f\n",ymin,ymax);
	fprintf(stderr,"vmin %f vmax %f\n",vmin,vmax);


	ixmin=0;
	iymin=0;

	ixmax=601;
	iymax=601;

	length=ixmax*iymax;
	data=(double *)cMalloc(length*sizeof(double),8203);
	if(!data)goto ErrorOut;
	zerol((char *)data,length*sizeof(double));

	didx=(double)(ixmax-ixmin)/(xmax-xmin);
	dxdi= 1.0/didx;
	djdy=(double)(iymax-iymin)/(ymax-ymin);
	dydj= 1.0/djdy;
	
	for(n=0;n<ne;++n){
		double ymn,ymx,xmn,xmx;
		double xs,ys;
		long i,j,ihi,ilo,jhi,jlo;

		xmn = 1e30;
		xmx =  -1e30;
		ymn = 1e30;
		ymx = -1e30;

		
		nc=0;
		stress=0;
		np=n*(s->nodesMax+1)+1;
		
		for(i=0;i<s->nodesMax;++i){
		    nn=elements[np+i];
		    if(nn < 0 || nn >= s->nodeCount)continue;
		    nc++;
			x=s->x[nn];
			y=s->y[nn];
			if(s->dxyzCount){
			  x += scale*s->dx[nn];
			  y += scale*s->dy[nn];
			}
			if(x < xmn)xmn=x;
			if(x > xmx)xmx=x;
			if(y < ymn)ymn=y;
			if(y > ymx)ymx=y;
			stress += s->v[nn];
		}
		
		stress /= (double)(nc);
		
		ymn=max(ymin,ymn);
		ymx=min(ymax,ymx);

		if(ymx < ymn)continue;

		xmn=max(xmin,xmn);
		xmx=min(xmax,xmx);


		if(xmx < xmn)continue;


		ilo=max(ixmin,(long)((xmn-xmin)*didx)+ixmin);
		ihi=min(ixmax-1,(long)((xmx-xmin)*didx)+ixmin);

		jlo=max(iymin,(long)((ymn-ymin)*djdy)+iymin);
		jhi=min(iymax-1,(long)((ymx-ymin)*djdy)+iymin);
/*
		jhi=min(iymax-1,(long)((ymn-ymin)*djdy)+iymax);
		jlo=max(iymin,(long)((ymx-ymin)*djdy)+iymax);
*/		
		if(s->nodesMax == 3 || elements[np+3] < 0){
		    n1=elements[np];
		    n2=elements[np+1];
		    n3=elements[np+2];
		    x1=s->x[n1];
		    if(s->dxyzCount)x1 += scale*s->dx[n1];
		    x2=s->x[n2];
		    if(s->dxyzCount)x2 += scale*s->dx[n2];
		    x3=s->x[n3];
		    if(s->dxyzCount)x3 += scale*s->dx[n3];
		    
		    y1=s->y[n1];
		    if(s->dxyzCount)y1 += scale*s->dy[n1];
		    y2=s->y[n2];
		    if(s->dxyzCount)y2 += scale*s->dy[n2];
		    y3=s->y[n3];
		    if(s->dxyzCount)y3 += scale*s->dy[n3];
		    
		    v1=s->v[n1];
		    v2=s->v[n2];
		    v3=s->v[n3];
		    
		    dx1=x2-x1;
		    dx2=x3-x1;

		    dy1=y2-y1;
		    dy2=y3-y1;
		    
			if((d=(dx1)*(dy2)-(dy1)*(dx2))){
						
				od=1.f/d;
				
				for(j=jlo;j<=jhi;++j){
					ys=ymin+0.5*dydj+(double)(j-iymin)*dydj;      /* locate y center of pixel */
					if(ys < ymn || ys > ymx)continue;
					for(i=ilo;i<=ihi;++i){
						xs=xmin+0.5*dxdi+(double)(i-ixmin)*dxdi;  /* locate x center of pixel */
						if(xs < xmn || xs > xmx)continue;
						
						/* (xs,ys) are the real coordinates of the raster points (i,j) */
						
						s1=((xs-x1)*dy2-(ys-y1)*dx2)*od;
						if(s1 < 0 || s1 > 1.0)continue;
						s2=(dx1*(ys-y1)-dy1*(xs-x1))*od;
						if(s2 < 0 || s2 > 1.0)continue;
						
						if((s1+s2) > 1.0)continue;
						
						data[i+ixmax*(iymax-j-1)]=v1+s1*(v2-v1)+s2*(v3-v1);
					}
				}
			}
		}else{
		    n1=elements[np];
		    n2=elements[np+1];
		    n3=elements[np+2];
		    n4=elements[np+3];
		    x1=s->x[n1];
		    if(s->dxyzCount)x1 += scale*s->dx[n1];
		    x2=s->x[n2];
		    if(s->dxyzCount)x2 += scale*s->dx[n2];
		    x3=s->x[n3];
		    if(s->dxyzCount)x3 += scale*s->dx[n3];
		    x4=s->x[n4];
		    if(s->dxyzCount)x4 += scale*s->dx[n4];
		    
		    y1=s->y[n1];
		    if(s->dxyzCount)y1 += scale*s->dy[n1];
		    y2=s->y[n2];
		    if(s->dxyzCount)y2 += scale*s->dy[n2];
		    y3=s->y[n3];
		    if(s->dxyzCount)y3 += scale*s->dy[n3];
		    y4=s->y[n4];
		    if(s->dxyzCount)y4 += scale*s->dy[n4];
		    
		    v1=s->v[n1];
		    v2=s->v[n2];
		    v3=s->v[n3];
		    v4=s->v[n4];
		    
		    dx1=x2-x1;
		    dx2=x3-x1;

		    dy1=y2-y1;
		    dy2=y3-y1;
		    
			if((d=(dx1)*(dy2)-(dy1)*(dx2))){
						
				od=1.f/d;
				
				for(j=jlo;j<=jhi;++j){
					ys=ymin+0.5*dydj+(double)(j-iymin)*dydj;     /* locate y center of pixel */
					if(ys < ymn || ys > ymx)continue;
					for(i=ilo;i<=ihi;++i){
						xs=xmin+0.5*dxdi+(double)(i-ixmin)*dxdi;  /* locate x center of pixel */
						if(xs < xmn || xs > xmx)continue;
						
						/* (xs,ys) are the real coordinates of the raster points (i,j) */
						
						s1=((xs-x1)*dy2-(ys-y1)*dx2)*od;
						if(s1 < 0 || s1 > 1.0)continue;
						s2=(dx1*(ys-y1)-dy1*(xs-x1))*od;
						if(s2 < 0 || s2 > 1.0)continue;
						
						if((s1+s2) > 1.0)continue;
						
						data[i+ixmax*(iymax-j-1)]=v1+s1*(v2-v1)+s2*(v3-v1);
					}
				}
			}
			
		    dx1=x3-x1;
		    dx2=x4-x1;

		    dy1=y3-y1;
		    dy2=y4-y1;
		    
			if((d=(dx1)*(dy2)-(dy1)*(dx2))){
						
				od=1.f/d;
				
				for(j=jlo;j<=jhi;++j){
					ys=ymin+0.5*dydj+(double)(j-iymin)*dydj;     /* locate y center of pixel */
					if(ys < ymn || ys > ymx)continue;
					for(i=ilo;i<=ihi;++i){
						xs=xmin+0.5*dxdi+(double)(i-ixmin)*dxdi;  /* locate x center of pixel */
						if(xs < xmn || xs > xmx)continue;
						
						/* (xs,ys) are the real coordinates of the raster points (i,j) */
						
						s1=((xs-x1)*dy2-(ys-y1)*dx2)*od;
						if(s1 < 0 || s1 > 1.0)continue;
						s2=(dx1*(ys-y1)-dy1*(xs-x1))*od;
						if(s2 < 0 || s2 > 1.0)continue;
						
						if((s1+s2) > 1.0)continue;
						
						data[i+ixmax*(iymax-j-1)]=v1+s1*(v3-v1)+s2*(v4-v1);
					}
				}
			}
			
			
			
		}
	}


	if(!strcmp(lastname,name)){
		count++;
	}else{
		strncpy(lastname,name,256);
		count=0;
	}

	if(nameFile){
		char buff[256];
		strncpy(buff,nameFile,256);
		strncat(buff,".s2d",256);
		sdsout.path=buff;
		sdsout.name=name;
	}else{
		sdsout.path=name;
		sdsout.name=name;
	}
	sdsout.ixmax=ixmax;
	sdsout.iymax=iymax;
	sdsout.data=data;
	sdsout.xmin=xmin;
	sdsout.xmax=xmax;
	sdsout.ymin=ymin;
	sdsout.ymax=ymax;
	sdsout.zmin=0;
	sdsout.zmax=0;
	sdsout.dt=dt;
	sdsout.n=count;
	sdsout.pioName=nameFile;

	writesds2d(&sdsout);

	return 0;
ErrorOut:
	return 1;
	

}
int WriteCellCenteredStresses(char *name,struct SHOW *s,FILE *out,long *elements,long ne)
{
	double xmin,xmax,ymin,ymax,vmin,vmax;
	static char lastname[256];
	static int count;
	double x,y,v;
	double didx,djdy,dxdi,dydj;
	long ixmin,iymin,ixmax,iymax,length;
	double *data;
	long n,np;
	int i;
	struct SDS2Dout sdsout;
	
	
	if(!name || !s || !elements)return 1;

		data=NULL;
	    xmin = 1e30;
	    xmax =  -1e30;
	    ymin = 1e30;
	    ymax = -1e30;
	    vmin = 1e30;
	    vmax = -1e30;
	    for(n=0;n<ne;++n){
	        np=n*(s->nodesMax+1)+1;
	        for(i=0;i<s->nodesMax;++i){
	            x=s->x[elements[np+i]];
	            y=s->y[elements[np+i]];
				if(s->dxyzCount){
				  x += scale*s->dx[elements[np+i]];
				  y += scale*s->dy[elements[np+i]];
				}
	            if(x < xmin)xmin=x;
	            if(x > xmax)xmax=x;
	            if(y < ymin)ymin=y;
	            if(y > ymax)ymax=y;
				v=s->v[elements[np+i]];
				if(v < vmin)vmin=v;
				if(v > vmax)vmax=v;
	        }
	    }
	    s->xmin=xmin;
	    s->ymin=ymin;
	    s->vmin=vmin;
	    s->zmin=0;
	    s->xmax=xmax;
	    s->ymax=ymax;
	    s->vmax=vmax;
	    s->zmax=0;
	    printf("xmin %f xmax %f\n",xmin,xmax);
	    printf("ymin %f ymax %f\n",ymin,ymax);
	    printf("vmin %f vmax %f\n",vmin,vmax);


	    ixmin=0;
	    iymin=0;

	    ixmax=401;
	    iymax=501;

	    ixmax=601;
	    iymax=601;

	    length=ixmax*iymax;
	    data=(double *)cMalloc(length*sizeof(double),8203);
	    if(!data)goto ErrorOut;
	    zerol((char *)data,length*sizeof(double));

	    didx=(double)(ixmax-ixmin)/(xmax-xmin);
	    dxdi= 1.0/didx;
	    djdy=(double)(iymin-iymax)/(ymax-ymin);
	    dydj= 1.0/djdy;
	    for(n=0;n<ne;++n){
                double ymn,ymx,xmn,xmx;
                double xs,ys;
                long i,j,ihi,ilo,jhi,jlo;

				xmn = 1e30;
				xmx =  -1e30;
				ymn = 1e30;
				ymx = -1e30;
	
				np=n*(s->nodesMax+1)+1;
	
				for(i=0;i<s->nodesMax;++i){
					x=s->x[elements[np+i]];
					y=s->y[elements[np+i]];
					if(s->dxyzCount){
					  x += scale*s->dx[elements[np+i]];
					  y += scale*s->dy[elements[np+i]];
					}
					if(x < xmn)xmn=x;
					if(x > xmx)xmx=x;
					if(y < ymn)ymn=y;
					if(y > ymx)ymx=y;
				}
				
                ymn=max(ymin,ymn);
                ymx=min(ymax,ymx);

                if(ymx < ymn)continue;

                xmn=max(xmin,xmn);
                xmx=min(xmax,xmx);


                if(xmx < xmn)continue;


                ilo=max(ixmin,(long)((xmn-xmin)*didx)+ixmin);
                ihi=min(ixmax-1,(long)((xmx-xmin)*didx)+ixmin);

                jhi=min(iymax-1,(long)((ymn-ymin)*djdy)+iymax);
                jlo=max(iymin,(long)((ymx-ymin)*djdy)+iymax);


	        v=s->stress[n];
/*
                if(n < 24){
                    printf("n %ld ne %ld nodesMax %ld %g\n",(long)n,(long)ne,(long)s->nodesMax,v);
                    printf("%ld %ld %ld %ld\n",(long)elements[np+0],(long)elements[np+1],(long)elements[np+2],(long)elements[np+3]);
                    printf("jlo %ld jhi  %ld ilo  %ld ihi  %ld \n",jlo,jhi,ilo,ihi);
                    printf("ymn %g ymx %g xmn %g xmx %g\n",ymn,ymx,xmn,xmx);
                    printf("ymin %g xmin %g\n",ymin,xmin);
                    printf("djdy %g didx %g\n",djdy,didx);
                    printf("\n");
                }
*/
	        v=s->stress[n];
                for(j=jlo;j<=jhi;++j){
                    ys=ymin+(double)(j-iymax)*dydj;
                    if(ys < ymn || ys > ymx)continue;
                    for(i=ilo;i<=ihi;++i){
                        xs=xmin+(double)(i-ixmin)*dxdi;
                        if(xs < xmn || xs > xmx)continue;
                        data[i+j*ixmax]=v;
                    }
                }
	    }


	    if(!strcmp(lastname,name)){
	        count++;
	    }else{
	        strncpy(lastname,name,256);
	        count=0;
	    }

	    if(nameFile){
	        char buff[256];
	        strncpy(buff,nameFile,256);
	        strncat(buff,".s2d",256);
	        sdsout.path=buff;
	        sdsout.name=name;
	    }else{
	        sdsout.path=name;
	        sdsout.name=name;
	    }
	    sdsout.ixmax=ixmax;
	    sdsout.iymax=iymax;
	    sdsout.data=data;
	    sdsout.xmin=xmin;
	    sdsout.xmax=xmax;
	    sdsout.ymin=ymin;
	    sdsout.ymax=ymax;
	    sdsout.zmin=0;
	    sdsout.zmax=0;
	    sdsout.dt=dt;
	    sdsout.n=count;
	    sdsout.pioName=nameFile;

	    writesds2d(&sdsout);
	    

	return 0;
ErrorOut:
	return 1;

}
int LoadBlock(char *name,struct SHOW *s)
{
        char buff[256];
	FILE *in;
	int ok;

	if(!name || !s)return 1;

	    ok=FALSE;

	    in=fopen(name,"rb");
	    if(in == NULL){
	       fprintf(stderr,"Could Not Open %s To Read\n",name);
	       goto Loop;
	    }

	    if(getString((unsigned char *)buff,4L,in))goto Loop;

	    buff[4]=0;

	    if(!strcmp(buff,"SHOW")){
	        ;
	    }else{
	        fprintf(stderr,"File %s Skipped Not A Show File (%s)\n",
	                name,buff);
	        goto Loop;
	    }


	    while(1){
	        if(getString((unsigned char *)buff,4L,in))goto Loop;

	        buff[4]=0;

	        if(!strcmp(buff,"BLK1")){
	            printf("BLK1\n");
	        }else if(!strcmp(buff,"END1")){
		    	if(SplitBlock(name,s))goto Loop;
	            if(FreeBlock(name,s))goto Loop;
	            printf("END1\n");
	            continue;
	        }else if(!strcmp(buff,"DONE")){
	            printf("DONE\n");
	            break;
	        }else{
	            fprintf(stderr,"File %s Skipped Not A BLK1 record (%s)\n",
	                name,buff);
	            goto Loop;
	        }


		if(getLongList(in,&s->nodeCount,1L))goto Loop;
		if(getLongList(in,&s->valueCount,1L))goto Loop;
		if(getLongList(in,&s->normalCount,1L))goto Loop;
		if(getLongList(in,&s->dxyzCount,1L))goto Loop;
		if(getLongList(in,&s->rsCount,1L))goto Loop;
		if(getLongList(in,&s->elementCount,1L))goto Loop;
		if(getLongList(in,&s->nodesMax,1L))goto Loop;
		if(getLongList(in,&s->localSystem,1L))goto Loop;
		if(getLongList(in,&s->stressCount,1L))goto Loop;


		    if(s->nodeCount){
				s->x=(double *)cMalloc(s->nodeCount*sizeof(double),8454);
				s->y=(double *)cMalloc(s->nodeCount*sizeof(double),8455);
				s->z=(double *)cMalloc(s->nodeCount*sizeof(double),8456);
				if(!s->x || !s->y || !s->z)goto Loop;
		        if(getFloatList(in,s->x,s->nodeCount))goto Loop;
		        if(getFloatList(in,s->y,s->nodeCount))goto Loop;
		        if(getFloatList(in,s->z,s->nodeCount))goto Loop;
		    }

		    if(s->valueCount){
				s->v=(double *)cMalloc(s->valueCount*sizeof(double),8457);
				if(!s->v )goto Loop;
		        if(getFloatList(in,s->v,s->valueCount))goto Loop;
		    }

		    if(s->normalCount){
				s->nx=(double *)cMalloc(s->normalCount*sizeof(double),8458);
				s->ny=(double *)cMalloc(s->normalCount*sizeof(double),8459);
				s->nz=(double *)cMalloc(s->normalCount*sizeof(double),8460);
				if(!s->nx || !s->ny || !s->nz)goto Loop;
		        if(getFloatList(in,s->nx,s->normalCount))goto Loop;
		        if(getFloatList(in,s->ny,s->normalCount))goto Loop;
		        if(getFloatList(in,s->nz,s->normalCount))goto Loop;
		    }

		    if(s->dxyzCount){
				s->dx=(double *)cMalloc(s->dxyzCount*sizeof(double),8461);
				s->dy=(double *)cMalloc(s->dxyzCount*sizeof(double),8462);
				s->dz=(double *)cMalloc(s->dxyzCount*sizeof(double),8463);
				if(!s->dx || !s->dy || !s->dz)goto Loop;
		        if(getFloatList(in,s->dx,s->dxyzCount))goto Loop;
		        if(getFloatList(in,s->dy,s->dxyzCount))goto Loop;
		        if(getFloatList(in,s->dz,s->dxyzCount))goto Loop;
		    }

		    if(s->rsCount){
				s->r=(double *)cMalloc(s->rsCount*sizeof(double),8464);
				s->s=(double *)cMalloc(s->rsCount*sizeof(double),8465);
				if(!s->r || !s->s)goto Loop;
		        if(getFloatList(in,s->r,s->rsCount))goto Loop;
		        if(getFloatList(in,s->s,s->rsCount))goto Loop;
		    }

		    if(s->elementCount){
				s->elements=(long *)cMalloc(s->elementCount*(s->nodesMax+1)*sizeof(long),8466);
				if(!s->elements)goto Loop;
		        if(getLongList(in,s->elements,s->elementCount*(s->nodesMax+1)))goto Loop;
		    }
		    if(s->stressCount){
				s->stress=(double *)cMalloc(s->stressCount*sizeof(double),8467);
				if(!s->stress)goto Loop;
		        if(getFloatList(in,s->stress,s->stressCount))goto Loop;
		    }
	    }
	    ok=1;
Loop:
	if(in)fclose(in); 
	in=NULL;
	if(ok)return 0;
	fprintf(stderr,"File %s Error\n",name);
	return 1;
}
int getFloatList(FILE *in,double *FloatData,unsigned long count)
{
	union test{
		unsigned long int n;
		unsigned char cn[4];
	}t;
	unsigned long length;
	int isPC;
	


	if(!in || !FloatData)return 1;

	length=count*sizeof(double);

	if(fread(FloatData,1,length,in) != length)return 1;
	
	t.n=1;
	isPC=t.cn[0];
	
	if(isPC){
		unsigned char c,*cp;
		long np;

	        for(np=0;np<count;np++){
	            cp=(unsigned char *)&FloatData[np];
				c=cp[0];
				cp[0]=cp[7];
				cp[7]=c;
				c=cp[1];
				cp[1]=cp[6];
				cp[6]=c;
				c=cp[2];
				cp[2]=cp[5];
				cp[5]=c;
				c=cp[3];
				cp[3]=cp[4];
				cp[4]=c;
	        }
	}

	return 0;
}
int getLongList(FILE *in,long *LongData,unsigned long count)
{
	unsigned long n;

	if(!in || !LongData)return 1;

	for(n=0;n<count;++n){
	   if(getLong(&LongData[n],in))return 1;
	}

	return 0;
}
int getLong(long *n,FILE *in)
{
    unsigned long c1,c2,c3,c4;
    unsigned char c[4];

	if(getString((unsigned char *)c,4L,in))return 1;

	c1=c[0];
	c2=c[1];
	c3=c[2];
	c4=c[3];
	
    *n = (long)((c1 << 24) | (c2 << 16) | (c3 << 8) | c4);

    return 0;
 
}
int getString(unsigned char *data,long length,FILE *out)
{
        if(!data || !out)return 1;
        if(fread(data,1,length,out) != length){
            return 1;
        }
        return 0;
}
char *cMalloc(unsigned long length,int tag)
{
/*      extern char *malloc(long); */
        char *ret;

        tag=tag;

        ret=(char *)malloc(length+4L);
        if(ret == NULL){
            fprintf(stderr,"cMalloc Out of Memory Requested (%ld)",length);
            return (char *)NULL;
        }
        return ret;
}
char *cRealloc(char *p,unsigned long r,int tag)
{
        char *pn;

        tag=tag;

        if(!(pn=(char *)realloc(p,r))){
            return (char *)NULL;
        }

        return pn;
}
int cFree(char *ptr)
{
        if(ptr){
            free(ptr);
            return 0;
        }
        return 1;
}
int showLong(FILE *out,long *node,unsigned long nn)
{

	if(!out || !node || !nn)return 1;

        if(showWrite(out,(unsigned char *)node,nn*sizeof(long)))return 1;

	return 0;
}
int showFloat(FILE *out,float *x,unsigned long nn)
{
	if(!out || !x || !nn)return 1;

        if(showWrite(out,(unsigned char *)x,nn*sizeof(float)))return 1;

	return 0;
}
int showWrite(FILE *out,unsigned char *data,unsigned long length)
{
	if(!out || !data || !length)return 1;

        if(fwrite(data,1,length,out) != length){
            return 1;
        }
	
	return 0;
}
int FreeBlock(char *name,struct SHOW *s)
{
	if(!name || !s)return 1;
	if(s->x)cFree((char *)s->x);
	s->x=NULL;
	if(s->y)cFree((char *)s->y);
	s->y=NULL;
	if(s->z)cFree((char *)s->z);
	s->z=NULL;
	if(s->v)cFree((char *)s->v);
	s->v=NULL;
	if(s->nx)cFree((char *)s->nx);
	s->nx=NULL;
	if(s->ny)cFree((char *)s->ny);
	s->ny=NULL;
	if(s->nz)cFree((char *)s->nz);
	s->nz=NULL;
	if(s->dx)cFree((char *)s->dx);
	s->dx=NULL;
	if(s->dy)cFree((char *)s->dy);
	s->dy=NULL;
	if(s->dz)cFree((char *)s->dz);
	s->dz=NULL;
	if(s->r)cFree((char *)s->r);
	s->r=NULL;
	if(s->s)cFree((char *)s->s);
	s->s=NULL;
	if(s->stress)cFree((char *)s->stress);
	s->stress=NULL;
	if(s->elements)cFree((char *)s->elements);
	s->elements=NULL;
	return 0;
}
int zerol(char *s,unsigned long n)
{
        unsigned long count;

        if(s == NULL)return 1;
        count=0;
        while(count++ < n)*s++ = 0;
        return 0;
}
int Warning(char *mess)
{
    if(!mess)return 1;
    fprintf(stderr,"%s\n",mess);
    return 0;
}
int writesds2d(struct SDS2Dout *sdsout)
{
	extern int DFerror;
	double vmin,vmax;
	char buff[256];
	int32 size[2];
	int lastref;
	double v;
	int rank;
	long n;

	if(!sdsout)return 1;

	rank=2;
	size[0]=(int)sdsout->iymax;
	size[1]=(int)sdsout->ixmax;

	if(sdsout->n == 0){
	    unlink((char *)sdsout->path);
	    DFSDclear();
	    DFSDrestart();
        DFSDsetNT(DFNT_FLOAT64);
        
	    if(DFSDputdata((char *)sdsout->path,rank,size,(void *)sdsout->data)){
	        sprintf(WarningBuff,"writesds2d DFSDputdata error %d",DFerror);
	        Warning(WarningBuff);
	        goto OutOfHere;
	    }
	}else{
	    if(DFSDadddata((char *)sdsout->path,rank,size,(void *)sdsout->data)){
	        sprintf(WarningBuff,"writesds2d DFSDadddata error %d",DFerror);
	        Warning(WarningBuff);
	        goto OutOfHere;
	    }
	}
	
	vmin =  1e33;
	vmax = -1e33;
	for(n=0;n<sdsout->ixmax*sdsout->iymax;++n){
	    v=sdsout->data[n];
	    if(v < vmin)vmin = v;
	    if(v > vmax)vmax = v;
	}

	lastref=DFSDlastref();
	if(lastref == -1){
	    sprintf(WarningBuff,"writesds2d DFSDlastref error %d",DFerror);
	    Warning(WarningBuff);
	    goto OutOfHere;
	}

	if(DFANputlabel(sdsout->path,DFTAG_SDG,lastref,sdsout->name) == -1){
	    sprintf(WarningBuff,"writesds2d DFANputlabel %s Name %s lastref %d error %d",
		               sdsout->path,sdsout->name,lastref,DFerror);
	    Warning(WarningBuff);
	    goto OutOfHere;
	}
	if(sdsout->pioName){
	    sprintf(buff,"xmin %g ymin %g zmin %g xmax %g ymax %g zmax %g vmin %g vmax %g time %g pioName \"%s\" ",
	            sdsout->xmin,sdsout->ymin,sdsout->zmin,
                sdsout->xmax,sdsout->ymax,sdsout->zmax,
                vmin,vmax,sdsout->dt*sdsout->n,sdsout->pioName);
	}else{
	    sprintf(buff,"xmin %g ymin %g zmin %g xmax %g ymax %g zmax %g vmin %g vmax %g time %g",
	            sdsout->xmin,sdsout->ymin,sdsout->zmin,
                sdsout->xmax,sdsout->ymax,sdsout->zmax,
                vmin,vmax,sdsout->dt*sdsout->n);
	}
	if(DFANputdesc(sdsout->path,DFTAG_SDG,lastref,(char *)buff,strlen((char *)buff)) == -1){
	    sprintf(WarningBuff,"DFANputdesc %s Name %s lastref %d DFerror %d",
		sdsout->path,sdsout->name,lastref,DFerror);
	    Warning(WarningBuff);
	    goto OutOfHere;
	}
	return 0;
OutOfHere:
	return 1;
}
