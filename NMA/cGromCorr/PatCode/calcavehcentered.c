#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double allatoms[4500][3]={0}; /*index1=atom number, index2=0,1,2,3,4 for atom number, atom type, x,y,z position*/
double oxygens[1500][3]={0};/*index 1=atom number, index 2=0,1,2, for x, y, z position*/
double densities[300][300][300][4]={0};/*index 1through3=boxindex, index 4=0,1,2,3 for x,y,z position,density*/
double avdensity[300][300][300]={0};
double h[1000][1000][3];/*index 1,2=cellindex, index 2=0,1,2 for x,y, h*/
double averageh[1000][1000][4];
long int Natoms=0;
double xmax=0;
double xmin=0;
double ymax=0;
double ymin=0;
double zmax=0;
double zmin=0;
double celllength=.5;
double s=3.1;
double s2=4.3;
double cutoff=10;
int timesteps=0;
double x,y,z=0;
double boxx, hboxx, boxy, hboxy, boxz, hboxz=0;
double volume;
double Nboxes;
double cellvolume;
double averagebulk=13.524816;
int i,j,k,l=0;
int c=0; /*keeps track of the number of oxygens*/
char trajfile[128];
char hpfile[128];
char havfile[128];
void createoxygenarray();
void discretizebox();
void recenter();
void computedensity();
void normalizedensity1();
void normalizedensity2();
void computeh();
void printgraph(const int);
void printaverage();
void average();

int main(int argc, char *argv[]){
  int t=0;
  int i1=0;
  char line[100];
  sscanf(argv[1], "%d", &Natoms);
  sscanf(argv[2], "%lf", &celllength);
  sscanf(argv[3], "%s", &trajfile);

  FILE *inputfile;
  inputfile=fopen(trajfile,"r");

  while((fgets(line, sizeof(line), inputfile)) !=NULL){
    if(t<20){
      fscanf(inputfile, "%s", &line);
      printf("%s\n", line);
      for (i=0; i<Natoms; i++){
	fscanf(inputfile,"%s %s %s %lf %lf %lf %s %s %s\n", &line, &line, &line, &allatoms[i][0], &allatoms[i][1], &allatoms[i][2], &line, &line, &line);

	// Convert to angstroms
	allatoms[i][0]=allatoms[i][0]*10; 
	allatoms[i][1]=allatoms[i][1]*10;
	allatoms[i][2]=allatoms[i][2]*10;
	if (i < 2 || i > Natoms - 2){
		printf("Atom %d: %f %f %f\n", i, allatoms[i][0], allatoms[i][1], allatoms[i][2]);
	}
      }
      fscanf(inputfile, " %lf %lf %lf\n", &boxx, &boxy, &boxz);
      printf("Box size %d: %f %f %f\n", t, boxx, boxy, boxz);
      t=t+1;
      boxx=boxx*10; boxy=boxy*10; boxz=boxz*10;
      xmax=boxx; ymax=boxy; zmax=boxz;
      printf("Discretizing\n");
      discretizebox();
      printf("Recentering\n");
      recenter();
      printf("Create Oxygening\n");
      createoxygenarray();
      printf("Compute Densitying\n");
      computedensity();
      printf("Averaging\n");
      average();
      timesteps = timesteps+1;
    }
  }

  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      for(k=0; k<=z-1; k++){
	printf("%f %f %f %f\n", densities[i][j][k][0], densities[i][j][k][1], densities[i][j][k][2], avdensity[i][j][k]/timesteps);
      }
    }
  }

  fclose(inputfile);

  return 0;
}


void recenter(){
  double x0, y0=0;
  x0=allatoms[Natoms-1][1];
  y0=allatoms[Natoms-1][2];

  for(i=0; i<=Natoms-1; i++){
    allatoms[i][1]=allatoms[i][1]-x0;
    allatoms[i][2]=allatoms[i][2]-y0;
  }
}


void createoxygenarray(){
  double comz=0;
  c=0;
  for (i=0; i<=Natoms-2; i++){
    if (i%3==0){
      oxygens[c][0]=allatoms[i][0];
      oxygens[c][1]=allatoms[i][1];
      oxygens[c][2]=allatoms[i][2];
      c = c+1;
    }
  }

}

void discretizebox(){

  x = (xmax-xmin)/celllength;
  y = (ymax-ymin)/celllength;
  z = (zmax-zmin)/celllength;
  Nboxes = volume/cellvolume; 
  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      for(k=0; k<=z-1; k++){
	densities[i][j][k][0] = xmin + i*celllength;
	densities[i][j][k][1] = ymin + j*celllength;
	densities[i][j][k][2] = zmin + k*celllength;
      }
    }
  }
}

void computedensity(){
  int m=0;
  double rx=0, ry=0, rz=0, r2=0;
  boxx= xmax-xmin;
  boxy= ymax-ymin;
  boxz= zmax-zmin;
  hboxx=boxx/2;
  hboxy=boxy/2;
  hboxz=boxz/2;

  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      for(k=0; k<=z-1; k++){
	densities[i][j][k][3]=0;
	for(m=0; m<=c-1; m++){
	  rx=oxygens[m][0]-densities[i][j][k][0];
	  if(rx>hboxx)
	    rx = -boxx+rx;
	  if(rx<-hboxx)
	    rx = boxx + rx;
	  ry=oxygens[m][1]-densities[i][j][k][1];
	  if(ry>hboxy)
	    ry = -boxy+ry;
	  if(ry<-hboxy)
	    ry = boxy+ry;
	  rz=oxygens[m][2]-densities[i][j][k][2];
	  if(rz>hboxz)
	    rz = -boxz+rz;
	  if(rz<-hboxz)
	    rz = boxz+rz;
	  r2=rx*rx + ry*ry + rz*rz;
	  if(r2<cutoff*cutoff){
	    densities[i][j][k][3]=densities[i][j][k][3] + exp(-r2/(2*s*s))-exp(-cutoff*cutoff/(2*s*s));
	  }
	}

      }
    }
  }


  

}

void average(){

  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      for(k=0; k<=z-1; k++){
	avdensity[i][j][k] +=densities[i][j][k][3];
      }
    }
  }

}

void normalizedensity1(){
  double m=0;
  int boxcenter=0;
  boxcenter=(int)(-zmin/celllength);
  averagebulk=0;

  for(i=0; i<=x; i++){
    for(j=0; j<=y; j++){
      averagebulk=averagebulk + densities[i][j][boxcenter][3];
      averagebulk=averagebulk + densities[i][j][boxcenter-1][3];
      averagebulk=averagebulk + densities[i][j][boxcenter-2][3];

      m=m+1;
    }
  }
  averagebulk=averagebulk/(3.0*m);

  printf("%f\n", averagebulk);
}

void normalizedensity2(){
  for(i=0; i<=x; i++){
    for(j=0; j<=y; j++){
      for(k=0; k<=z; k++){
	densities[i][j][k][3]=densities[i][j][k][3]/averagebulk;
      }
    }
  }
}


void computeh(){
  double n=0;
  int u=0;
  double slope=0;
  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      for(k=z; k>0; k=k-1){
	u=k;
	n=densities[i][j][k][3];
	if (n > .5)
	  break;
      }
      h[i][j][0]=densities[i][j][u][0];
      h[i][j][1]=densities[i][j][u][1];
      slope=(densities[i][j][u][2]-densities[i][j][u+1][2])/(densities[i][j][u][3]-densities[i][j][u+1][3]);
      //printf("%f %f %f %f\n", densities[i][j][u][2], densities[i][j][u+1][2], densities[i][j][u][3], densities[i][j][u+1][3]);
      h[i][j][2]=slope*.5 - slope*densities[i][j][u+1][3]+densities[i][j][u+1][2];

    }
  }
}

void printgraph(const int ts){
  FILE *hptr;
  char filename[128];
  sprintf(filename, "%s%d.dat", hpfile, ts);
  hptr = fopen(filename,"w");
  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      fprintf(hptr, "%.4f %.4f %.4f\n", h[i][j][0], h[i][j][1], h[i][j][2]);
    }
  }
  fclose(hptr);
}



void printaverage(){
  int i2=0;
  FILE *havptr;
  havptr=fopen(havfile,"w");
  for(i=0; i<=x-1; i++){
    for(j=0; j<=y-1; j++){
      averageh[i][j][2]=averageh[i][j][2]/timesteps;
      averageh[i][j][3]=sqrt(averageh[i][j][3]/timesteps-averageh[i][j][2]*averageh[i][j][2]);

      fprintf (havptr,"%.4f %.4f %.4f %.4f\n", h[i][j][0], h[i][j][1], averageh[i][j][2], averageh[i][j][3]);
    }
  }
  fclose(havptr);
}



