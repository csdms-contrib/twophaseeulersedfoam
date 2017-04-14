#include <stdio.h>
#include <string.h>


int main()
{
FILE *myfile,*mywfile,*myfolders;
float x[15000],y[15000],z[15000],u[15000],v[15000],w[15000],u_mean,v_mean,w_mean,u0,v0,w0;
int i, N;
char c[300],filename[300],filename1[300];

system("ls -1R|grep : | sed 's/:$//g' | sed 's/..//'>times");
mywfile=fopen("Uensemble.dat","w");
myfolders=fopen("times","r");

fscanf(myfolders,"%s",filename);
fscanf(myfolders,"%s",filename);

while(fscanf(myfolders,"%s",filename)==1)
{
fprintf(mywfile,"%s\t",filename);
strcat (filename,"/U_interpolatedPlane.raw");
printf("%s\n",filename);
myfile=fopen(filename,"r");
for (i=0;i<11;i++) {fscanf(myfile,"%s\n",c);printf("%s\n",c);}
u0=0;
v0=0;
w0=0;
 for (i=0;!feof(myfile);i++)
{
fscanf(myfile,"%e %e %e %e %e %e\n",&x[i],&y[i],&z[i],&u[i],&v[i],&w[i]);
u0=u0+u[i];
v0=v0+v[i];
w0=w0+w[i];
N=i+1;
}

u_mean=u0/(double)N;
v_mean=v0/(double)N;
w_mean=w0/(double)N;

fprintf(mywfile,"%e %e %e %e\n",y[1],u_mean,v_mean,w_mean);

fclose(myfile);
}

fclose(mywfile);

}
