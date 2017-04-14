#include <stdio.h>
#include <string.h>


int main()
{
FILE *myfile,*mywfile1,*mywfile2;
float x[15000],y[15000],z[15000],u[15000],v[15000],w[15000];
int i, N;
char c[300];

mywfile1=fopen("xfixed.dat","w");
mywfile2=fopen("zfixed.dat","w");

myfile=fopen("2520/U_interpolatedPlane.raw","r");
for (i=0;i<11;i++)      //need to be changed
    {
		fscanf(myfile,"%s\n",c);
		printf("%s\n",c);
	}
	
for (i=0;!feof(myfile);i++)
 {
   fscanf(myfile,"%e %e %e %e %e %e\n",&x[i],&y[i],&z[i],&u[i],&v[i],&w[i]);
   if (x[i]>1.99 && x[i]<2)
   {
	fprintf(mywfile1,"%e %e %e %e\n",z[i],u[i],v[i],w[i]);
   }
   if (z[i]>0.99 && z[i]<1)
   {
	fprintf(mywfile2,"%e %e %e %e\n",x[i],u[i],v[i],w[i]);
   }
 }

fclose(myfile);
fclose(mywfile1);
fclose(mywfile2);

}
