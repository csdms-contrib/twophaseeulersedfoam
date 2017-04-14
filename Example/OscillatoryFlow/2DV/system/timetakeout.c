#include <stdio.h>
#include <string.h>


int main()
{
FILE *myfile,*mywfile,*myfolders;
 float x[15000],y[15000],z[15000],u[15000],v[15000],w[15000];
 int i, N;
char c[300],filename[300];

system("ls -1R|grep : | sed 's/:$//g' | sed 's/..//'>times");
mywfile=fopen("Utimeseries.dat","w");
myfolders=fopen("times","r");

fscanf(myfolders,"%s",filename);
fscanf(myfolders,"%s",filename);
	
	
while(fscanf(myfolders,"%s",filename)==1)
{
fprintf(mywfile,"%s\t",filename);
	strcat (filename,"/centerPatch_U.xy");
printf("%s\n",filename);
myfile=fopen(filename,"r");
	
 for (i=0;!feof(myfile);i++)
{
fscanf(myfile,"%e %e %e %e %e %e\n",&x[i],&y[i],&z[i],&u[i],&v[i],&w[i]);
   if (x[i]<=2 && x[i]>=1.99)
   {
	fprintf(mywfile,"%e %e %e\n",u[i],v[i],w[i]);
   }
}

fclose(myfile);
}

fclose(mywfile);

}
