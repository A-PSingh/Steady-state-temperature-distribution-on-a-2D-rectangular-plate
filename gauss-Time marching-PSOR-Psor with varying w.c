#include<stdio.h>
#include<math.h>
float pi = 3.14;
int main(){
	int M=21,N=41,i,j,iteration1=0;
	float x,y,dx=0.05,dy=0.05,temp;
	float psi[M][N],A[M][N];
    float error1 = 1.0;
    //FILE OPENING
	FILE *fp1,*fp2,*fp3,*fp4;
	fp1=fopen("profile1.dat","w");
	fp2=fopen("profile2.txt","w");fp3=fopen("profile3.txt","w"); fp4=fopen("profile4.txt","w");
	
	//define the commom terms
	float b,aa,ab;
	b=dx/dy;
	ab=b*b;
	aa=0.5/(1+ab);//ac=w*0.5 / (1+ab);

	
	//Boundary conditions
	for(i=0;i<M;i++)
	{ for(j=0;j<N;j++)
	    {  
		  if(j==N-1)   //Top boundary condition
	      psi[i][j]=0.0;
	      else if(i==0)  //Left boundary condition
	      psi[i][j]=0.0;
	      else if(i==M-1) //Right boundary condition
	      psi[i][j]=0.0;
	      else if(j==0) //Bottom boundary condition
	      psi[i][j]=100.0;
	      else           //interior conditions
		  psi[i][j]=0.0;
		}
	}
	//STORING ALL GRID POINT VALUES IN ANOTHER ARRAY A[][]
	for(i=0;i<M;i++)
	  { for(j=0;j<N;j++)
	      {   A[i][j]=psi[i][j];
		  }
	  }

  //POINT GAUSS SEIDEL METHOD
	printf("\nBy Point Gauss seidel Method:\n");
	
	do{
	error1=0.0;
	for(i=1;i<M-1;i++)
	  { for(j=1;j<N-1;j++)
	     {temp= psi[i][j]; 
		 psi[i][j]= aa*( ab*(psi[i][j-1] + psi[i][j+1]) + psi[i-1][j] + psi[i+1][j]);
		 error1 = error1 + fabs(-temp+psi[i][j]);
		 }
	  }
	  iteration1++;
	  printf("\niteratitons = %d\t", iteration1 );
	  printf("Error = %lf",error1);	  
}while(error1>0.01);
	for(i=0;i<M;i++)
	 { for(j=0;j<N;j++)
	   { fprintf(fp1,"%d %d %lf\n",i,j,psi[i][j]); 
	   }
	 }
	fclose(fp1);
	
	
	//TIME MARCHING METHOD
	printf("\n\nFor Timr Marching Method:\n");
	for(i=0;i<M;i++)
	  { for(j=0;j<N;j++)
	      {   psi[i][j]=A[i][j];
		  }
	  }
	  float phi[M][N];
	  
	int iteration2=0;
	double rx=0.25,ry=0.25,error2=1.0;
		do{
	error2=0.0;
	for(i=1;i<M-1;i++)
	  { for(j=1;j<N-1;j++)
	     {temp= psi[i][j]; 
		 phi[i][j]= psi[i][j] + rx*( psi[i+1][j] - 2*psi[i][j] + psi[i-1][j]) + ry*( psi[i][j+1] - 2*psi[i][j] +psi[i][j-1]);
	 
		 error2 = error2 + fabs(temp-phi[i][j]);
		 }
	  }
	   for(i=1;i<M-1;i++)
	  { for(j=1;j<N-1;j++)
	      {   psi[i][j]=phi[i][j];
		  }
	  }	
	  iteration2++;
	  printf("\niteratitons = %d\t", iteration2 );
	  printf("Error = %lf",error2);	  
}while(error2>0.01);
	 
	for(i=1;i<M-1;i++)
	 { for(j=1;j<=N-1;j++)
	   { fprintf(fp2,"%d %d %lf\n",i,j,psi[i][j]); 
	   }
	 }
	fclose(fp2);
	
	// PSOR
	
	  printf("\n\nFor Point Successive over relaxation Method\n");
		for(i=0;i<M;i++)
	  { for(j=0;j<N;j++)
	      {   psi[i][j]=A[i][j];
		  }
	  }	
	int iteration3=0;
	double al,w,r,s,t,u,v,wl,ac,error3=1.0;
	r=cos(pi/(M-1)); s= cos(pi/(N-1)); t=(r+ab*s) / (1+ab); 
	al= pow(t,2.0);
	v=sqrt(1-al);
	w= (2 - 2*v) / al;printf("\nw optimum=%lf",w);
	wl=1-w; ac=w*0.5 / (1+ab);
	
	do{
	error3=0.0;
	for(i=1;i<M-1;i++)
	  { for(j=1;j<N-1;j++)
	     {temp= psi[i][j]; 
		 psi[i][j]= wl*psi[i][j] +ac*( ab*(psi[i][j-1] + psi[i][j+1]) + psi[i-1][j] + psi[i+1][j]);
		 error3 = error3 + fabs(temp-psi[i][j]);
		 }
	  }
	  iteration3++;
	  printf("\niteratitons = %d\t", iteration3 );
	  printf("Error = %lf",error3);	  
}while(error3>0.01);

	for(i=1;i<M-1;i++)
	 { for(j=1;j<=N-1;j++)
	   { fprintf(fp3,"%d %d %lf\n",i,j,psi[i][j]); 
	   }
	 }
fclose(fp3);

//PSOR WITH CHANGING w

	for(i=0;i<M;i++)
	  { for(j=0;j<N;j++)
	      {   psi[i][j]=A[i][j];
		  }
	  }
	  
	  	int iteration4;
	double error4=1.0;
	
	printf("\n\nPSOR with varying w values:");
	 for(w=0.8;w<=2.0;w=w+0.1)
	 {  	for(i=1;i<M-1;i++)
	 { for(j=1;j<N-1;j++)
	    { 
	     psi[i][j]=A[i][j];
		}
	 }
	 iteration4=0;
	 printf("\n\nFor w=%lf:\n",w);
	wl=1-w; ac=(w*0.5) / (1+ab);
		do{
	error4=0.0;
	for(i=1;i<M-1;i++)
	  { for(j=1;j<N-1;j++)
	    {temp= psi[i][j]; 
		 psi[i][j]= wl*psi[i][j] +ac*( ab*(psi[i][j-1] + psi[i][j+1]) + psi[i-1][j] + psi[i+1][j]);
		 error4 = error4 + fabs(temp-psi[i][j]);
		 }
	   }
	    iteration4++;
	  printf("\niteratitons = %d\t", iteration4 );
	  printf("Error = %lf",error4);	  
}while(error4>0.01);
fprintf(fp4,"%d %lf\n",iteration4,w);
	for(i=0;i<M;i++)
	  { for(j=0;j<N;j++)
	      {   psi[i][j]=A[i][j];
		  }
	  }
}
	fclose(fp4);
	
	
}
