#include <stdio.h>
#include<math.h>
void add(unsigned long int a[], unsigned long int b[], unsigned long int addr_31[]){
	unsigned long int carry, c[3];    
	c[0] = a[0]+b[0];
	c[1] = a[1]+b[1];
	carry = c[0]>>31;	c[0] = c[0] & 0x7FFFFFFF;	c[1] = c[1]+carry;
	carry = c[1]>>31;	c[1] = c[1] & 0x7FFFFFFF;	c[2] = carry;
    addr_31[0]=c[0]+(c[2]<<1); addr_31[1]=c[1];
}
void point_negation(unsigned long  int a[],unsigned long  int negation_31[]){
	unsigned long long int c;
	c=a[0]+(a[1]<<31);
	c=0x1FFFFFFFFFFFFFFF-c;
 	negation_31[1]=c>>31;
    negation_31[0] = c & 0x7FFFFFFF;
}
void point_substraction(unsigned long int a[], unsigned long int b[], unsigned long int substraction_31[]){
	unsigned long int negation_311[2] ;
	point_negation(b,negation_311);
	add(a,negation_311,substraction_31);
}
void multr(unsigned long int a[], unsigned long int b[], unsigned long int mult_31[]){
	unsigned long int carry,c[4];
	c[0] = a[0]*b[0];
	c[2] = a[1]*b[1];
	c[1] = (a[0]+a[1])*(b[0]+b[1])-c[0]-c[2];
    mult_31[0]=c[0]+(c[2]>>1);
    mult_31[1]=c[1];
	carry = c[0]>>31;	c[0] = c[0] & 0x7FFFFFFF;	c[1] = c[1]+carry;
	carry = c[1]>>31;	c[1] = c[1] & 0x7FFFFFFF;	c[2] = c[2]+carry;
	carry = c[2]>>31;	c[2] = c[2] & 0x7FFFFFFF;	c[3] = carry;
	mult_31[0]=c[0]+(c[2]<<1);mult_31[1]=c[1]+(c[3]<<1);
  	carry=mult_31[0]>>31; mult_31[0] = mult_31[0] & 0x7FFFFFFF; mult_31[1]=mult_31[1]+carry;
	carry=mult_31[1]>>31;mult_31[1] = mult_31[1] & 0x7FFFFFFF;
    mult_31[0]=mult_31[0]+(carry<<1);
}

void dubling(unsigned long int P_x[], unsigned long int P_z[],
                unsigned long int D_x[],unsigned long int D_z[]){
                unsigned long int S1[2],S2[2],reduced1[2],reduced2[2],A1[2]={ 1447112650,898208916};//A1=(A+2)/4 MOD P
				add(P_x,P_z,reduced2); multr(reduced2,reduced2,S1);
				point_substraction(P_x,P_z,reduced2);multr(reduced2,reduced2,S2);
				point_substraction(S1,S2,reduced1);	multr(S1,S2,reduced2);D_x[0]=reduced2[0];D_x[1]=reduced2[1];
				multr(A1,reduced1,S1);add(S2,S1,reduced2);multr(reduced1,reduced2,S2);D_z[0]=S2[0];D_z[1]=S2[1];
}
void FIELD_inversion(unsigned long int RZ[],unsigned long int F[]){
	unsigned long int CC[2],FF[2];
	FF[0]=RZ[0];FF[1]=RZ[1];
	multr(FF,FF,CC);
	for(int i=0;i<59;i++) {multr(CC,CC,CC);multr(FF,CC,FF);}
	F[0]=FF[0];F[1]=FF[1];


}


void point_addition_EC(unsigned long int P_x[], unsigned long int P_z[],
                     unsigned long int Q_x[],unsigned long int Q_z[],unsigned long int R_x[],unsigned long int R_z[]){

                    unsigned long int S1[2], S2[2], reduced1[2], reduced2[2];
					point_substraction(Q_x,Q_z,reduced1);add(P_x,P_z,reduced2);multr(reduced1,reduced2,S1);
					point_substraction(P_x,P_z,reduced1);add(Q_x,Q_z,reduced2);multr(reduced1,reduced2,S2);
					add(S1,S2,reduced1);multr(P_z,reduced1,reduced2);
					R_x[0]=reduced2[0];R_x[1]=reduced2[1];
					point_substraction(S1,S2,reduced1);multr(P_x,reduced1,reduced2);
					R_z[0]=reduced2[0];R_z[1]=reduced2[1];
            }

void Montgomery_ladder(unsigned long int P_x[], unsigned long int P_z[],
                     unsigned long int R_x[],unsigned long int R_z[],unsigned long int n,unsigned long int l){
                        R_x[0]=P_x[0];R_x[1]=P_x[1];R_z[0]=P_z[0];R_z[1]=P_z[1];
						unsigned long int k,S_x[2], S_z[2],C_x[2], C_z[2],D_x[2],D_z[2],F_x[2], F_z[2], G_x[2], G_z[2];

                        dubling(P_x,P_z,S_x,S_z);
						 int i,ni;
                     	for(i=l-2;i>=0;i--){
                        	ni= n & (1<<i);							
                        	point_addition_EC( S_x, S_z, R_x, R_z, C_x, C_z);dubling( S_x, S_z, D_x, D_z); dubling( R_x, R_z, F_x, F_z);
                        	F_x[0]=(1-ni)*F_x[0];F_x[1]=(1-ni)*F_x[1];F_z[0]=(1-ni)*F_z[0]; F_z[1]=(1-ni)*F_z[1];
                        	D_x[0]=(ni)*D_x[0];D_x[1]=(ni)*D_x[1];D_z[0]=(ni)*D_z[0];D_z[1]=(ni)*D_z[1];
                        	G_x[0]=(1-ni)*C_x[0];G_x[1]=(1-ni)*C_x[1];G_z[0]=(1-ni)*C_z[0]; G_z[1]=(1-ni)*C_z[1]; 
                        	C_x[0]=(ni)*C_x[0];C_x[1]=(ni)*C_x[1];C_z[0]=(ni)*C_z[0];C_z[1]=(ni)*C_z[1];
                        	add(G_x,  D_x, S_x);add(G_z,  D_z, S_z);
							add(F_x,  C_x, R_x);add(F_z,  C_z, R_z);							                     
                     }
                     }

int main(){
	unsigned long int P_x[2]={ 486145363,268485549},P_z[2]={1,0}, RX[2],RZ[2],RX_RZ[2],Z[2],n;
	int i,l;
	n=57646075231573330;
	l=(log(n)/log(2))+1;
	Montgomery_ladder(P_x, P_z,RX,RZ, n,l);
	printf("\n\n X coordinate of nP in projective space :  ");
    for(i=0;i<2;i++)
    printf("%lu ",RX[i]);
	printf("\n\n");

	printf(" Z coordinate of nP in projective space :  ");
    for(i=0;i<2;i++)
    printf("%lu ",RZ[i]);
	printf("\n\n");



	FIELD_inversion(RZ,Z);
	multr(RX,Z,RX_RZ);
		printf("The affine coordinate of nP_X is :  ");
    for(i=0;i<2;i++)
    printf("%lu ",Z[i]);
	printf("\n\n");

}
