
/*                           Group - 05 (Cryptographic and Security Implementation)
                            Asif Kalam     - CrS2001
                            Subir Das      - CrS2010
                            Pousali Dey    - Crs2023
                            Swagata Sasmal - CrS2024                      */
#include <stdio.h>

void swap(int,int,int);

int main()
{   int x1 , x2 , n;
    printf("\n Enter x1 : ");
    scanf("%d", &x1);
    printf("\n Enter x2 : ");
    scanf("%d", &x2);
    
    printf("\n x1 = %d ",x1);
    printf("\n x2 = %d ",x2);
    
    printf("\n\n Enter 1 if wish to swap, else 0 : ");
    scanf("%d", &n);
    
    swap(x1,x2,n);  
    
    return 0;
}

void swap(int x1, int x2, int n){
    
    int t;
    t = (x1 ^ x2) & (-n);
    x2 = x2 ^ t;
    x1 = x1 ^ t;
    
    printf("\n x1 = %d ",x1);
    printf("\n x2 = %d ",x2);
}