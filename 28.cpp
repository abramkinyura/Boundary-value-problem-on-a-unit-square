#include <math.h>
#include <stdio.h>
//static double pi=atan(1.)*4.; 
int main()
{
int N=80;
int m,i,j,n,k;


/* Как создаются двумерные динамические массивы
int n, m;//n и m – количество строк и столбцов матрицы 
float **matr; //указатель для массива указателей
matr = new float * [n]; //выделение динамической памяти 
                          под массив указателей
for (int i=0; i<n; i++)
  matr[i] = new float [m]; //выделение динамической памяти 
                           для массива значений
*/

//double h=0.015625;

double h=0.0125;

//Объявляем динамические массивы

double *alpha=new double[N];

double *beta=new double[N];

double *normsinpiky=new double[N];

double **scalarproduct=new double*[N];
for (int i=0;i<N;i++)
scalarproduct[i]=new double [N];

double *x=new double[N];
//a11=(double*)malloc(N*sizeof(double));
double *y=new double[N];

double **f=new double*[N];
for (int i=0;i<N;i++)
f[i]=new double [N];

double *a=new double[N];
double *b=new double[N];

double *c=new double[N];

double **u=new double*[N];
for (int i=0;i<N;i++)
u[i]=new double [N];

double **usolution=new double*[N];
for (int i=0;i<N;i++)
usolution[i]=new double [N];

double *lambda=new double[N];
double **coeff=new double*[N];
for (int i=0;i<N;i++)
coeff[i]=new double [N];

FILE *out=fopen("data.txt","w");
//Закончили объявлять динамические массивы



for(i=0;i<N;i++)
{
for(j=0;j<N;j++)
//f[i][j]=(i*h-0.5)*(i*h-0.5)+(j*h*j*h); //Эллиптический параболоид
//f[i][j]=(2*(j*h)*(1-j*h)+2*sin(M_PI*i*h)*sin(M_PI*i*h)*(i*h)*(1-i*h));
//f[i][j]=((i*h)*(j*h));
//f[i][j]=100*j*h;
f[i][j]=1/(i*h);
}


for(k=1;k<=(N-2);k++)
lambda[k]=4*sin(M_PI*k*h/2)*sin(M_PI*k*h/2)/(h*h);

for(n=0;n<N;n++)
x[n]=n*h;

for(m=0;m<N;m++)
y[m]=m*h;

for(k=1;k<=(N-1);k++)
{
for(m=1;m<=(N-1);m++)
normsinpiky[k]+=h*sin(M_PI*k*m*h)*sin(M_PI*k*m*h);

for(n=1;n<=(N-1);n++)
for(j=1;j<=(N-1);j++)
scalarproduct[k][n]+=h*f[n][j]*sin(M_PI*k*y[j]);
}

for(n=1;n<=(N-1);n++)
for(k=1;k<=(N-1);k++)
coeff[k][n]=scalarproduct[k][n]/normsinpiky[k]; 

//Всё необходимое есть, теперь можно писать метод прогонки для трёхдиагональной матрицы





for (k=1;k<=(N-2);k++)
{
for (i=0;i<=(N-1);i++)
b[i]=1.0/(h*h);

for (i=0;i<=(N-1);i++)
c[i]=2.0/(h*h)+sin(M_PI*x[i])*sin(M_PI*x[i])*lambda[k];

for (i=0;i<=(N-1);i++)
a[i]=1.0/(h*h);

for (i=0;i<=(N-1);i++)
{alpha[i]=0;  //Обнулил на всякий случай, прежде чем считать
 beta[i]=0;}

alpha[0]=0;

beta[0]=0;

for (i=0;i<=(N-2);i++)
alpha[i+1]=a[i]/(c[i]-b[i]*alpha[i]);

for (i=0;i<=(N-2);i++)
beta[i+1]=(coeff[k][i]+b[i]*beta[i])/(c[i]-b[i]*alpha[i]);

u[k][N-1]=0; //Так как на границе u равно нулю
//На всякий случай по-другому напишу
//u[k][N-1]=(coeff[k][N-1]+alpha[N-1]*beta[N-1])/(c[N-1]-a[N-1]*alpha[N-1]); //Формула из книги

for(i=(N-2);i>=0;i--)
u[k][i]=alpha[i+1]*u[k][i+1]+beta[i+1]; //Альфа и бета всё равно пересчитываются каждый раз, поэтому для
//них хватит одномерного массива, а u[k][i] надо сохранять. Поэтому добавляем дополнительную переменную k


}

for (n=1;n<=(N-2);n++)
for (m=1;m<=(N-2);m++)
for (k=1;k<=(N-2);k++)
usolution[n][m]+=u[k][n]*sin(M_PI*k*y[m]);


//Записываем условия, что на границе u обращается в ноль
for (n=0;n<=(N-1);n++)
usolution[n][0]=0;

for (m=0;m<=(N-1);m++)
usolution[0][m]=0;


for (i=0;i<N;i++)
for (j=0;j<N;j++)
fprintf(out,"{%f,%f,%f},\n",i*h,j*h, usolution[i][j]);


return 0;
}
