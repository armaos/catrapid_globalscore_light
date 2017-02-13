#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <new>


using namespace std;

int main( int argc, char *argv[])
{
  if ((argv[1]==NULL)||(argv[2]==NULL)||(argv[3]==NULL))
    {
      printf ("Usage:\n");
      cout<<("Seed + Update Time  + Fourier Coefficients")<<endl;
      exit (0);
    }

  int seed, update, fourier;

  sscanf( argv[1], "%i", &seed);
  sscanf( argv[2], "%i", &update);
  sscanf( argv[3], "%i", &fourier);


  //initial
  //const char *baitprey[] = {"P3", "P3", "R", "R"};
  //nega3 is introduced instead of nega2 to speed things up
        char  baitprey[2][20] = {"posi", "nega4"};
        char  set[3][30]      = {"params/coefficients.pro.txt", "params/coefficients.rna.txt", "params/coefficients.txt"};
  double   ** B ,   ** P,   **B2, ** P2, ** B3 ,   ** P3,   **B4, ** P4 ;
  int      ** L;
  double   ** C ,   *C_rna, *C_pro, C_pro_old, C_rna_old;
  int      i,j,k,l,m,ii;

  // parameters
  int  W=fourier, WMAX,WR;

  // functions
  double        X(    double  **,    double **,  double **,  int,int);
  double       DX(    double  **,    double **,  double,     int, int, double,int, int, double);
  double     corr(    double ***,    double **,  double **,  int, int,    int);
  int        rows(    char *);
  int     columns(    char *);
  double    * FR,    * FRN,   * FP,* FPN, * FN, * FNN;
  double    KT=1;

  cout<<"# Random seed        "<<seed<<endl;
  cout<<"# Updated every      "<<update<<" step(s)"<<endl;
  cout<<"# Parameters         "<<set[2]<<endl;


  int       n=1,r1,r3,r4,r8;
  double    S,D,c,r2,r5,R, oldC, newC,fn,fp,ftot,T[10],SS[200],FFF;
  int       H,H2,H3,H4,H5;


  WR=rows(set[0]);
  cout<<"# Variables          "<<WR<<endl;
  cout<<"# uploading...       "<<endl;


  // BAIT POSI     (protein)
  char *upload = new char [200];
  //char upload[500];

  *upload = '\0'; strcat(upload, "./data/bait."); strcat(upload, baitprey[0]);
  //strcat(upload, ".txt");
   //sizes
   WMAX=columns(upload); if(W>WMAX){W=WMAX;}
      H=rows(upload);
  //
  cout<<"# Positives          "<<H<<endl;
  cout<<"# Coefficients       "<<W<<" / "<<WMAX<<endl;
  //
  ifstream bait(upload);                                                                              cout<<"# "<<upload<<endl;
                           B    = new double * [H] ;         B3    =  new double * [int(H)]  ;
  for (i=0; i<H;      i++){B[i] = new double[WMAX] ;         B3[i] =  new         double[WMAX]  ; }
  for (i=0; i<H;      i++){k=0; for(l=0;l<WMAX;l++){   bait>>B[i][l]; if(B[i][l]!=0){k++;}}       }
   bait.close();

  //PREY  POSI     (rna)
  *upload = '\0'; strcat(upload, "./data/prey."); strcat(upload, baitprey[0]);
  //strcat(upload, ".txt");
  ifstream prey(upload);                                                                              cout<<"# "<<upload<<endl;
                           P    = new double  * [H] ;        P3    = new double  * [int(H)] ;
  for (i=0; i<H;      i++){P[i] = new  double[WMAX] ;        P3[i] = new         double[WMAX]  ; }
  for (i=0; i<H;      i++){k=0; for(l=0;l<WMAX;l++){   prey>>P[i][l]; if(P[i][l]!=0){k++;}}      }
  prey.close();





  // BAIT NEGA
  *upload = '\0'; strcat(upload, "./params/bait."); strcat(upload, baitprey[1]);
  //strcat(upload, ".txt");
  //size
  H2=rows(upload);
  cout<<"# Negatives          "<<H2<<" vs "<<H<<endl;
  if(columns(upload)!=WMAX){cout<<"The Files is Corrupted!"<<endl; exit(0);}
  //
  ifstream bait2(upload);                                                                              cout<<"# "<<upload<<endl;
                            B2    = new double * [H2] ;      B4    = new double * [int(H2)] ;
  for (i=0; i<H2;      i++){B2[i] = new double[WMAX]  ;      B4[i] = new         double[WMAX]  ; }
  for (i=0; i<H2;      i++){k=0; for(l=0;l<WMAX;l++){ bait2>>B2[i][l]; if(B2[i][l]!=0){k++;}}    }
  bait2.close();


  //PREY  NEGA
   *upload = '\0'; strcat(upload, "./params/prey."); strcat(upload, baitprey[1]);
   //strcat(upload, ".txt");
  ifstream prey2(upload);                                                                              cout<<"# "<<upload<<endl;
                            P2  = new double  * [H2] ;       P4    = new double * [int(H2)] ;
  for (i=0; i<H2;      i++){P2[i] = new double[WMAX] ;       P4[i] = new          double[WMAX] ; }
  for (i=0; i<H2;      i++){k=0; for(l=0;l<WMAX;l++){ prey2>>P2[i][l]; if(P2[i][l]!=0){k++;}}    }
  prey2.close();

 FP  = new double  [H] ;  FPN  = new double  [H] ;  FN  = new double  [H2] ;  FNN  = new double  [H2] ;

  //lengths
  if(H> H2){H5=H; }
  if(H<=H2){H5=H2;}


  cout<<"# upload done."<<endl;
  cout<<"# checking coefficients..."<<endl;

  //INTERACTION MATRIX
  C = new double * [W] ;
  //size
  if(columns(set[2])!=W){cout<<"The Files is Corrupted!"<<endl; exit(0);}

  ifstream coeff(set[2]);
  for(i=0; i<W; i++){ C[i] = new double [W] ;}
  for(l=0; l<W; l++){ for(m=0;m<W;m++){ coeff>>C[l][m];}}
  coeff.close();
  cout<<"# OK! running montecarlo...  "<<endl;

   //PROTEIN COEFFICIENTS
  ifstream coeffp(set[0]);
  C_pro = new double [WR] ;
  for(l=0; l<WR; l++){ coeffp>>C_pro[l];}
  coeffp.close();

  //RNA COEFFICIENTS
  C_rna = new double [WR] ;
  ifstream coeffr(set[1]);
  for(l=0; l<WR; l++){ coeffr>>C_rna[l];}
  coeffr.close();

  double TOT,score, scoreT;
  D=-100;
  srand(seed);
  int    w1,w2,k1,k2;
  double w3,w4,C_rmax, C_pmax;

  //MAIN LOOP
  for (j=0; j<n; j++){


    // MOVES
    // rna and protein coefficients
    if(j%(1)==0){
    //cout<<"# move"<<" "<<j<<endl;

    w1=rand()%WR; w2=rand()%WR;
    w3=.05*(( (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5);
    w4=.05*(( (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5);
    // no changes when initializing
    if (j==0){w3=0; w4=0;}
    //
    C_rna_old=C_rna[w2];
    C_pro_old=C_pro[w1];
    //
    C_rna[w2]=C_rna[w2]+w4;
    C_pro[w1]=C_pro[w1]+w3;

    // one coefficient = no move
    if(WR==1){C_rna[w1]=1; C_pro[w2]=1;}

    // normalizes C_rna and C_pro
    C_pmax=-1000; for(m=0;m<WR;m++){ if(C_pmax<=sqrt(C_pro[m]*C_pro[m])){C_pmax=sqrt(C_pro[m]*C_pro[m]);}}
    C_rmax=-1000; for(m=0;m<WR;m++){ if(C_rmax<=sqrt(C_rna[m]*C_rna[m])){C_rmax=sqrt(C_rna[m]*C_rna[m]);}}
                  for(m=0;m<WR;m++){  C_rna[m]=C_rna[m]/C_rmax;          C_pro[m]=C_pro[m]/C_pmax;       }


     // updates rna / protein coefficients
     //       for(i=0;i<H3;    i++)   {for(l=0;l<WMAX;l++){ B3[H3][l]=0; }}
     //       for(i=0;i<H3;    i++)   {for(l=0;l<WMAX;l++){ P3[H3][l]=0; }}
     //       for(i=0;i<H4;    i++)   {for(l=0;l<WMAX;l++){ B4[H4][l]=0; }}
     //       for(i=0;i<H4;    i++)   {for(l=0;l<WMAX;l++){ P4[H4][l]=0; }}
     // H3=0; for(i=0;i<H ;    i=i+WR){for(l=0;l<WMAX;l++){ for(m=0;m<WR;m++){B3[H3][l]=B3[H3][l]+C_pro[m]* B[i+m][l];}} H3++; }
     // H3=0; for(i=0;i<H ;    i=i+WR){for(l=0;l<WMAX;l++){ for(m=0;m<WR;m++){P3[H3][l]=P3[H3][l]+C_rna[m]* P[i+m][l];}} H3++; }
     // H4=0; for(i=0;i<H2;    i=i+WR){for(l=0;l<WMAX;l++){ for(m=0;m<WR;m++){B4[H4][l]=B4[H4][l]+C_pro[m]*B2[i+m][l];}} H4++; }
     // H4=0; for(i=0;i<H2;    i=i+WR){for(l=0;l<WMAX;l++){ for(m=0;m<WR;m++){P4[H4][l]=P4[H4][l]+C_rna[m]*P2[i+m][l];}} H4++; }
     H3=0; for(i=0;i<H ;    i=i+WR){for(l=0;l<WMAX;l++){ B3[H3][l]=0; for(m=0;m<WR;m++){B3[H3][l]=B3[H3][l]+C_pro[m]* B[i+m][l];}} H3++; }
     H3=0; for(i=0;i<H ;    i=i+WR){for(l=0;l<WMAX;l++){ P3[H3][l]=0; for(m=0;m<WR;m++){P3[H3][l]=P3[H3][l]+C_rna[m]* P[i+m][l];}} H3++; }
     H4=0; for(i=0;i<H2;    i=i+WR){for(l=0;l<WMAX;l++){ B4[H4][l]=0; for(m=0;m<WR;m++){B4[H4][l]=B4[H4][l]+C_pro[m]*B2[i+m][l];}} H4++; }
     H4=0; for(i=0;i<H2;    i=i+WR){for(l=0;l<WMAX;l++){ P4[H4][l]=0; for(m=0;m<WR;m++){P4[H4][l]=P4[H4][l]+C_rna[m]*P2[i+m][l];}} H4++; }
               }

    if(j==0){cout<<"# Datasets           "<<H3<<" vs "<<H4<<endl; }
    // protein / rna interaction
    if(j%1==0){
    r1=rand()%W;  r3=rand()%W;
    r2=2*(( (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5);
    // no changes when initializing
    if(j==0){r2=0;}
    for(i=0;i<H3; i++){FP[i]=X(B3, P3,  C, W, i);     }
    for(i=0;i<H4; i++){FN[i]=X(B4, P4,  C, W, i);     }

    oldC=C[r1][r3];
    //
    C[r1][r3]=C[r1][r3]+r2;
    //symmetry
    //C[r3][r1]=C[r1][r3];
    //
    newC=C[r1][r3];

    score=0; fn=0; fp=0; ftot=0;
    for(i=0;i<H3; i++){FPN[i]=DX(B3,  P3,    newC, r1, r3,  FP[i],             i,i,oldC); }
    for(i=0;i<H4; i++){FNN[i]=DX(B4,  P4,    newC, r1, r3,  FN[i],             i,i,oldC);
    //for(i=0;i<H3; i++){FPN[i]=X(B3, P3,  C, W,L, i);     }
    //for(i=0;i<H4; i++){FNN[i]=X(B4, P4,  C, W,L, i);
     for(l=0;l<H3 ; l++){ ftot++; if(FPN[l]-FNN[i]>0){fp++;}  }                              }
    }

    T[0]=2*((fp)/ftot-0.5);
    // T[0]=fp/ftot;
    scoreT=T[0];
    TOT=(T[0]+(scoreT))/2;


    //optimization
    //force gradient
    //if(j==0){if(j%W==0){  KT=0;}}
    //KT=0;

      if(TOT<=D){
      //montecarlo
      r5=(double)rand() / ((double)(RAND_MAX)+(double)(1));
      if(j>=25){
      if(TOT-D> log(r5)*KT){ D=TOT;                                                                               }
      if(TOT-D<=log(r5)*KT){ C[r1][r3]=oldC; /*C[r3][r1]=oldC;*/
                             C_rna[w2]=C_rna_old; C_pro[w1]=C_pro_old;                                              }

      //temperature update
      if( j%W==0){KT=sqrt((TOT-D)*(TOT-D))/0.69+0.000001/0.69;


                  }
               }
      if(j< 25){              C[r1][r3]=oldC; /*C[r3][r1]=oldC;*/
                              C_rna[w2]=C_rna_old; C_pro[w1]=C_pro_old;                                             }
      //output
      cout<<setiosflags(ios::fixed) << setprecision(3)<<setw(9)<<j<<" "<<D<<" "<<T[0]<<" "<<setw(3)<<r1<<" "<<setw(3)<<r3;
      cout<<" "<<setw(6 )<<r2<<" "<<w1<<" "<<setw(6)<<w3<<" "<<w2<<" "<<setw(6)<<w4<<endl;
       }


    //gradient
    if(TOT> D){D=TOT;
      //output

      cout<<setiosflags(ios::fixed) << setprecision(3)<<setw(9)<<j<<" "<<D<<" "<<T[0]<<" "<<setw(3)<<r1<<" "<<setw(3)<<r3;
      cout<<" "<<setw(6)<<r2<<" "<<w1<<" "<<setw(6)<<w3<<" "<<w2<<" "<<setw(6)<<w4<<endl;



      //SAVES
      if (j%update==0)
	{
          //protein/rna
	  ofstream matrix(set[2]);
	  for(l=0;l<W;l++){ for(m=0;m<W;m++){ matrix<<" "<<C[l][m];};  matrix<<""<<endl;}
	  matrix.close();
          //rna
          ofstream matrixr(set[1]);
	  for(l=0;l<WR;l++){ matrixr<<C_rna[l]<<endl;}
	  matrixr.close();
          //protein
          ofstream matrixp(set[0]);
	  for(l=0;l<WR;l++){matrixp<<C_pro[l]<<endl;}
	  matrixp.close();
	  //posi
          ofstream matrip("positives.txt");
	  for(l=0;l<H3;l++){matrip<<FP[l]<<endl;}
	  matrip.close();
	  //posi
          ofstream matrin("negatives.txt");
	  for(l=0;l<H4;l++){matrin<<FN[l]<<endl;}
	  matrin.close();


	}
    }
  }


  return 0;
}



























//FUNCTIONS------------------------------------------------------



//calculates the increment
double X(double **b, double **p, double **c,  int w,int j)
{
  int l, m;
  double r;
  r=0;
  //m>l
  for(l=0; l<w; l++){
    //m=0 (non-symmetric)or m=l  (symmetric)
    for(m=0; m<w; m++)
      {
	// signal
  
	if(c[l][m]!=0){	r=r+double(b[j][l])*(c[l][m])*double(p[j][m]);}
      }
  }
  return r;

}


double DX(double **b, double **p, double newC, int rx, int ry, double f,int j,int k, double oldC)
{
  int r,rrr ;

  // symmetric: rrr>r
  //if (rx>ry){r=ry; rrr=rx;} if(rx<=ry){r=rx; rrr=ry;}
  //f=f+double(b[j][r])*(newC*p[k][rrr]-oldC*p[j][rrr]);

  // symmetric (alternative)
  // if(rx!=ry){ f=f+double(b[j][ry])*(newC-oldC)*double(p[j][rx])+double(b[j][rx])*(newC-oldC)*double(p[j][ry]); }
  // symmetric (alterantive) diagonal
  if(rx==ry){ f=f+double(b[j][ry])*(newC-oldC)*double(p[j][rx]); }
  // non-symmeteric
  if (rx!=ry){ f=f+double(b[j][rx])*(newC-oldC)*double(p[j][ry]); }


  return f;
}



//calculates the correlation
double corr(double ***aggregation, double **aggregation_rate, double **C, int pos, int N, int W)
{
  int t,i,l,m,count;
  double c,s,r,m1,m2,s1,s2,rho;
  double pred_rate[N];
  // initialises
  t=0; c=0; m1=0; m2=0,s1=0; s2=0;
  for (i=0; i<N; i++)
    {
      s=0; r=0; count=0;
      for (l=0; l<W; l++){ for(m=0; m<W; m++){
	  if((aggregation[pos][i][l]!=0)&&(aggregation[pos][i][m]!=0)&&(C[l][m]!=0))
	    {
	      s=s+double(aggregation[pos][i][l])*C[l][m]*double(aggregation[pos][i][m]);
	      count++ ;
	    }}}
      pred_rate[i]=s/(W*W);
    }
  //averages
  for(i=0;i<N;i++){
    m1=m1+pred_rate[i];
    m2=m2+aggregation_rate[pos][i];
  }
  //deviations
  for(i=0;i<N;i++){
    s1=s1+(pred_rate[i]-m1/N)*(pred_rate[i]-m1/N)/N;
    s2=s2+( aggregation_rate[pos][i]-m2/N)*(aggregation_rate[pos][i] -m2/N)/N;
  }
  //correlations
  rho=0;
  for(i=0;i<N;i++){
    rho=rho+(pred_rate[i]-m1/N)/sqrt(s1)*(aggregation_rate[pos][i]-m2/N)/sqrt(s2)/(N);
  }
  return double(rho);
}


double rate(double ***aggregation, int i, double **C, int pos, int N, int W)
{
  int t,l,m,count;
  double c,s,r,m1,m2,s1,s2,rho;
  double rate;
  // initialises
  t=0; c=0; m1=0; m2=0,s1=0; s2=0;

  s=0; r=0; count=0;
  for (l=0; l<W; l++){ for(m=0; m<W; m++){
      if((aggregation[pos][i][l]!=0)&&(aggregation[pos][i][m]!=0)&&(C[l][m]!=0))
	{
	  s=s+double(aggregation[pos][i][l])*C[l][m]*double(aggregation[pos][i][m]);
	  count++ ;
	}}}
  rate=s/(W*W);

  return double(rate);
}

int rows(char *matrix)
{
 ifstream mat5(matrix);
  int nrows=0, ncols=0,col_count=0;
  const char EOL = '\n'; float    temp;
  while   ( mat5.peek()!=EOL) {         mat5>>temp;         ncols++;                    }
  while   (!mat5.eof())       {     if (mat5.peek()==EOL){  mat5.ignore(3,EOL); nrows++; }
			  else{  while (mat5.peek()!=EOL){  mat5 >>temp;col_count++; if (mat5.eof()) break;}}}
      	    mat5.close();
	   return nrows;
}

int columns(char *matrix)
{
 ifstream mat6(matrix);
  int nrows=0, ncols=0,col_count=0;
  const char EOL = '\n'; float    temp;
  while   ( mat6.peek()!=EOL) {         mat6>>temp;         ncols++;                    }
  while   (!mat6.eof())       {     if (mat6.peek()==EOL){  mat6.ignore(3,EOL); nrows++; }
			  else{  while (mat6.peek()!=EOL){  mat6 >>temp;col_count++; if (mat6.eof()) break;}}}
      	    mat6.close();
	   return ncols;
}
