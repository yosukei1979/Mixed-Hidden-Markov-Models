libname data ".\data";
libname ModDir ".\module_stationary";
proc iml;
reset storage=ModDir.Mod_DM;  /* set location for storage */
load module=_all_;
reset storage=ModDir.Mod_SAMPLE;  /* set location for storage */
load module=_all_;
COV_X={
1 0 1,
2 0 1.2,
3 0 1.4,
4 0 1.6,
5 0 1.8,
6 1 1.1,
7 1 1.3,
8 1 1.5,
9 1 1.7,
10 1 1.9};
nsbj=nrow(COV_X);
%*m : Number of states;
m=2;
c=2;
obs=10;
beta0={1 1.5};

l1=-1;
l2=1;
lambda0=l1||l2;
gamma0={1.5 0.5};*�Ίp�����̂ݎw�肵�A�c��͓��o;
*SIGMA0=0.25;
SIGMA0=2;

opt={1,5,1,1,1};
cons={
. . 0 0 0 0 0 0 . . 0.001 . .,
. . 1 1 1 1 1 1 . . .     . .,
. . 1 1 . . . . . . .     0 1,
. . . . 1 1 . . . . .     0 1,
. . . . . . 1 1 . . .     0 1
};
cons={
. . . . . . . -10 ,
. . . . . . . .
};        

        pvt=lambda0||gamma0||beta0||SIGMA0;
_nsbj=1;
    use anal_data;
    read all into x var{COV1 COV2 count} where(SUBJECT=_nsbj);
    close anal_data; 
   
run tran_m_inv(pvt,lambda,gamma,beta,delta,sigma);
a=MLLK_S(2,lambda,gamma,delta,beta);
print a;
b=mllk(pvt);
print b;
quit;
