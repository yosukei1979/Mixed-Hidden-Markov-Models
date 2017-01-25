/*******************************************************************************
    NAME    : Poi_MHMM_COV_DM.sas
    TITLE   : Parameter estimation for Poisson_HMM by direct mazimization
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2012/10/12
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter;
/*******************************************************************************
    èâä˙ê›íË
*******************************************************************************/
/*directory*/
libname data ".\data";
libname ModDir ".\MODULE";
%********** 2-state **********;
proc iml;
reset storage=ModDir.Mod_DM;  /* set location for storage */
load module=_all_;
reset storage=ModDir.Mod_SAMPLE;  /* set location for storage */
load module=_all_;
COV_X={
1 0 10,
2 0 12,
3 0 14,
4 0 16,
5 0 18,
6 1 11,
7 1 13,
8 1 15,
9 1 17,
10 1 17};
nsbj=nrow(COV_X);
%*m : Number of states;
m=2;
c=2;
obs=10;
beta0={1 0.01};

l1=log(1);
l2=log(5);
lambda0=l1||l2;
gamma0={0.4 0.6,
       0.6 0.4};
delta0={0.5 0.5};
SIGMA0=0.25;
*SIGMA0=1;

opt={1,2};
cons={
. . 0 0 0 0 0 0 . . 0.001 . .,
. . 1 1 1 1 1 1 . . .     . .,
. . 1 1 . . . . . . .     0 1,
. . . . 1 1 . . . . .     0 1,
. . . . . . 1 1 . . .     0 1
};
run _simulation_(551,599);
quit;
/********************************* pgm end! *********************************/
