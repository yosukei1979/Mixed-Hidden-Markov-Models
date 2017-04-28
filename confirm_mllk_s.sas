libname data ".\data";
libname ModDir ".\module";
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
beta0={0.5 0.5};

l1=1;
l2=2;
lambda0=l1||l2;
gamma0={1.5 0.5};*対角成分のみ指定し、残りは導出;
delta0={0.5};*一個だけ指定して導出;
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
seed=1;
run _simulation_data_(obs/*生成する行数*/
                    ,m/*状態数*/
                    ,COV_X/*入力行列 1行N列、1:1=SUBJIID 2:N=共変量の値*/
                    ,lambda0,gamma0,delta0/*HMMのパラメータ*/
                    ,beta0/*共変量の係数*/
                    ,SIGMA0/*ランダム効果*/
                    ,x_all/*出力*/
                    ,seed/*シード*/);

_nsbj=1;
    use anal_data;
    read all into x var{COV1 COV2 count} where(SUBJECT=_nsbj);
    close anal_data;

        pvt=lambda0||gamma0||delta0||beta0||SIGMA0;
        _pvt=lambda0||gamma0||delta0||beta0||SIGMA0;
a=MLLK_S(2);
print a;
b=mllk(_pvt);
print b;
quit;
