/*******************************************************************************
    NAME    : Poi_MHMM_COV_DM.sas
    TITLE   : Parameter estimation for Poisson_HMM by direct mazimization
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2012/10/12
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter;
options CPUCOUNT=4 THREADS ;

/*******************************************************************************
    初期設定
*******************************************************************************/
/*directory*/
libname data ".\data";
libname ModDir ".\module";
%********** 2-state **********;

%Macro sim_MHMM(start,end,output);
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
SIGMA0=0.5;

opt={1,5,1,1,2};
*cons={
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
run _simulation_(&start.,&end.);
quit;
data data.&output;
    set boot_dm;
run;
%Mend sim_MHMM;
%sim_MHMM(1,1,boot_dm_1);
/*%sim_MHMM(101,200,boot_dm_2);*/
/*%sim_MHMM(201,300,boot_dm_3);*/
/*%sim_MHMM(301,400,boot_dm_4);*/
/*%sim_MHMM(401,500,boot_dm_5);*/
/*%sim_MHMM(501,600,boot_dm_6);*/
/*%sim_MHMM(601,700,boot_dm_7);*/
/*%sim_MHMM(701,800,boot_dm_8);*/
/*%sim_MHMM(801,900,boot_dm_9);*/
/*%sim_MHMM(901,1000,boot_dm_10);*/
/**/
/*%sim_MHMM(1001,1100,boot_dm_11);*/
/*%sim_MHMM(1101,1200,boot_dm_12);*/
/*%sim_MHMM(1201,1300,boot_dm_13);*/
/*%sim_MHMM(1301,1400,boot_dm_14);*/
/*%sim_MHMM(1401,1500,boot_dm_15);*/
/*%sim_MHMM(1501,1600,boot_dm_16);*/
/*%sim_MHMM(1601,1700,boot_dm_17);*/
/*%sim_MHMM(1701,1800,boot_dm_18);*/
/*%sim_MHMM(1801,1900,boot_dm_19);*/
/*%sim_MHMM(1901,2000,boot_dm_20);*/
/*%sim_MHMM(2001,2100,boot_dm_21);*/
/*%sim_MHMM(2101,2200,boot_dm_22);*/
/*%sim_MHMM(2201,2300,boot_dm_23);*/
/*%sim_MHMM(2301,2400,boot_dm_24);*/
/*%sim_MHMM(2401,2500,boot_dm_25);*/
/*%sim_MHMM(2501,2600,boot_dm_26);*/
/*%sim_MHMM(2601,2700,boot_dm_27);*/
/*%sim_MHMM(2701,2800,boot_dm_28);*/
/*%sim_MHMM(2801,2900,boot_dm_29);*/
/*%sim_MHMM(2901,3000,boot_dm_30);*/
/*%sim_MHMM(3001,3100,boot_dm_31);*/
/*%sim_MHMM(3101,3200,boot_dm_32);*/
/*%sim_MHMM(3201,3300,boot_dm_33);*/
/*%sim_MHMM(3301,3400,boot_dm_34);*/
/*%sim_MHMM(3401,3500,boot_dm_35);*/
/*%sim_MHMM(3501,3600,boot_dm_36);*/
/*%sim_MHMM(3601,3700,boot_dm_37);*/
/*%sim_MHMM(3701,3800,boot_dm_38);*/
/*%sim_MHMM(3801,3900,boot_dm_39);*/
/*%sim_MHMM(3901,4000,boot_dm_40);*/
/*%sim_MHMM(4001,4100,boot_dm_41);*/
/*%sim_MHMM(4101,4200,boot_dm_42);*/
/*%sim_MHMM(4201,4300,boot_dm_43);*/
/*%sim_MHMM(4301,4400,boot_dm_44);*/
/*%sim_MHMM(4401,4500,boot_dm_45);*/
/*%sim_MHMM(4501,4600,boot_dm_46);*/
/*%sim_MHMM(4601,4700,boot_dm_47);*/
/*%sim_MHMM(4701,4800,boot_dm_48);*/
/*%sim_MHMM(4801,4900,boot_dm_49);*/
/*%sim_MHMM(4901,5000,boot_dm_50);*/
/*%sim_MHMM(5001,5100,boot_dm_51);*/
/*%sim_MHMM(5101,5200,boot_dm_52);*/
/*%sim_MHMM(5201,5300,boot_dm_53);*/
/*%sim_MHMM(5301,5400,boot_dm_54);*/
/*%sim_MHMM(5401,5500,boot_dm_55);*/
/*%sim_MHMM(5501,5600,boot_dm_56);*/
/*%sim_MHMM(5601,5700,boot_dm_57);*/
/*%sim_MHMM(5701,5800,boot_dm_58);*/
/*%sim_MHMM(5801,5900,boot_dm_59);*/
/*%sim_MHMM(5901,6000,boot_dm_60);*/
/*%sim_MHMM(6001,6100,boot_dm_61);*/
/*%sim_MHMM(6101,6200,boot_dm_62);*/
/*%sim_MHMM(6201,6300,boot_dm_63);*/
/*%sim_MHMM(6301,6400,boot_dm_64);*/
/*%sim_MHMM(6401,6500,boot_dm_65);*/
/*%sim_MHMM(6501,6600,boot_dm_66);*/
/*%sim_MHMM(6601,6700,boot_dm_67);*/
/*%sim_MHMM(6701,6800,boot_dm_68);*/
/*%sim_MHMM(6801,6900,boot_dm_69);*/
/*%sim_MHMM(6901,7000,boot_dm_70);*/
/*%sim_MHMM(7001,7100,boot_dm_71);*/
/*%sim_MHMM(7101,7200,boot_dm_72);*/
/*%sim_MHMM(7201,7300,boot_dm_73);*/
/*%sim_MHMM(7301,7400,boot_dm_74);*/
/*%sim_MHMM(7401,7500,boot_dm_75);*/
/*%sim_MHMM(7501,7600,boot_dm_76);*/
/*%sim_MHMM(7601,7700,boot_dm_77);*/
/*%sim_MHMM(7701,7800,boot_dm_78);*/
/*%sim_MHMM(7801,7900,boot_dm_79);*/
/*%sim_MHMM(7901,8000,boot_dm_80);*/
/*%sim_MHMM(8001,8100,boot_dm_81);*/
/*%sim_MHMM(8101,8200,boot_dm_82);*/
/*%sim_MHMM(8201,8300,boot_dm_83);*/
/*%sim_MHMM(8301,8400,boot_dm_84);*/
/*%sim_MHMM(8401,8500,boot_dm_85);*/
/*%sim_MHMM(8501,8600,boot_dm_86);*/
/*%sim_MHMM(8601,8700,boot_dm_87);*/
/*%sim_MHMM(8701,8800,boot_dm_88);*/
/*%sim_MHMM(8801,8900,boot_dm_89);*/
/*%sim_MHMM(8901,9000,boot_dm_90);*/
/*%sim_MHMM(9001,9100,boot_dm_91);*/
/*%sim_MHMM(9101,9200,boot_dm_92);*/
/*%sim_MHMM(9201,9300,boot_dm_93);*/
/*%sim_MHMM(9301,9400,boot_dm_94);*/
/*%sim_MHMM(9401,9500,boot_dm_95);*/
/*%sim_MHMM(9501,9600,boot_dm_96);*/
/*%sim_MHMM(9601,9700,boot_dm_97);*/
/*%sim_MHMM(9701,9800,boot_dm_98);*/
/*%sim_MHMM(9801,9900,boot_dm_99);*/
/*%sim_MHMM(9901,10000,boot_dm_100);*/


/*data data.boot_dm;*/
/*  set data.boot_dm_1-data.boot_dm_10;*/
/*run;*/
/********************************* pgm end! *********************************/
