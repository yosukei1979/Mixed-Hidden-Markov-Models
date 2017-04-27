/*******************************************************************************
    NAME    : Poi_MHMM_COV_DM.sas
    TITLE   : Parameter estimation for Poisson_HMM with Covariate by direct mazimization
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2016/01/17
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter nosource;
/********** モジュールのstore先 **********/
libname ModDir "C:\Users\biostat\Documents\MHMM＿Paper\Program\module";

proc iml;

start logistic(x);
res=1/(1+EXP(-x));
return(res);
finish;
/********** データに対しpoisson distributionに基づき確率を計算するモジュール **********/
Start calc_res(lambda,beta,_R) global(X,c);
NR=NROW(X);
NL=NCOL(lambda);
do I=1 to NR;
    X1=X[I,c+1];*データを見てendpointの列に変更する;
/*    X1_1=X1-1;*/
    do J=1 to NL;
/*        Lambda1=lambda[1,J];*/
        Lambda1=exp(lambda[1,J] + x[I,1:c]*beta` + _R);
        Pois1=PDF('POISSON',X1,lambda1);
/*        if X1^=0 thenX1=0の時エラーになるのを回避*/
/*        Pois1=poisson(lambda1,X1)-poisson(lambda1,X1_1);*/
/*        else Pois1=poisson(lambda1,X1);*/
        POIS=POIS||Pois1;
    End;
    X2=X2||POIS;
    free POIS;

    RES=RES//X2;
    free X2;
end;
/*RES=X||RES;*/
/*print RES;*/
return(RES);

finish calc_res;

/********** 0に対して-9999999999を返すlog関数"log_"の作成 **********/
start log_(MTX);
    tmp1=choose(MTX<=0,.,MTX);
    tmp2=log(tmp1);
    result=choose(tmp2=.,-9999999999,tmp2);
    return(result);
finish log_;
/********** 正方行列を行毎に分解して横結合してベクトルにするモジュール **********/
start tran_m(mtx);/*mtx:正方行列*/
N=nrow(mtx);

mvect=j(1,N*N,1);
do y_=1 to N;
    do z_=1 to N;
    mvect[((y_-1)*N)+z_]=mtx[y_,z_];
    end;
end;
return(mvect);
finish tran_m;
/********** 行ベクトルからλ、Γ,βを再構成するモジュール **********/
/*20170426修正：Γ、deltaにlogit変換・足して1になる変数は導出等してconsを指定しないで済むように変更*/
start tran_m_inv(rv,lambda,gamma,delta,beta,sigma) global(m,c);/*rv:行ベクトル*/
lambda=j(1,m,0);
gamma=j(m,m,0);
delta=j(1,m,0);
beta=j(1,c,0);
sigma=j(1,1,0);

*一般化するの大変なのでm=2の場合に固定、ベタ書き対応;
do i=1 to m;
/*    do j=1 to m;*/
/*    */
/*    if j=1 then*/
/*    gamma[i,j]=rv[,(m*i+j**1)];*/
/*    else gamma[i,j]=rv[,(m*i+j**1)];*/
/*    end;*/
    lambda[,i]=rv[,i];
    delta[,i]=rv[,(m*(m+1)+i)];
end;

gamma11=logistic(rv[,3]);
gamma12=1-gamma11;
gamma22=logistic(rv[,4]);
gamma21=1-gamma22;

gamma[1,1]=gamma11;
gamma[1,2]=gamma12;
gamma[2,1]=gamma21;
gamma[2,2]=gamma22;
delta1=logistic(rv[,5]);
delta2=1-delta1;
delta[,1]=delta1;
delta[,2]=delta2;
do _k=6 to 7;
beta[,_k-5]=rv[,_k];
end;

sigma[1,1]=rv[,8];

finish tran_m_inv;

/********** Log-likelihood **********/
start mllk_s(_x) global(pvt,x,m,c);
_flg=0;
    run tran_m_inv(pvt,lambda,gamma,delta,beta,sigma);
    _R=PDF ('NORMAL', _x, 0, exp(sigma*2));/*ランダム効果*/
    allprobs=calc_res(lambda,beta,_x);
    *print allprobs;
    n=NROW(allprobs);/*出力系列の数*/

    foo=delta#allprobs[1,];
    sumfoo=foo[,+];
    if sumfoo=0 then do;
        sumfoo=1;*divided by 0 防止;
/*        print "w a r n i n g : replacing sumfoo by 1 because sumfoo equals 0";*/
        foo=foo+1/m;
    end;
    lscale=log_(sumfoo);
    foo=foo/sumfoo;
    lalpha=log_(foo)+lscale;
    do _I=2 to n;
        foo= foo*gamma#allprobs[_I,];
        sumfoo=foo[,+];
    if sumfoo=0 then do;
        sumfoo=1;
/*        print "w a r n i n g : replacing sumfoo by 1 because sumfoo equals 0";*/
        foo=foo+1/m;
    end;
        lscale=lscale+log_(sumfoo);
        foo=foo/sumfoo;
        lalpha2=log_(foo)+lscale;
        lalpha=lalpha//lalpha2;
    end;
    lalpha=t(lalpha);
    _c=max(lalpha[,n]);
/*    llk=log_(exp(_c+log_(sum(exp(lalpha[,n]-_c))))) +log_(_R);*/
    llk=_c+log_(sum(exp(lalpha[,n]-_c)));
	llk2=exp(llk)*_R;
    return(llk2);
	skip:
finish mllk_s;

start mllk(_pvt) global(pvt,m,c,x,nsbj,obs);
pvt=_pvt;
*RANG={-10 10};
RES=j(1,nsbj,0);
do _nsbj=1 to nsbj;
    use anal_data;
    read all into x var{COV1 COV2 count} where(SUBJECT=_nsbj);
    close anal_data; 
/*    print _nsbj;*/
/*    print x;*/
/*X=x_all[(_nsbj*obs-9):(_nsbj*obs),];*/
    
    _integ=j(1,41,0);
	_in2=0;
    do _in=-2 to 2 by 0.1;
	    _in2=_in2+1;
	    _integ[1,_in2]=mllk_S(_in);
    end;
	z=_integ[,+];
/*	print _integ;*/
/*	RANG={-1 1};*/
/*call quad(z,"mllk_S",rang);*/
/*logsum ver*/
/*	_cz=max(_integ);*/
/*	_z2=_cz+log_(sum(exp(_integ-_cz)));*/
/*	z=_z2+log(4)-log(41);*/

    RES[1,_nsbj]=z;
/*    print RES;*/
end;
zz=RES[,+];
const=10*log(4/41);
lzz=log_(zz)+const;
/*print zz;*/
return(lzz);
finish mllk;

reset storage=ModDir.Mod_DM;  /* set location for storage */
store module=_all_;* log_ lalphabeta EM viterbi sampling generate_sample tran_m bootstrap;

quit;

/********************************* pgm end! *********************************/

