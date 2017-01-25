/*******************************************************************************
    NAME    : Poi_MHMM_COV_DM.sas
    TITLE   : Parameter estimation for Poisson_HMM with Covariate by direct mazimization
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2016/01/17
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter nosource;
/********** モジュールのstore先 **********/
libname ModDir ".\MODULE";

proc iml;
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
start tran_m_inv(rv,lambda,gamma,delta,beta,sigma) global(m,c);/*rv:行ベクトル*/
lambda=j(1,m,0);
gamma=j(m,m,0);
delta=j(1,m,0);
beta=j(1,c,0);
sigma=j(1,1,0);

do i=1 to m;
    do j=1 to m;
    gamma[i,j]=rv[,(m*i+j)];/*gamma*/
    end;
    lambda[,i]=rv[,i];/*lambda*/
    delta[,i]=rv[,(m*(m+1)+i)];
end;
do _k=(m + m*m + m)+1 to (m+m*m+m)+c;
beta[,_k-(m + m*m + m)]=rv[,_k];
end;

sigma[1,1]=rv[,(m+m*m+m)+c+1];

finish tran_m_inv;

/********** Log-likelihood **********/
start mllk_s(_x) global(pvt,x,m,c);
    run tran_m_inv(pvt,lambda,gamma,delta,beta,sigma);

    _R=PDF ('NORMAL', _x, 0, sigma);/*ランダム効果*/
    allprobs=calc_res(lambda,beta,_x);
    *print allprobs;
    n=NROW(allprobs);/*出力系列の数*/

    foo=delta#allprobs[1,];
    sumfoo=foo[,+];
    if sumfoo=0 then do;
        sumfoo=1;*divided by 0 防止;
        print "warn"||"ing : replacing sumfoo by 1 because sumfoo equals 0";
        foo=foo+1/m;
    end;
    lscale=log_(sumfoo);
    foo=foo/sumfoo;
    lalpha=log_(foo)+lscale;
    do _I=2 to n;
        foo= foo*gamma#allprobs[_I,];
        sumfoo=foo[,+];
    if sumfoo=0 then do;
        /*if sumfoo=0 then sumfoo=1;divided by 0 防止*/
        sumfoo=1;*divided by 0 防止;
        print "warn"||"ing : replacing sumfoo by 1 because sumfoo equals 0";
        foo=foo+1/m;
    end;
        lscale=lscale+log_(sumfoo);
        foo=foo/sumfoo;
        lalpha2=log_(foo)+lscale;
        lalpha=lalpha//lalpha2;
    end;
    
    lalpha=t(lalpha);
    _c=max(lalpha[,n]);
    llk=_c+log_(sum(exp(lalpha[,n]-_c))) + log_(_R);
    return(llk);
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
    
    _integ=j(1,21,0);
    do _in=-1 to 1 by 0.1;
    if _in=-1 then _in2=1;
    else _in2=_in2+1;
        _integ[1,_in2]=mllk_S(_in)*0.1;
    end;
    z=_integ[,+];
    *call quad(z,"mllk_S",rang);

    RES[1,_nsbj]=z;
/*    print RES;*/
end;
zz=RES[,+];
return(zz);
finish mllk;


/*******************************************************************************
    parametric bootstrapによる信頼区間の算出
*******************************************************************************/

/********** 有限離散分布生成関数 **********/
start sampling(m,prob,seed);/*prob:各項目を選び出す確率のベクトル(和が1になる)*/
u=j(1,1,0);
call randseed(seed);
call randgen(u,"uniform");
dens=j(1,m,0);
result=1;
do I=1 to m;

    dens1=prob[,(1:I)];
    dens2=dens1[,+];
    if u<=dens2 then do;
        result=I;
        goto skip;
    end;
    end;
skip:
return(result);
finish sampling;
/********** モデルからのsampling関数 **********/
start generate_sample(n,m,rv,samp,seed);
run tran_m_inv(rv,lambda,gamma,delta);
call randseed(seed);/* set random number seed */
mvect=j(1,m,1);/*ベクトル{1,2,...m}を作成*/
do I=1 to m;
mvect[i]=I;
end;
state=j(n,1,0);
state[1]=sampling(m,delta,seed);
do I=2 to n;
s=state[(I-1),];
t=gamma[s,];
state[I]=sampling(m,t,seed);
end;
samp=j(n,2,0);
do J=1 to n;
ss=state[J];
call randgen(x,'POISSON',lambda[,ss]); 
samp[J,1]=ss;
samp[J,2]=x;
end;
finish generate_sample;

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

/********** bootstrap本体 **********/
start bootstrap(bootiter,rv) global(x,pv,opt,cons);
    do k_=1 to bootiter;
        seed=112345-k_;

        run generate_sample(107,3,rv,x,seed);
        CALL NLPNRA(rc, result, "mllk", rv, opt,cons);

        result_dm=result_dm//result;
    end;

    names={'lambda1' 'lambda2' 'lambda3' 'gamma11' 'gamma12' 'gamma13' 'gamma21' 'gamma22' 'gamma23' 'gamma31' 'gamma32' 'gamma33' 'delta1' 'delta2' 'random_sigma'};

    create boot_dm from result_dm[colname=names];
    append from result_dm;
    close result_dm; 
finish bootstrap;

reset storage=ModDir.Mod_DM;  /* set location for storage */
store module=_all_;* log_ lalphabeta EM viterbi sampling generate_sample tran_m bootstrap;

quit;

/********************************* pgm end! *********************************/

