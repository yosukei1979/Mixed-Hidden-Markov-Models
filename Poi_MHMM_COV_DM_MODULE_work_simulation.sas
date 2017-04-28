/*******************************************************************************
    NAME    : Poi_MHMM_COV_DM.sas
    TITLE   : Parameter estimation for Poisson_HMM with Covariate by direct mazimization
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2016/01/17
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter nosource;
/********** ���W���[����store�� **********/
libname ModDir ".\module";

proc iml;

start logistic(x);
res=1/(1+EXP(-x));
return(res);
finish;
/********** �f�[�^�ɑ΂�poisson distribution�Ɋ�Â��m�����v�Z���郂�W���[�� **********/
Start calc_res(lambda,beta,_R) global(X,c);
NR=NROW(X);
NL=NCOL(lambda);
do I=1 to NR;
    X1=X[I,c+1];*�f�[�^������endpoint�̗�ɕύX����;
/*    X1_1=X1-1;*/
    do J=1 to NL;
/*        Lambda1=lambda[1,J];*/
        Lambda1=exp(lambda[1,J] + x[I,1:c]*beta` + _R);
        Pois1=PDF('POISSON',X1,lambda1);
/*        if X1^=0 thenX1=0�̎��G���[�ɂȂ�̂����*/
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

/********** 0�ɑ΂���-9999999999��Ԃ�log�֐�"log_"�̍쐬 **********/
start log_(MTX);
    tmp1=choose(MTX<=0,.,MTX);
    tmp2=log(tmp1);
    result=choose(tmp2=.,-9999999999,tmp2);
    return(result);
finish log_;
/********** �����s����s���ɕ������ĉ��������ăx�N�g���ɂ��郂�W���[�� **********/
start tran_m(mtx);/*mtx:�����s��*/
N=nrow(mtx);

mvect=j(1,N*N,1);
do y_=1 to N;
    do z_=1 to N;
    mvect[((y_-1)*N)+z_]=mtx[y_,z_];
    end;
end;
return(mvect);
finish tran_m;
/********** �s�x�N�g������ɁA��,�����č\�����郂�W���[�� **********/
/*20170426�C���F���Adelta��logit�ϊ��E������1�ɂȂ�ϐ��͓��o������cons���w�肵�Ȃ��ōςނ悤�ɕύX*/
start tran_m_inv(rv,lambda,gamma,delta,beta,sigma) global(m,c);/*rv:�s�x�N�g��*/
lambda=j(1,m,0);
gamma=j(m,m,0);
delta=j(1,m,0);
beta=j(1,c,0);
sigma=j(1,1,0);

*��ʉ�����̑�ςȂ̂�m=2�̏ꍇ�ɌŒ�A�x�^�����Ή�;
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
    _R=PDF ('NORMAL', _x, 0, exp(sigma*2));/*�����_������*/
    allprobs=calc_res(lambda,beta,_x);
    *print allprobs;
    n=NROW(allprobs);/*�o�͌n��̐�*/

    foo=delta#allprobs[1,];
    sumfoo=foo[,+];
    if sumfoo=0 then do;
        sumfoo=1;*divided by 0 �h�~;
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
    lalpha=t(lalpha[n,]);
/*    lalpha=lalpha//log_(_R);*/
    _c=max(lalpha);
/*    lalpha2=exp(lalpha-_c);*/
    llk=_c+log_(sum(exp(lalpha-_c)));
	log_R=log_(_R);
	_cllkr=max(llk,log_R);
	_cllkr_v=llk//log_R;
    llkr=_cllkr+log_(sum(exp(_cllkr_v-_cllkr)));
    return(llkr);
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
    
    _integ=j(1,61,0);
    _in2=0;
    do _in=-3 to 3 by 0.1;
        _in2=_in2+1;
        _integ[1,_in2]=mllk_S(_in);
    end;
/*    z=_integ[,+];*/
/*  print _integ;*/
/*  RANG={-1 1};*/
/*call quad(z,"mllk_S",rang);*/
/*logsum ver*/
  _cz=max(_integ);
  _z2=_cz+log_(sum(exp(_integ-_cz)));
  z=_z2;

    RES[1,_nsbj]=z;
/*    print RES;*/
end;
/*print RES;*/
zz=RES[,+];
const=10*log(6/61);
lzz=zz+const;
/*print lzz;*/
return(lzz);
finish mllk;

reset storage=ModDir.Mod_DM;  /* set location for storage */
store module=_all_;* log_ lalphabeta EM viterbi sampling generate_sample tran_m bootstrap;

quit;

/********************************* pgm end! *********************************/

