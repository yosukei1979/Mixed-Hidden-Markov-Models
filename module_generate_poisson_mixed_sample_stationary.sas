/*******************************************************************************
    NAME    : generate_poisson_mixed_sample.sas
    TITLE   : generate sample data
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2016/07/25
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter nosource;
/********** ���W���[����store�� **********/
libname ModDir ".\module_stationary";
/********** 2-state **********/
proc iml;
/********** �L�����U���z�����֐� **********/
start sampling(m,prob,seed);/*prob:�e���ڂ�I�яo���m���̃x�N�g��(�a��1�ɂȂ�)*/
u=j(1,1,0);
call randseed(seed);
call randgen(u,"uniform");
dens=j(1,m,0);
result=1;
do m1_I=1 to m;

    dens1=prob[,(1:m1_I)];
    dens2=dens1[,+];
    if u<=dens2 then do;
        result=m1_I;
        goto skip;
    end;
    end;
skip:
return(result);
free result;
finish sampling;
/********** ���f�������sampling�֐� **********/
start generate_sample(
                     n/*��������s��*/
                    ,m/*��Ԑ�*/
                    ,input/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                    ,lambda,gamma/*HMM�̃p�����[�^*/
                    ,beta/*���ϗʂ̌W��*/
                    ,u/*�����_������*/
                    ,result/*�o��*/
                    ,seed/*�V�[�h*/);
call randseed(seed);/* set random number seed */
mvect=j(1,m,1);/*�x�N�g��{1,2,...m}���쐬*/
    do m2_I=1 to m;
        mvect[m2_I]=m2_I;
    end;

gamma1=j(2,2,0);
/*print gamma;*/
gamma11=gamma[1,1];
gamma12=gamma[1,2];
gamma1[1,1]=logistic(gamma11);
gamma1[1,2]=1-gamma1[1,1];
gamma1[2,2]=gamma12;
gamma1[2,1]=1-gamma1[2,2];
/*print gamma1;*/

state=j(n,1,0);
_I_DELTA=I(2);
___d=(_I_DELTA-GAMMA1+1)`;
___e={1,1};
___x=solve(___d,___e);
delta=___x`;

state=j(n,1,0);

state[1]=sampling(m,delta,seed);

    do m3_I=2 to n;
        s=state[(m3_I-1),];
        t=gamma1[s,];
        *if m3_I=2 then print t;
        state[m3_I]=sampling(m,t,seed);
    end;
result=j(n,7,0);
h_mu=input[,2:ncol(input)];
*print h_mu;
    do M2_J=1 to n;
        ss=state[M2_J];
        lambda2=lambda[,ss];
        cov=(h_mu#beta)[,+];
        _dummy=exp(lambda2+cov+u);
        call randgen(x,'POISSON',exp(lambda2+cov+u)); 
        result[M2_J,1]=input[,1];
        result[M2_J,2]=input[,2];
        result[M2_J,3]=input[,3];
        result[M2_J,4]=M2_J;
        result[M2_J,5]=ss;
        result[M2_J,6]=x;
        result[M2_J,7]=u;
    end;
/*print lambda,gamma1,delta,beta,u;*/
finish generate_sample;

start _simulation_data_(  obs/*��������s��*/
                    ,m/*��Ԑ�*/
                    ,X/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                    ,lambda,gamma/*HMM�̃p�����[�^*/
                    ,beta/*���ϗʂ̌W��*/
                    ,_u_sigma/*�����_�����ʂ̕��U*/
                    ,result/*�o��*/
                    ,seed/*�V�[�h*/);
do m4_j=1 to nrow(X);

call randseed(seed+m4_j); 
_u=randnormal(1,0,exp(_u_sigma*2));
    input=X[m4_j,];
    run generate_sample(
                         obs/*��������s��*/
                        ,m/*��Ԑ�*/
                        ,input/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                        ,lambda,gamma/*HMM�̃p�����[�^*/
                        ,beta/*���ϗʂ̌W��*/
                        ,_u/*�����_������*/
                        ,_result/*�o��*/
                        ,seed/*�V�[�h*/);
    result=result//_result;
names={'Subject' 'COV1' 'COV2' 'Visit' 'state' 'COUNT' 'random_effect'};

create anal_data from result[colname=names];
append from result;
close anal_data; 
end;
free result;
finish _simulation_data_;

/********** simulation�{�� **********/
start _simulation_(start_iter,bootiter) global(m,obs,nsbj,lambda0,gamma0,beta0,SIGMA0,COV_X,opt,cons,names,c);

    do k_=start_iter to bootiter;
      submit k_;
      %put Iteration= &k_;
      endsubmit;


/*        seed=112345-k_;*/
        seed=k_;

        run _simulation_data_(  obs/*��������s��*/
                    ,m/*��Ԑ�*/
                    ,COV_X/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                    ,lambda0,gamma0/*HMM�̃p�����[�^*/
                    ,beta0/*���ϗʂ̌W��*/
                    ,SIGMA0/*�����_������*/
                    ,x_all/*�o��*/
                    ,seed/*�V�[�h*/);
        *tgamma0=tran_m(gamma0);
*        pvt=lambda0||tgamma0||delta0||beta0||SIGMA0;
        pvt=lambda0||gamma0||beta0||SIGMA0;
*        pvt=lambda0||tgamma0||delta0||beta0||SIGMA0;
       CALL NLPNRA(rc,result_, "mllk", pvt, opt,cons);
/*     CALL NLPQN(rc,result_, "mllk", pvt, opt);*/
	   result_=result_||rc;
        result_dm=result_dm//result_;
    end;
names={'lambda1' 'lambda2' 'gamma11'  'gamma22' 'beta1' 'beta2' 'random_sigma' 'return_code'};
    create boot_dm from result_dm[colname=names];
    append from result_dm;
    close boot_dm; 

finish _simulation_;


reset storage=ModDir.Mod_SAMPLE;  /* set location for storage */
store module=_all_;* log_ lalphabeta EM viterbi sampling generate_sample tran_m bootstrap;

quit;
/********************************* pgm end! *********************************/
