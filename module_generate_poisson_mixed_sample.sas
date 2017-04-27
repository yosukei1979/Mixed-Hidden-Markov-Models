/*******************************************************************************
    NAME    : generate_poisson_mixed_sample.sas
    TITLE   : generate sample data
    PRODUCT : SAS R9.3
    AUTHOR  : Y.Inaba
    DATE    : 2016/07/25
*******************************************************************************/
options ls=132 ps=8000 nonumber nodate noxwait noxsync nocenter nosource;
/********** ���W���[����store�� **********/
libname ModDir "C:\Users\biostat\Documents\MHMM�QPaper\Program\module";
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
                    ,lambda,gamma,delta/*HMM�̃p�����[�^*/
                    ,beta/*���ϗʂ̌W��*/
                    ,u/*�����_������*/
                    ,result/*�o��*/
                    ,seed/*�V�[�h*/);
call randseed(seed);/* set random number seed */
mvect=j(1,m,1);/*�x�N�g��{1,2,...m}���쐬*/
    do m2_I=1 to m;
        mvect[m2_I]=m2_I;
    end;
state=j(n,1,0);

delta1=j(1,2,0);
delta1[1,1]=1/(1+EXP(-delta[1,1]));
delta1[1,2]=1-1/(1+EXP(-delta[1,1]));

state[1]=sampling(m,delta1,seed);

gamma1=j(2,2,0);
gamma1[1,1]=1/1+EXP(-gamma[1,1]);
gamma1[1,2]=1-(1/1+EXP(-gamma[1,1]));
gamma1[2,1]=1-(1/1+EXP(-gamma[1,2]));
gamma1[2,2]=1-(1/1+EXP(-gamma[1,2]));



    do m3_I=2 to n;
        s=state[(m3_I-1),];
        t=gamma1[s,];
        state[m3_I]=sampling(m,t,seed);
    end;
result=j(n,6,0);
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
    end;

finish generate_sample;

start _simulation_data_(  obs/*��������s��*/
                    ,m/*��Ԑ�*/
                    ,X/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                    ,lambda,gamma,delta/*HMM�̃p�����[�^*/
                    ,beta/*���ϗʂ̌W��*/
                    ,_u_sigma/*�����_�����ʂ̕��U*/
                    ,result/*�o��*/
                    ,seed/*�V�[�h*/);
call randseed(seed); 
_u=randnormal(1,0,exp(_u_sigma*2));
do m4_j=1 to nrow(X);
    input=X[m4_j,];
    run generate_sample(
                         obs/*��������s��*/
                        ,m/*��Ԑ�*/
                        ,input/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                        ,lambda,gamma,delta/*HMM�̃p�����[�^*/
                        ,beta/*���ϗʂ̌W��*/
                        ,_u/*�����_������*/
                        ,_result/*�o��*/
                        ,seed/*�V�[�h*/);
    result=result//_result;
names={'Subject' 'COV1' 'COV2' 'Visit' 'state' 'COUNT'};

create anal_data from result[colname=names];
append from result;
close anal_data; 
end;
free result;
finish _simulation_data_;

/********** simulation�{�� **********/
start _simulation_(start_iter,bootiter) global(m,obs,nsbj,lambda0,gamma0,delta0,beta0,SIGMA0,COV_X,opt,cons,names,c);

    do k_=start_iter to bootiter;
      submit k_;
      %put Iteration= &k_;
      endsubmit;


/*        seed=112345-k_;*/
        seed=k_;

        run _simulation_data_(  obs/*��������s��*/
                    ,m/*��Ԑ�*/
                    ,COV_X/*���͍s�� 1�sN��A1:1=SUBJIID 2:N=���ϗʂ̒l*/
                    ,lambda0,gamma0,delta0/*HMM�̃p�����[�^*/
                    ,beta0/*���ϗʂ̌W��*/
                    ,SIGMA0/*�����_������*/
                    ,x_all/*�o��*/
                    ,seed/*�V�[�h*/);
        *tgamma0=tran_m(gamma0);
*        pvt=lambda0||tgamma0||delta0||beta0||SIGMA0;
        pvt=lambda0||gamma0||delta0||beta0||SIGMA0;
*        pvt=lambda0||tgamma0||delta0||beta0||SIGMA0;
       CALL NLPNRA(rc,result_, "mllk", pvt, opt,cons);
/*	   CALL NLPQN(rc,result_, "mllk", pvt, opt);*/
        result_dm=result_dm//result_;
    end;
names={'lambda1' 'lambda2' 'gamma11'  'gamma22' 'delta1' 'beta1' 'beta2' 'random_sigma'};
    create boot_dm from result_dm[colname=names];
    append from result_dm;
    close boot_dm; 

finish _simulation_;


reset storage=ModDir.Mod_SAMPLE;  /* set location for storage */
store module=_all_;* log_ lalphabeta EM viterbi sampling generate_sample tran_m bootstrap;

quit;
/********************************* pgm end! *********************************/
