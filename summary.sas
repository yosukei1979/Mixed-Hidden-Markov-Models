libname library "C:\Users\102650\Documents\git\MHMM\data_timepoint";
data test;
    set library.boot_dm_tp20:;
run;
proc univariate data=test plots;run; 
