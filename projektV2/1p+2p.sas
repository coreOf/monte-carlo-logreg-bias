%LET SEED = 1235;
%LET NREP = 10;

/* REGRESSION COEFS */
%LET BETA_10 =  3; %LET BETA_11 = -3;
%LET BETA_20 =  4; %LET BETA_21 = -2; %LET BETA_22 = -2;

%LET DISTRIBUTION_NAMES = NORMAL GAMMA CONTAMINATED_NORMAL UNIFORM;

/* N(0, 1) */
%LET MU    = 0;
%LET SIGMA = 1;

/* GAMMA(0.5, 1) */
%LET SHAPE = 0.5;

/* 0.9 N(0, 2) & 0.1 N(20, 5)  */
%LET MU_1    =  0; %LET SIGMA_1 =  2;
%LET MU_2    = 20; %LET SIGMA_2 =  5;
%LET RATIO = 0.9;

/* UNIFORM(-5, 5) */
%LET A = -5;
%LET B =  5;

%MACRO GENERATE(SAMPLE_N, DISTRIBUTION, RETURN_DATA);
	DATA &RETURN_DATA;
		CALL STREAMINIT(&SEED);
		%IF &DISTRIBUTION = 1 %THEN %DO; /* NORMAL */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("NORMAL", &MU, &SIGMA);
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 2 %THEN %DO; /* GAMMA */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("GAMMA", &SHAPE) - &SHAPE;
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 3 %THEN %DO; /* CONT NORMAL */
			DO I = 1 TO &SAMPLE_N;
				IF I <= &RATIO*&SAMPLE_N THEN
					X = RAND ("NORMAL", &MU_1, &SIGMA_1);
				ELSE
					X = RAND ("NORMAL", &MU_2, &SIGMA_2);
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 4 %THEN %DO; /* UNIFORM */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("UNIFORM", &A, &B);
				OUTPUT;
			END;
		%END;
		DROP I;
	RUN;
%MEND;

%MACRO REPLICATE(SAMPLE_N, DIST_NAME, N_REP, REP_IDX);
	/* GENERATE DATA FOR EACH DISTRIBUTION */
	
	%LET DIST_NAME_CHAR = %SCAN(&DISTRIBUTION_NAMES, &DIST_NAME);  /* Fix za text*/

	%GENERATE(&SAMPLE_N, &DIST_NAME, X_TEMP);
	
	/* CALCULATE BERNOUILLI RESPONSE */
	DATA TEMP_DATA;
		SET X_TEMP;
		ETA = &BETA_10 + &BETA_11 * X;
		PROB = EXP (ETA) / (1 + EXP (ETA));
		Y = RAND ("BERNOULLI", PROB);
	RUN;

	/* FIT LOGREG */
	ODS SELECT NONE;
	PROC LOGISTIC DATA = TEMP_DATA;
		MODEL Y = X;
		ODS OUTPUT ParameterEstimates = PE_TEMP;
	RUN;
	ODS SELECT ALL;

	/* CHECK IF PE_TEMP EXISTS AND HAS DATA */
	%LET PE_EXISTS = 0;
	PROC SQL NOPRINT;
		SELECT COUNT(*) INTO :PE_EXISTS
		FROM DICTIONARY.TABLES
		WHERE LIBNAME = 'WORK' AND MEMNAME = 'PE_TEMP';
	QUIT;

	%IF &PE_EXISTS > 0 %THEN %DO;
		/* EXTRACT ESTIMATES */
		DATA PE_FINAL;
			SET PE_TEMP;
			REPLICATIONS = &N_REP;
			SAMPLE_SIZE = &SAMPLE_N;
			DISTRIBUTION = "&DIST_NAME_CHAR";   /* FIX */

			REP_NUMBER = &REP_IDX;
			
			IF VARIABLE = "Intercept" THEN BETA_0_EST = ESTIMATE;
			ELSE IF VARIABLE = "X" THEN BETA_1_EST = ESTIMATE;
			
			IF VARIABLE IN ("Intercept", "X");
			KEEP REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER VARIABLE ESTIMATE BETA_0_EST BETA_1_EST;
		RUN;
		
		PROC TRANSPOSE DATA = PE_FINAL OUT = PE_TRANSPOSED PREFIX = BETA_;
			BY REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER;
			ID VARIABLE;
			VAR ESTIMATE;
		RUN;
		
		DATA PE_FINAL_CLEAN;
			SET PE_TRANSPOSED;
			BETA_0_EST = BETA_Intercept;
			BETA_1_EST = BETA_X;
			KEEP REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER BETA_0_EST BETA_1_EST;
		RUN;
		
		/* APPEND TO RESULTS */
		PROC APPEND BASE = MC_RESULTS DATA = PE_FINAL_CLEAN FORCE;
		RUN;
	%END;
	%ELSE %DO;
		%PUT WARNING: PROC LOGISTIC failed for sample size &SAMPLE_N, distribution &DIST_NAME, replication &REP_IDX;
	%END;

	/* REMOVE TEMP */
	PROC DATASETS LIBRARY = WORK NOLIST;
		DELETE X_TEMP TEMP_DATA PE_TEMP PE_FINAL PE_TRANSPOSED PE_FINAL_CLEAN;
	RUN;
	QUIT;
%MEND;

%MACRO SIMULATE();
	/* INITIALIZE DATASET */
	DATA MC_RESULTS;
		LENGTH DISTRIBUTION $20;
		REPLICATIONS = .; SAMPLE_SIZE = .; REP_NUMBER = .;
		BETA_0_EST = .; BETA_1_EST = .;
		STOP;
	RUN;

	%DO N_IDX = 1 %TO 3;
		%LET N = %SCAN(50 100 200, &N_IDX);
		%DO DISTRIBUTION_IDX = 1 %TO 4;
			%LET DISTRIBUTION_NAME = %SCAN(&DISTRIBUTION_NAMES, &DISTRIBUTION_IDX);
			%DO REP = 1 %TO &NREP;
				%REPLICATE(&N, &DISTRIBUTION_IDX, &NREP, &REP);
				
				/* TRACK PROGRESS */
				%IF %EVAL(&REP/100) = %EVAL(&REP/100) %THEN %DO;
					%PUT Completed &REP replications for N=&N, Distribution=&DISTRIBUTION_NAME;
				%END;
			%END;
		%END;
	%END;
%MEND;

%SIMULATE();

/* REPEAT TRUE PARAMETERS */
%LET TRUE_BETA_0 = 3;   /* &BETA_10 */
%LET TRUE_BETA_1 = -3;  /* &BETA_11 */

/* CALCULATE BIAS */
PROC SQL;
    CREATE TABLE BIAS_SUMMARY AS
    SELECT 
        SAMPLE_SIZE,
        DISTRIBUTION,
        COUNT(*) AS N_REPLICATIONS,
        
        /* BETA_10 */
        MEAN(BETA_0_EST) AS MEAN_BETA_0_EST,
        STD(BETA_0_EST) AS SD_BETA_0_EST,
        
        /* BETA_11 */
        MEAN(BETA_1_EST) AS MEAN_BETA_1_EST,
        STD(BETA_1_EST) AS SD_BETA_1_EST,
        
        /* ABSOLUTE BIAS */
        MEAN(BETA_0_EST) - &TRUE_BETA_0 AS BIAS_BETA_0,
        MEAN(BETA_1_EST) - &TRUE_BETA_1 AS BIAS_BETA_1,
        
        /* RELATIVE BIAS */
        ((MEAN(BETA_0_EST) - &TRUE_BETA_0) / &TRUE_BETA_0) * 100 AS REL_BIAS_BETA_0_PCT,
        ((MEAN(BETA_1_EST) - &TRUE_BETA_1) / &TRUE_BETA_1) * 100 AS REL_BIAS_BETA_1_PCT,
        
        /* MEAN SQUARE ERROR */
        MEAN((BETA_0_EST - &TRUE_BETA_0)**2) AS MSE_BETA_0,
        MEAN((BETA_1_EST - &TRUE_BETA_1)**2) AS MSE_BETA_1
        
    FROM MC_RESULTS
    GROUP BY SAMPLE_SIZE, DISTRIBUTION
    ORDER BY SAMPLE_SIZE, DISTRIBUTION;
QUIT;

/* REPORT */
PROC PRINT DATA = BIAS_SUMMARY ROUND;
    TITLE "Monte Carlo Bias Analysis Summary";
    TITLE2 "True beta_10 = &TRUE_BETA_0, True beta_11 = &TRUE_BETA_1";
    VAR SAMPLE_SIZE DISTRIBUTION N_REPLICATIONS 
        MEAN_BETA_0_EST BIAS_BETA_0 REL_BIAS_BETA_0_PCT
        MEAN_BETA_1_EST BIAS_BETA_1 REL_BIAS_BETA_1_PCT;
    FORMAT MEAN_BETA_0_EST BIAS_BETA_0 MEAN_BETA_1_EST BIAS_BETA_1 8.4
           REL_BIAS_BETA_0_PCT REL_BIAS_BETA_1_PCT 8.2;
RUN;

proc sort data=BIAS_SUMMARY;
    by sample_size distribution;
run;

/*************************************************************************
* 2. Pomoćni makro za panel–graf (jedan graf = jedan β)                 *
*************************************************************************/
%macro panel_bias(beta=, label=);
    proc sgpanel data=BIAS_SUMMARY noautolegend;
        title "Relativna pristranost &label. po distribuciji i veličini uzorka";
        panelby sample_size / columns=3 rows=1 novarname;  /* n=50|100|200 u 3 stupca */
        vbar distribution / response=&beta. dataskin=crisp
                            datalabel datalabelattrs=(size=7); /* prikaži brojke na stupcima */
        colaxis label="Distribucija prediktora" valueattrs=(size=8) fitpolicy=rotatethin;
        rowaxis label="Relativna pristranost (%)";
    run;
%mend;

/* --- Crtaj za sva tri koeficijenta --- */
%panel_bias(beta=REL_BIAS_BETA_0_PCT , label=β₀);
%panel_bias(beta=REL_BIAS_BETA_1_PCT , label=β₁);
%panel_bias(beta=REL_BIAS_BETA_2_PCT , label=β₂);

/*************************************************************************
* 3. Prosječna apsolutna relativna pristranost po distribuciji (β₀–β₂)   *
*************************************************************************/
proc sql;
    create table ABS_BIAS_AVG as
    select  distribution,
            mean(abs(REL_BIAS_BETA_0_PCT)) as abs_b0 format=8.2,
            mean(abs(REL_BIAS_BETA_1_PCT)) as abs_b1 format=8.2,
            mean(abs(REL_BIAS_BETA_2_PCT)) as abs_b2 format=8.2,
            /* prosjek kroz sva tri koeficijenta */
            mean( mean(abs(REL_BIAS_BETA_0_PCT)),
                  mean(abs(REL_BIAS_BETA_1_PCT)),
                  mean(abs(REL_BIAS_BETA_2_PCT)) ) as avg_all format=8.2
    from    BIAS_SUMMARY
    group   by distribution;
quit;

/* Bar-graf – koja distribucija u prosjeku nosi najveću pristranost? */
proc sgplot data=ABS_BIAS_AVG;
    title "Prosj. apsolutna relativna pristranost (β₀, β₁, β₂) po distribuciji";
    vbar distribution / response=avg_all datalabel dataskin=matte;
    yaxis label="Prosj. |rel. pristranost| (%)";
    xaxis label="Distribucija";
run;

/* 2. GRAFIKONI RELATIVNE PRISTRANOSTI */
ods graphics on / width=6in height=4in;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β1 (X)";
  series x=sample_size y=rel_bias_beta_1_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β2 (X2)";
  series x=sample_size y=rel_bias_beta_2_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β0 (Intercept)";
  series x=sample_size y=rel_bias_beta_0_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

ods graphics off;

%LET SEED = 1235;

/* REGRESSION COEFS */
%LET BETA_20 =  4; %LET BETA_21 = -2; %LET BETA_22 = -2;

%LET DISTRIBUTION_NAMES = NORMAL GAMMA CONTAMINATED_NORMAL UNIFORM;

/* N(0, 1) */
%LET MU    = 0;
%LET SIGMA = 1;

/* GAMMA(0.5, 1) */
%LET SHAPE = 0.5;

/* 0.9 N(0, 2) & 0.1 N(20, 5)  */
%LET MU_1    =  0; %LET SIGMA_1 =  2;
%LET MU_2    = 20; %LET SIGMA_2 =  5;
%LET RATIO = 0.9;

/* UNIFORM(-5, 5) */
%LET A = -5;
%LET B =  5;

%MACRO GENERATE(SAMPLE_N, DISTRIBUTION, RETURN_DATA);
	DATA &RETURN_DATA;
		CALL STREAMINIT(&SEED);
		%IF &DISTRIBUTION = 1 %THEN %DO; /* NORMAL */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("NORMAL", &MU, &SIGMA);
				X2 = RAND("NORMAL", &MU, &SIGMA);
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 2 %THEN %DO; /* GAMMA */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("GAMMA", &SHAPE) - &SHAPE;
				X2 = RAND("GAMMA", &SHAPE) - &SHAPE;
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 3 %THEN %DO; /* CONT NORMAL */
			DO I = 1 TO &SAMPLE_N;
				IF I <= &RATIO*&SAMPLE_N THEN DO;
					X = RAND ("NORMAL", &MU_1, &SIGMA_1);
					X2 = RAND("NORMAL", &MU_1, &SIGMA_1);
				END;
				ELSE DO;
					X = RAND ("NORMAL", &MU_2, &SIGMA_2);
					X2 = RAND("NORMAL", &MU_2, &SIGMA_2);
				END;
				OUTPUT;
			END;
		%END;
		%ELSE %IF &DISTRIBUTION = 4 %THEN %DO; /* UNIFORM */
			DO I = 1 TO &SAMPLE_N;
				X = RAND ("UNIFORM", &A, &B);
				X2 = RAND("UNIFORM", &A, &B);
				OUTPUT;
			END;
		%END;
		DROP I;
	RUN;
%MEND;

%MACRO REPLICATE(SAMPLE_N, DIST_NAME, N_REP, REP_IDX);
	%LET DIST_NAME_CHAR = %SCAN(&DISTRIBUTION_NAMES, &DIST_NAME);
	/* GENERATE DATA FOR EACH DISTRIBUTION */
	%GENERATE(&SAMPLE_N, &DIST_NAME, X_TEMP);
	
	/* CALCULATE BERNOUILLI RESPONSE */
	DATA TEMP_DATA;
		SET X_TEMP;
		ETA = &BETA_20 + &BETA_21 * X + &BETA_22 * X2;
		PROB = EXP (ETA) / (1 + EXP (ETA));
		Y = RAND ("BERNOULLI", PROB);
	RUN;

	/* FIT LOGREG */
	ODS SELECT NONE;
	PROC LOGISTIC DATA = TEMP_DATA;
		MODEL Y = X X2;
		ODS OUTPUT ParameterEstimates = PE_TEMP;
	RUN;
	ODS SELECT ALL;

	/* CHECK IF PE_TEMP EXISTS AND HAS DATA */
	%LET PE_EXISTS = 0;
	PROC SQL NOPRINT;
		SELECT COUNT(*) INTO :PE_EXISTS
		FROM DICTIONARY.TABLES
		WHERE LIBNAME = 'WORK' AND MEMNAME = 'PE_TEMP';
	QUIT;

	%IF &PE_EXISTS > 0 %THEN %DO;
		/* EXTRACT ESTIMATES */
		DATA PE_FINAL;
			SET PE_TEMP;
			REPLICATIONS = &N_REP;
			SAMPLE_SIZE = &SAMPLE_N;
			DISTRIBUTION = "&DIST_NAME_CHAR";
			REP_NUMBER = &REP_IDX;
			
			IF VARIABLE = "Intercept" THEN BETA_0_EST = ESTIMATE;
			ELSE IF VARIABLE = "X" THEN BETA_1_EST = ESTIMATE;
			ELSE IF VARIABLE = "X2"   THEN BETA_2_EST = ESTIMATE;
			
			IF VARIABLE IN ("Intercept", "X", "X2");

			KEEP REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER VARIABLE ESTIMATE BETA_0_EST BETA_1_EST BETA_2_EST;
		RUN;
		
		PROC TRANSPOSE DATA = PE_FINAL OUT = PE_TRANSPOSED PREFIX = BETA_;
			BY REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER;
			ID VARIABLE;
			VAR ESTIMATE;
		RUN;
		
		DATA PE_FINAL_CLEAN;
			SET PE_TRANSPOSED;
			BETA_0_EST = BETA_Intercept;
			BETA_1_EST = BETA_X;
			BETA_2_EST = BETA_X2; 
			KEEP REPLICATIONS SAMPLE_SIZE DISTRIBUTION REP_NUMBER BETA_0_EST BETA_1_EST BETA_2_EST;
		RUN;
		
		/* APPEND TO RESULTS */
		PROC APPEND BASE = MC_RESULTS DATA = PE_FINAL_CLEAN FORCE;
		RUN;
	%END;
	%ELSE %DO;
		%PUT WARNING: PROC LOGISTIC failed for sample size &SAMPLE_N, distribution &DIST_NAME, replication &REP_IDX;
	%END;

	/* REMOVE TEMP */
	PROC DATASETS LIBRARY = WORK NOLIST;
		DELETE X_TEMP TEMP_DATA PE_TEMP PE_FINAL PE_TRANSPOSED PE_FINAL_CLEAN;
	RUN;
	QUIT;
%MEND;



%MACRO SIMULATE();
	/* INITIALIZE DATASET */
	DATA MC_RESULTS;
		LENGTH DISTRIBUTION $20;
		REPLICATIONS = .; SAMPLE_SIZE = .; REP_NUMBER = .; 
		BETA_0_EST = .; BETA_1_EST = .; BETA_2_EST = .;
		STOP;
	RUN;

	%DO N_IDX = 1 %TO 3;
		%LET N = %SCAN(50 100 200, &N_IDX);
		%DO DISTRIBUTION_IDX = 1 %TO 4;
			%LET DISTRIBUTION_NAME = %SCAN(&DISTRIBUTION_NAMES, &DISTRIBUTION_IDX);
			%DO REP = 1 %TO &NREP;
				%REPLICATE(&N, &DISTRIBUTION_IDX, &NREP, &REP);
				
				/* TRACK PROGRESS */
				%IF %EVAL(&REP/100) = %EVAL(&REP/100) %THEN %DO;
					%PUT Completed &REP replications for N=&N, Distribution=&DISTRIBUTION_NAME;
				%END;
			%END;
		%END;
	%END;
%MEND;

%SIMULATE();

/* REPEAT TRUE PARAMETERS */
%LET TRUE_BETA_0 = 4;   /* &BETA_20 */
%LET TRUE_BETA_1 = -2;  /* &BETA_21 */
%LET TRUE_BETA_2 = -2;

/* CALCULATE BIAS */
PROC SQL;
    CREATE TABLE BIAS_SUMMARY AS
    SELECT 
        SAMPLE_SIZE,
        DISTRIBUTION,
        COUNT(*) AS N_REPLICATIONS,
        
        /* BETA_20 */
        MEAN(BETA_0_EST) AS MEAN_BETA_0_EST,
        STD(BETA_0_EST) AS SD_BETA_0_EST,
        
        /* BETA_21 */
        MEAN(BETA_1_EST) AS MEAN_BETA_1_EST,
        STD(BETA_1_EST) AS SD_BETA_1_EST,
        
        MEAN(BETA_2_EST) AS MEAN_BETA_2_EST,    
    	STD(BETA_2_EST)  AS SD_BETA_2_EST,
        
        /* ABSOLUTE BIAS */
        MEAN(BETA_0_EST) - &TRUE_BETA_0 AS BIAS_BETA_0,
        MEAN(BETA_1_EST) - &TRUE_BETA_1 AS BIAS_BETA_1,
        
        MEAN(BETA_2_EST)-&TRUE_BETA_2 AS BIAS_BETA_2, 
        
        /* RELATIVE BIAS */
        ((MEAN(BETA_0_EST) - &TRUE_BETA_0) / &TRUE_BETA_0) * 100 AS REL_BIAS_BETA_0_PCT,
        ((MEAN(BETA_1_EST) - &TRUE_BETA_1) / &TRUE_BETA_1) * 100 AS REL_BIAS_BETA_1_PCT,
        
        ((MEAN(BETA_2_EST)-&TRUE_BETA_2)/&TRUE_BETA_2)*100 AS REL_BIAS_BETA_2_PCT,
        
        /* MEAN SQUARE ERROR */
        MEAN((BETA_0_EST - &TRUE_BETA_0)**2) AS MSE_BETA_0,
        MEAN((BETA_1_EST - &TRUE_BETA_1)**2) AS MSE_BETA_1,
        
        MEAN((BETA_2_EST-&TRUE_BETA_2)**2) AS MSE_BETA_2
        
    FROM MC_RESULTS
    GROUP BY SAMPLE_SIZE, DISTRIBUTION
    ORDER BY SAMPLE_SIZE, DISTRIBUTION;
QUIT;

/* REPORT */
PROC PRINT DATA = BIAS_SUMMARY ROUND;
    TITLE "Monte Carlo Bias Analysis Summary";
    TITLE2 "True BETA_20 = &TRUE_BETA_0, True BETA_21 = &TRUE_BETA_1, True BETA_22 = &TRUE_BETA_2";
    VAR SAMPLE_SIZE DISTRIBUTION N_REPLICATIONS 
        MEAN_BETA_0_EST BIAS_BETA_0 REL_BIAS_BETA_0_PCT
        MEAN_BETA_1_EST BIAS_BETA_1 REL_BIAS_BETA_1_PCT
        MEAN_BETA_2_EST BIAS_BETA_2 REL_BIAS_BETA_2_PCT;
    FORMAT MEAN_BETA_0_EST BIAS_BETA_0 MEAN_BETA_1_EST BIAS_BETA_1 8.4
           REL_BIAS_BETA_0_PCT REL_BIAS_BETA_1_PCT 8.2
           MEAN_BETA_2_EST BIAS_BETA_2 REL_BIAS_BETA_2_PCT 8.2;
           
RUN;



proc sort data=BIAS_SUMMARY;
    by sample_size distribution;
run;

/*************************************************************************
* 2. Pomoćni makro za panel–graf (jedan graf = jedan β)                 *
*************************************************************************/
%macro panel_bias(beta=, label=);
    proc sgpanel data=BIAS_SUMMARY noautolegend;
        title "Relativna pristranost &label. po distribuciji i veličini uzorka";
        panelby sample_size / columns=3 rows=1 novarname;  /* n=50|100|200 u 3 stupca */
        vbar distribution / response=&beta. dataskin=crisp
                            datalabel datalabelattrs=(size=7); /* prikaži brojke na stupcima */
        colaxis label="Distribucija prediktora" valueattrs=(size=8) fitpolicy=rotatethin;
        rowaxis label="Relativna pristranost (%)";
    run;
%mend;

/* --- Crtaj za sva tri koeficijenta --- */
%panel_bias(beta=REL_BIAS_BETA_0_PCT , label=β₀);
%panel_bias(beta=REL_BIAS_BETA_1_PCT , label=β₁);
%panel_bias(beta=REL_BIAS_BETA_2_PCT , label=β₂);

/*************************************************************************
* 3. Prosječna apsolutna relativna pristranost po distribuciji (β₀–β₂)   *
*************************************************************************/
proc sql;
    create table ABS_BIAS_AVG as
    select  distribution,
            mean(abs(REL_BIAS_BETA_0_PCT)) as abs_b0 format=8.2,
            mean(abs(REL_BIAS_BETA_1_PCT)) as abs_b1 format=8.2,
            mean(abs(REL_BIAS_BETA_2_PCT)) as abs_b2 format=8.2,
            /* prosjek kroz sva tri koeficijenta */
            mean( mean(abs(REL_BIAS_BETA_0_PCT)),
                  mean(abs(REL_BIAS_BETA_1_PCT)),
                  mean(abs(REL_BIAS_BETA_2_PCT)) ) as avg_all format=8.2
    from    BIAS_SUMMARY
    group   by distribution;
quit;

/* Bar-graf – koja distribucija u prosjeku nosi najveću pristranost? */
proc sgplot data=ABS_BIAS_AVG;
    title "Prosj. apsolutna relativna pristranost (β₀, β₁, β₂) po distribuciji";
    vbar distribution / response=avg_all datalabel dataskin=matte;
    yaxis label="Prosj. |rel. pristranost| (%)";
    xaxis label="Distribucija";
run;

/* 2. GRAFIKONI RELATIVNE PRISTRANOSTI */
ods graphics on / width=6in height=4in;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β1 (X)";
  series x=sample_size y=rel_bias_beta_1_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β2 (X2)";
  series x=sample_size y=rel_bias_beta_2_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

proc sgplot data=BIAS_SUMMARY;
  title "Relativna pristranost za β0 (Intercept)";
  series x=sample_size y=rel_bias_beta_0_pct / group=distribution markers;
  xaxis label="Veličina uzorka";
  yaxis label="Relativna pristranost (%)";
run;

ods graphics off;

