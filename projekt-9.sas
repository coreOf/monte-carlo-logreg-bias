%LET SEED = 1235;
%LET NREP = 1000;

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
			DISTRIBUTION = &DIST_NAME;
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
