*! version 1.02 16 December 2015
*! version 1.01 22 September 2015
*! version 1.0 6 May 2015
*! Authors: Monica Hernandez and Steve Pudney

/***********************************************************/
/*     Generalized bivariate ordinal regression model      */
/***********************************************************/


program bicop
	version 13.1
	if replay() {
		if ("`e(cmd)'" != "bicop") {
			noi di in red "results for bicop not found"
			exit 301
		}
		Replay `0'
	}
	else Estimate `0'
end




program Replay, eclass
	syntax [, Level(cilevel)]
	
	local copula `e(copula)'
	local mixture `e(mixture)'
	
	if "`copula'" == "indep" & "`mixture'" == "none" {
		ml display, level(`level') nolstretch
	}
	else {
		ml display, level(`level') nolstretch plus 
	}
	
	//	Set up diparm for dependency parameter to be reported in natural form
	local temp = -invnormal((1-(`level'/100))/2) // use level to calculate z value for ci
	
	if "`copula'" == "gaussian" {
		_diparm depend, tanh prob notab
		 output_line "theta" `r(est)' `r(se)' `r(z)' `r(p)' `temp'		 
	}
    if "`copula'" == "frank" {
       	_diparm depend, f(@) derivative(1) prob notab
		output_line "theta" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
	}
    if "`copula'" == "clayton" {
      	_diparm depend, f(exp(@)) derivative(exp(@)) prob notab 
		output_line "theta" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
	}
    if "`copula'" == "joe"|"`copula'" == "gumbel" {
       	_diparm depend, f(exp(@)+1) derivative(exp(@)) prob notab
		output_line "theta" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
	}


	if "`copula'" != "indep" {	
		matrix theta = r(est)
		mat colnames theta = theta
	}
	
	if "`copula'" == "indep" {
		matrix extpar = 1
	}
	else {
		matrix extpar = theta
	}
	
	// set up diparm for mixture parameters to be reported in natural form
	
	if "`mixture'" == "both" | "`mixture'" == "mix1" | "`mixture'" == "equal" {
		_diparm pu1, invlogit prob notab
		output_line "pi_u_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = r(est)
		_diparm pu1, func(1/(1+exp(@))) der(-exp(@)/((1+exp(@))^2)) prob notab
		output_line "pi_u_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = u , r(est)
		_diparm pu1 mu2, func(-(@2)/exp(@1)) der(@2/exp(@1) (-1/exp(@1))) prob notab
		output_line "mean_u_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = u , r(est)
		_diparm mu2, func(@) der(1) prob notab
		output_line "mean_u_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = u , r(est)
		_diparm pu1 mu2 su2, func((1+exp(@1)-((@3^2)+(@2^2)))/exp(@1)-(-@2/exp(@1))^2) der(1-(1+exp(@1)-((@3^2)+(@2^2)))/exp(@1)+2*(@2/exp(@1))^2 ((@2/exp(@1))*(-2-exp(-@1))) (-2*@3/exp(@1))) prob notab			
		output_line "var_u_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = u , r(est)
		_diparm su2, func(@^2) der(2*@) prob notab
		output_line "var_u_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix u = u , r(est)
		mat colnames u = pi_u_1 pi_u_2 mean_u_1 mean_u_2 var_u_1 var_u_2
		matrix extpar = extpar,u
	}
		
	if "`mixture'" == "both" | "`mixture'" == "mix2" {
		_diparm pv1, invlogit prob notab
		output_line "pi_v_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = r(est)
		_diparm pv1, func(1/(1+exp(@))) der(-exp(@)/((1+exp(@))^2)) prob notab
		output_line "pi_v_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = v , r(est)
		_diparm pv1 mv2, func(-(@2)/exp(@1)) der(@2/exp(@1) (-1/exp(@1)))prob notab
		output_line "mean_v_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = v , r(est)
		_diparm mv2, func(@) der(1) prob notab
		output_line "mean_v_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = v , r(est)
		_diparm pv1 mv2 sv2, func((1+exp(@1)-((@3^2)+(@2^2)))/exp(@1)-(-@2/exp(@1))^2) der(1-(1+exp(@1)-((@3^2)+(@2^2)))/exp(@1)+2*(@2/exp(@1))^2 ((@2/exp(@1))*(-2-exp(-@1))) (-2*@3/exp(@1)))prob notab
		output_line "var_v_1" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = v , r(est)
		_diparm sv2, func(@^2) der(2*@) prob notab
		output_line "var_v_2" `r(est)' `r(se)' `r(z)' `r(p)' `temp'
		matrix v = v , r(est)
		mat colnames v = pi_v_1 pi_v_2 mean_v_1 mean_v_2 var_v_1 var_v_2
		matrix extpar = extpar,v
	}
	
	if  "`copula'" != "indep" {
		ereturn matrix extpar extpar
		di as text "{hline 13}{c BT}{hline 64}"
	} 
	
	else if  "`copula'" == "indep"  & colsof(extpar) > 1  {
		matrix extpar = extpar[1, 2.. colsof(extpar) ]
		ereturn matrix extpar extpar	
		di as text "{hline 13}{c BT}{hline 64}"
	}
	
	
	
	
	if "`diparm'" != "" {
		matrix auxpar=r(table)
		matrix auxpar = auxpar[1, "_diparm1:"]
		matrix colnames auxpar = _:
		ereturn matrix auxpar auxpar
	}
	
	if `"$bc_waldtest"' != ""  {
		display `"$bc_waldtest"'
	}
	
	if `"$bc_waldindp"' != ""  {
		display `"$bc_waldindp"'
	}
end




program Estimate, eclass sortpreserve

	local cmdline `"bicop `0'"'
/* drop global variables to be used */
foreach global in offo1 offo2 Nthr1 Nthr2 copu mix kx1 kx2 waldtest waldindp {
	macro drop bc_`global'
}

/* syntax taken from biprobit and modified*/
/* two syntax, handle one at a time */

	gettoken first : 0, match(paren)	
	if "`paren'" == "" {
	
/* syntax 1, bivariate model with common covariate vector 
										- parenthesis has not been found */
										
		gettoken dep1 0:0, parse (" =,[")
		_fv_check_depvar `dep1' /* check that dep var isn't a factor variable */
		tsunab dep1 : `dep1' /* Expand a variable with time series operators */
		rmTS `dep1' /* deals with time series operators in the dependent variable*/
		confirm variable `r(rmTS)'
		gettoken dep2 0:0, parse (" =,[")
		_fv_check_depvar `dep2'
		tsunab dep2 : `dep2'
		rmTS `dep2'
		confirm variable `r(rmTS)'
		gettoken junk left :0, parse ("=")
		if "`junk'" == "=" {
			local 0 "`left'"
		}
		syntax [varlist(default=none ts fv)] [if] [in] [pw fw iw] /*
			*/[, COPula(string) MIXture(string) /* 
			*/ Robust CLuster(varname)   	/*
			*/ offset1(varname) offset2(varname)    /*
			*/ Level(cilevel) 		/* 
			*/ noLog FROM(string) SEarch(string) Repeat(string)		/*
			*/ /*ITERate(passthru)*/ VCE(passthru) Constraints(string) * ]
		
		local ind1 `varlist'
		local ind2 `varlist'

		
        local dep1n "`dep1'"
        local dep1n : subinstr local dep1 "." "_"
        local dep2n "`dep2'"
        local dep2n : subinstr local dep2 "." "_"

        local e_neq = 2
		
		//mlopts mlopts, `options'
		local option0 `options'
		marksample touse
		markout `touse' `dep1' `dep2' `offset1' `offset2', strok
	
	}
	else {
// syntax 2, different covariate vectors - with parentheses
// get first equation 
		gettoken first 0:0, parse(" ,[") match(paren)
		local left "`0'"
		local junk: subinstr local first ":" ":", count(local number) /* This 
											picks up equation name if given */
		if "`number'" == "1" {
			gettoken dep1n first: first, parse(":")
			gettoken junk first: first, parse(":")
		}
		local first : subinstr local first "=" " "
		gettoken dep1 0: first, parse(" ,[") 
		_fv_check_depvar `dep1'
		tsunab dep1: `dep1'
		rmTS `dep1' 
		confirm variable `r(rmTS)'
		if "`dep1n'" == "" {
			local dep1n "`dep1'"
		}
		syntax [varlist(default=none ts fv)] [, OFFset(varname numeric)]
		local ind1 `varlist'
		local offset1 `offset' 
		
// get second equation 
		local 0 "`left'"
		gettoken second 0:0, parse(" ,[") match(paren)
		if "`paren'" != "(" {
			dis in red "two equations required"
			exit 110
		}
		local left "`0'"
		local junk : subinstr local second ":" ":", count(local number)
		if "`number'" == "1" {
			gettoken dep2n second: second, parse(":")
			gettoken junk second: second, parse(":")
		}
		local second : subinstr local second "=" " "
		gettoken dep2 0: second, parse(" ,[") 
		_fv_check_depvar `dep2'
		tsunab dep2: `dep2' 
		rmTS `dep2'
		confirm variable `r(rmTS)'
		if "`dep2n'" == "" {
			local dep2n "`dep2'"
		}
		syntax [varlist(default=none ts fv)] [, OFFset(varname numeric) ]
		local ind2 `varlist'
		local offset2 `offset' 

// remaining options 
		local 0 "`left'"
		syntax [if] [in] [pw fw iw] [, COPula(string) MIXture(string)   /*
			*/ Robust Cluster(varname) /*ITERate(passthru)*/ /*
			*/ Level(cilevel)  	   /*
			*/ noLOG FROM(string) SEarch(string) Repeat(string) VCE(passthru) * ]

		local option0 `options'
		marksample touse
		markout `touse' `dep1' `dep2' `ind1' `ind2' `offset1'  `offset2'
	}


//	Various options...	
	local wtype `weight'
    local wtexp `"`exp'"'
    if "`weight'" != "" { 
		local wgt `"[`weight'`exp']"'  
	}
	_get_diopts diopts option0, `option0'
	mlopts mlopts, `option0'
	local cns `s(constraints)'	
	if "`cluster'" ! = "" { 
		local clopt "cluster(`cluster')" 
	}
	_vce_parse, argopt(CLuster) opt(Robust oim opg) old: ///
                [`weight'`exp'], `vce' `clopt' `robust'	
        local cluster `r(cluster)'
        local robust `r(robust)'
	if "`cluster'" ! = "" { 
		local clopt "cluster(`cluster')" 
	}
	local cns `s(constraints)'
	
	if `"`robust'"' != "" {
		local crtype crittype("log pseudolikelihood")
	}


// test collinearity 
	local coll `s(collinear)'
	_rmcoll `ind1' `wgt' if `touse', `coll'
	local ind1 "`r(varlist)'"
	_rmcoll `ind2' `wgt' if `touse', `coll'

	
	local ind2 "`r(varlist)'"
//	significance level
	if "`level'" != "" {
		local level "level(`level')"
	}
//	offsets

	if "`offset1'" != "" {
		local offo1 "offset(`offset1')" 
		qui gen double bc_offo1 = `offset1'
	}
    if "`offset2'" != "" { 
		local offo2 "offset(`offset2')" 
		qui gen double bc_offo2 = `offset2'
	}
//	noisy or quiet
	qui {
		if "`log'" == "" {
			local log "noisily"
        }
        else    local log "quietly"
//	construct recoded dependent variables taking values 1,2,...
		tempvar y1 y2
		egen `y1'=group(`dep1')
		egen `y2'=group(`dep2')

//	Check sample size - return error if no observations
		count if `touse' 
		local N = r(N)
		if r(N) == 0 {
			error 2000
		}
// Return error if only one category 
		tab `y1' if `touse'
		global bc_Nthr1 = r(r) -1
		tab `y2' if `touse'
		global bc_Nthr2 = r(r) -1
		if $bc_Nthr1 == 0 {
			di in red "There is no variation in `dep1'"
			exit 2000
		}
		if $bc_Nthr2 == 0 {
			di in red "There is no variation in `dep2'"
			exit 2000
		}
		
//	Check validity of copula declaration & set up global for copula name
		capture assert 	"`copula'" == "" | "`copula'" == "gaussian"  /*
			*/ |"`copula'" == "frank" |"`copula'" == "clayton"  /*
			*/ | "`copula'" == "joe"|"`copula'" == "gumbel" |"`copula'" == "indep" 
		if _rc == 9 {
			di in red "Copula must be indep, gaussian, frank, clayton, gumbel or joe"
			exit 198
		}
		if "`copula'" == "" {
			local copula "gaussian"
		}	
		global bc_copu = "`copula'"

//	Check validity of mixture option & set up global for mixture choice
		capture assert 	"`mixture'" == "" | "`mixture'" == "none" |  "`mixture'" == "equal" | /*
		*/ "`mixture'" == "mix1" | "`mixture'" == "mix2" | "`mixture'" == "both" 
		if _rc == 9 {
			di in red "Mixture must be none, mix1, mix2, equal or both"
			exit 198
		}
		if "`mixture'" == "" {
			local mixture "none"			
		}	
		global bc_mix = "`mixture'"	 
	}


// handling from	
	if "`from'" == "" {
		tempname a from0 mixname
//	Ordered probits to supply starting values
/*
		noi display in green "*************************************************"
		noi display in green "Univariate ordered probits for starting values"
		noi display in green "*************************************************"
*/
		qui oprobit `dep1' `ind1' if `touse' `wgt', `offo1' /* // add quietly
			*/ iter(`=min(1000,c(maxiter))') /*`mlopts'*/ constraints(`constraints')
			
		scalar ll_ind=e(ll)
		if _rc == 0 {
			tempname cb1 ct1 tmp 
			global bc_kx1=e(k)-e(k_aux)
			mat `tmp' = e(b)
			if `"`ind1'"' != ""  {
				mat `cb1'=`tmp'[1,"`dep1':"] // ml model does not allow
					//to pass an equation with ,nocons and no variables 
					// so the model with no variables needs to be handled
					// separately
			}
			
			local nthr1=$bc_Nthr1
			mat `ct1'=`tmp'[1,"cut1:_cons".."cut`nthr1':_cons"]
			local cuts
			forvalues k=1/`nthr1' {
				local cuts "`cuts' /cuteq1_`k'"
				local ceqs "`ceqs' cuteq1_`k'"
			}
			mat coleq `ct1'=`ceqs' 
			mat colnames `ct1'=_cons
		}	
				
		qui oprobit `dep2' `ind2' if `touse' `wgt', `offo2' /*  
			*/ iter(`=min(1000,c(maxiter))') /*`mlopts'*/ constraints(`constraints')
		scalar ll_ind=ll_ind+e(ll)
		if _rc == 0 {
			tempname cb2 ct2 tmp
			global bc_kx2=e(k)-e(k_aux)
			mat `tmp' = e(b)
			if `"`ind2'"' != ""  {
				mat `cb2'=`tmp'[1,"`dep2':"]
			}	
			local nthr2=$bc_Nthr2
			mat `ct2'=`tmp'[1,"cut1:_cons".."cut`nthr2':_cons"]
			local cuts 
			local ceqs
			forvalues k=1/`nthr2' {
				local cuts "`cuts' /cuteq2_`k'"
				local ceqs "`ceqs' cuteq2_`k'"
			}
			mat coleq `ct2'=`ceqs'
			mat colnames `ct2'=_cons
		}
		global ll_ind=ll_ind
		noi di "LogL for independent ordered probit model " ll_ind  
		
//	Insert into vector of start values

		if "`ind1'" != "" & "`ind2'" != "" {
			//mat `from0' = `cb1',`cb2',`ct1',`ct2'
			mat `from0' = `cb1',`cb2', `ct1',`ct2'
		}
		else if ("`ind1'" == "") & ("`ind2'" == "") {
			mat `from0' = `ct1',`ct2'
		}
		else if (("`ind1'" != "") & ("`ind2'" == ""))  { 
			mat `from0' = `cb1', `ct1',`ct2'
		}
		else if (("`ind1'" == "") & ("`ind2'" != "")) { 
			mat `from0' = `cb2', `ct1',`ct2'
			
		}
		
	
		
		//	Starting value for transformed dependency parameter if not independent
		if "`copula'" != "indep" {
			mat `a' = (1.2)
			mat colnames `a' = depend:_cons
			mat `from0' = `from0',`a'
		}
		if "`mixture'" == "equal" | "`mixture'" == "mix1" {
			mat `mixname'=-0.01,-0.01,1.01
			mat colnames `mixname' = pu1:_cons mu2:_cons su2:_cons
			mat `from0'=`from0',`mixname'
		}	
 		if "`mixture'" == "mix2" {
			mat `mixname'=-0.01,-0.01,1.01
			mat colnames `mixname' = pv1:_cons mv2:_cons sv2:_cons
			mat `from0'=`from0',`mixname'
		}	
      	if "`mixture'" == "both" {	
			mat `mixname'=0.8,-0.35,1.4,0.8,-0.35,1.4
			mat colnames `mixname' = pu1:_cons mu2:_cons su2:_cons pv1:_cons /// 
									mv2:_cons sv2:_cons
			mat `from0'=`from0',`mixname'
		}
		local from="`from0'"
		mat `from'=`from0'
		
	}

	else {
		local nn : word count `ind1'
		global bc_kx1= `nn'
		local nn : word count `ind2'
		global bc_kx2= `nn'
		
		/*running oprobits to get ll_ind and starting values if needed*/
		qui oprobit `dep1' `ind1' if `touse' `wgt', `offo1' /* // add quietly
			*/ iter(`=min(1000,c(maxiter))') /*`mlopts'*/ constraints(`constraints')
		scalar ll_ind=e(ll)
		qui oprobit `dep2' `ind2' if `touse' `wgt', `offo2' /*  // add quietly
			*/ iter(`=min(1000,c(maxiter))') /*`mlopts'*/ constraints(`constraints')
		scalar ll_ind=ll_ind+e(ll)
		global ll_ind=ll_ind
		noi di "LogL for independent ordered probit model " ll_ind  
		
	}
//	***  PARAMETER NAMES  ***

//	threshold parameters
	forvalues  i = 1(1)$bc_Nthr1 {
		local thr10 `"`thr10' /cuteq1_`i'"'
	}
	forvalues  i = 1(1)$bc_Nthr2 {
		local thr20 `"`thr20' /cuteq2_`i'"'
	}



//	coefficient names

	if "`ind1'"!="" & "`ind2'"!="" {
		local equs0 `"`equs0' (`dep1n': `y1' = `ind1', `offo1'  nocons  )"'
	    local equs0 `"`equs0' (`dep2n': `y2' = `ind2', `offo2'  nocons  )"'
		local equs0 `"`equs0' `thr10' `thr20'"'
		
	}
	else if "`ind1'"=="" & "`ind2'"=="" {
		local equs0 `"`equs0' (cuteq1_1: `y1' = , `offo1'  )"'
		gettoken drop thr10 : thr10
		local equs0 `"`equs0' `thr10'"'
		local equs0 `"`equs0' (cuteq2_1: `y2' = , `offo2'  )"'
		gettoken drop thr20 : thr20
		local equs0 `"`equs0' `thr20'"'
	}
	else if ("`ind1'"=="") & ("`ind2'"!="") {
		local equs0 `"`equs0' (`dep2n': `y2' = `ind2', `offo2'  nocons  )"'
		local equs0 `"`equs0' (cuteq1_1: `y1' = , `offo1'  )"'
		gettoken drop thr10 : thr10
		local equs0 `"`equs0' `thr10' `thr20'"'
	}
	else if ("`ind1'"!="") & ("`ind2'"=="") {
		local equs0 `"`equs0' (`dep1n': `y1' = `ind1', `offo1'  nocons  )"'
		local equs0 `"`equs0' `thr10'"'
	    local equs0 `"`equs0' (cuteq2_1: `y2' = , `offo2'  )"'
		gettoken drop thr20 : thr20
		local equs0 `"`equs0' `thr20'"'
	}
	
	
	
	if "`copula'" != "indep" {
		local equs0 `"`equs0'  /depend"'	/* transformed dependency parameter
										   is always called "depend" */
	}
	
	
//	mixture parameters

       if "`mixture'" == "equal" | "`mixture'" == "mix1" {	
		local equs0 `"`equs0'  /pu1 /mu2 /su2"'
	}	
       if "`mixture'" == "mix2" {	
		local equs0 `"`equs0'  /pv1 /mv2 /sv2"'
	}	
	if "`mixture'" == "both" {	
		local equs0 `"`equs0'  /pu1 /mu2 /su2 /pv1 /mv2 /sv2 "'
	}	
	


//	Call ML here, waldtest restricts all coefficients in the first two equations to zero
		`log' ml model lf1 bicop_lf1 `equs0'	if `touse' `wgt', /// 
				noconstant collinear missing maximize difficult `robust' `mlopts' ///
				init(`"`from'"') search(`search') repeat(`repeat') constraints(`constraints') `clopt' waldtest(-2) //gradient
			
// Wald test of common coefficients
	capture test[`dep1n' = `dep2n'], common	
	if _rc == 0 {   // there are common variables
		local wdf = r(df)
		local wchi2 = string(r(chi2),"%10.3f")
		local wp = string(r(p),"%4.3f")
		global bc_waldtest = `"Wald test of equality of coefficients chi2(df = `wdf')= `wchi2' [p-value=`wp']"'
	}
	else {
		global bc_waldtest = ""    // there are no variables in common			
	}
	
// Wald test of independence for Gaussian, Clayton and Frank
	quietly {
		if "$bc_copu"=="gaussian" {
			capture testnl tanh([depend]:_cons) = 0
			if _rc == 0 { 	
				scalar wdfi = r(df)
				scalar wchi2i = round(r(chi2),.001)
				scalar wpi = round(r(p),.001)
				local wdfib = r(df)
				local wchi2ib = string(r(chi2),"%10.3f")
				local wpib = string(r(p),"%4.3f")
				global bc_waldindp = `"Wald test of independence chi2(df = `wdfib')= `wchi2ib' [p-value=`wpib']"'
				ereturn scalar df_ind = wdfi 
				ereturn scalar chi2_ind = wchi2i
				ereturn scalar p_ind = wpi
			}
			else {
				global bc_waldindp = ""
			}
		}
		else if "$bc_copu"=="frank" {
			capture test [depend]:_cons
			if _rc == 0 {
				scalar wdfi = r(df)
				scalar wchi2i = round(r(chi2),.001)
				scalar wpi = round(r(p),.001)
				local wdfib = r(df)
				local wchi2ib = string(r(chi2),"%10.3f")
				local wpib = string(r(p),"%4.3f")
				global bc_waldindp = `"Wald test of independence chi2(df = `wdfib')= `wchi2ib' [p-value=`wpib']"'
				ereturn scalar df_ind = wdfi 
				ereturn scalar chi2_ind = wchi2i
				ereturn scalar p_ind = wpi
			}
			else {
				global bc_waldindp = ""
			}	
		
		} 
		else if "$bc_copu"=="clayton" {
			capture testnl exp([depend]:_cons) = 0
			if _rc == 0 {
				scalar wdfi = r(df)
				scalar wchi2i = round(r(chi2),.001)
				scalar wpi = round(r(p),.001)
				local wdfib = r(df)
				local wchi2ib = string(r(chi2),"%10.3f")
				local wpib = string(r(p),"%4.3f")
				global bc_waldindp = `"Wald test of independence chi2(df = `wdfib')= `wchi2ib' [p-value=`wpib']"'
				ereturn scalar df_ind = wdfi 
				ereturn scalar chi2_ind = wchi2i
				ereturn scalar p_ind = wpi
			}
			else {
				global bc_waldindp = ""
			}	
		} 
		else {
			global bc_waldindp = ""    			
		}
	}

		
				

/************************************************************************************/
	// Return extra results
	ereturn local cmd bicop
	ereturn local cmdline `cmdline'
	ereturn local copula `copula'
	ereturn local mixture `mixture'	
	ereturn local depvar "`dep1' `dep2'"	
	ereturn local indepvars1 `ind1'	
	ereturn local indepvars2 `ind2'	
	local title /// 
			"Generalized bivariate ordinal regression model (copula: `copula', mixture: `mixture')"
	ereturn local title `title'
	
	capture confirm variable bc_offo1
	if !_rc {
		ereturn local offset1 `offset1'
		qui drop bc_offo1 
	}
	capture confirm variable bc_offo2
	if !_rc {
		ereturn local offset2 `offset2'
		qui drop bc_offo2 
	}
	
	ereturn local predict "bicop_p"  
	
	
	if "$bc_copu"=="indep" {
		local adjind = 1
	}
	else {
		local adjind = 0
	}


	if "`mixture'" == "none" { 
		ereturn scalar k_aux = $bc_Nthr1 +$bc_Nthr2+ 1 -`adjind'
	}
    if "`mixture'" == "equal" |"`mixture'" == "mix1" | "`mixture'" == "mix2" {	
		ereturn scalar k_aux = $bc_Nthr1 +$bc_Nthr2+ 1 + 3 - `adjind'
	}	
    if "`mixture'" == "both" {	
		ereturn scalar k_aux = $bc_Nthr1 +$bc_Nthr2+ 1 + 6 - `adjind'
	}	
	
	ereturn scalar ll_c = $ll_ind  //loglikelihood for the comparison model
	Replay, `level'
end


program define rmTS, rclass

	local tsnm = cond( match("`0'", "*.*"),  		/*
			*/ substr("`0'", 			/*
			*/	  (index("`0'",".")+1),.),     	/*
			*/ "`0'")
	return local rmTS `tsnm'
end



program output_line
args vname estp sep zval pval zci

display as text %12s "`vname'" " {c |}" /*
*/ as result /*
*/ _col(17) %9.0g `estp' /*
*/ _col(28) %9.0g `sep' 
end
