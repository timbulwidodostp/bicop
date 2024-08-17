*! version 1.02 16 December 2015
*! version 1.01 18 September 2015
*! version 1.0 6 May 2015
*! Authors: Monica Hernandez and Steve Pudney

/***********************************************************/
/*      Generalized bivariate ordinal regression model     */
/*               Predictions                               */
/***********************************************************/


program define bicop_p
	version 13.1
	syntax [anything] [if] [in], [, SCores * ]
	
	
	if ("`e(cmd)'" != "bicop") {
			error 301
			display in red  "bicop was not the last command"
		}	
	/*     If scores option: calculate scores and exit  */	
	if `"`scores'"' != "" {
		ml score `0'
		exit
	}
	
	
	/* place command-unique options in local myopts  */
	
	local myopts PR XB1 XB2 STDP1 STDP2 PCOND1 PCOND2
	local myopts `myopts'  Outcome(string)
	/* call _pred_se */
	_pred_se "`myopts'" `0'
	
	if (`s(done)') exit
	local vtyp  `s(typ)'
	local varn `s(varn)'
	local 0 `"`s(rest)'"'	
	
	/*  parse the syntax */
	syntax [if] [in] [, `myopts' noOFFset]
	
	/*concatenate switch options together */
	
	local type `pr' `xb1' `xb2' `stdp1' `stdp2' 
	local type `type' `pcond1' `pcond2'
	
	
	tokenize `e(depvar)'
	local dep1 `1'
	local dep2 `2'
	
	local mixture `e(mixture)'
	
	
	tempvar xbee1 xbee2 tmpv tmpv1 bott1 topp1 bott2 topp2 res1 res2
	tempname depend theta pu1 pu2 pv1 pv2 mu1 mu2 mv1 mv2 su1 su2 sv1 sv2 
	
	if (missing("`type'") & missing("`outcome'")) | (("`type'"== "pr") & missing("`outcome'")) {
		local outcome "#1, #1"
		local type "pr"
		noisily display as text "(option pr outcome(#1, #1) assumed)"
	}
	else if  (missing("`type'")){
		local type "pr"	
	}
	else if ((("`type'"== "pcond1") | ("`type'"== "pcond2") ) & missing("`outcome'")) {
		local outcome "#1, #1"
		noisily display as text "(option outcome(#1, #1) assumed)"
	}
	
	
	tempname off1 off2
	tempvar xbee1 xbee2 offset1 offset2
	local a = e(offset1)
	qui gen double  `offset1' = `a'
	qui summ `offset1'
	local a = r(N)	
	local bc_kx1 : word count `e(ind1)'
	
	if (`a' >0) & (`bc_kx1' == 0) {
		qui gen double `xbee1' = `offset1'	
	}
	else {
		qui _predict double `xbee1' `if' `in', eq(#1) `offset'
	}
	
	local a = e(offset2)
	qui gen double `offset2' = `a'
	qui summ `offset2'
	local a = r(N)
	local bc_kx2 : word count `e(ind2)'

	if (`a' >0) & (`bc_kx2' == 0) {
		qui gen double `xbee2' = `offset2'		
	}
	else {
		qui _predict double `xbee2' `if' `in', eq(#2) `offset'
	}	
	
	/*  Mark the prediction sample */
	
	marksample touse
	
	//adapted from bioprobit
	
	if ~missing("`outcome'") {
		// get cuttoff points
		_parse comma out1 out2: outcome
		local out2 : subinstr local out2 "," "", all
		local out2 = trim("`out2'")
		forvalues i = 1/2 {
			tempvar sy
			egen `sy' = group(`dep`i'')
			capture confirm numeric variable `dep`i''	
			scalar rc = _rc
			if substr(`"`out`i''"',1,1)=="#" {
				local out = substr(`"`out`i''"',2,.)
				Chk confirm integer number `out'
				Chk assert `out' >= 1
				quietly tabulate `dep`i'' if `touse'
				if (`out'>r(r)) {
					noisily display as error "there is no outcome #`out'" _n ///
										 "there are only `r(r)' categories in `dep`i''"
					exit 111
				}
				else local c =`out'
			}
			else {
				if (rc) local out`i' "`"`out`i''"'"
				if missing(`out`i'') {
					local c=-99
					if ("`type'" == "pcond1") | ("`type'" == "pcond1") {
						display in red "Option pcond requires a pair of numerical outcomes."
						exit 498
					}
				}
				else {
					Chk confirm number `out`i''
					summarize `sy' if `touse' & `dep`i''==`out`i'', meanonly 
					local c = r(mean)
					if missing(`c') {
						noisily display as error "`dep`i'' never takes value `out`i''"
						exit 111
					}
				}
			}
			if (`c' == -99) {
				local top`i'   = maxfloat()
				local bot`i' = minfloat()
			}
			else {
				quietly summarize `sy' if `touse'
				if (`c'==r(max)) local top`i' = maxfloat() /** top threshold = inf if highest group */
				else             local top`i' = _b[cuteq`i'_`c':_cons] /**  top threshold **/
				local --c
				if (`c'==0) local bot`i' = minfloat()     /*** bottom threshold = -inf if lowest group **/
				else 		local bot`i' = _b[cuteq`i'_`c':_cons]  /*** bottom threshold ***/
			}
			local c`i' = `c'
			
		} // end forvalues
	

	//	dependency parameter theta	
		
		local copu = e(copula)
		if "`copu'"=="indep" {
			local depend = 0
			scalar `theta'= 0
		}
		else {
			local depend = _b[depend:_cons]
		}
	
	
		if "`copu'"=="gaussian" {
			scalar `theta'=tanh(`depend')
		}
		if "`copu'"=="frank" {
			scalar `theta'=`depend'
		} 
		if "`copu'"=="clayton" {
			scalar `theta'=exp(`depend')
		} 
		if "`copu'"=="joe"|"`copu'"=="gumbel" {
			scalar `theta'=1+exp(`depend')
		} 


//	mixing distribution
		local mix = e(mixture)
		if "`mix'"=="equal" | "`mix'"=="mix1" | "`mix'"=="both" {
		
			scalar `pu1' = _b[pu1:_cons]
			scalar `pu1'= exp(`pu1')/(1+exp(`pu1'))
			scalar `pu2'=1-`pu1'
		
			scalar `mu2' = _b[mu2:_cons]
			scalar `mu1' = -`pu2'*`mu2'/`pu1'
		
			scalar `su2' = _b[su2:_cons]
			scalar `su2' = `su2'^2
			scalar `su1' = (1-`pu2'*(`su2'+(`mu2'^2)))/`pu1'  -  (`mu1'^2)
			if `su1'>=0 scalar `su1'=sqrt(`su1')	
			scalar `su2'=sqrt(`su2')
			if "`mix'"=="equal" {
				scalar `pv1'=`pu1'
				scalar `pv2'=`pu2'
				scalar `mv1'=`mu1'
				scalar `mv2'=`mu2'
				scalar `sv1'=`su1'
				scalar `sv2'=`su2'
			}
		}
		if "`mix'"=="both" | "`mix'"=="mix2" {
		
			scalar `pv1' = _b[pv1:_cons]
			scalar `pv1'= exp(`pv1')/(1+exp(`pv1'))
			scalar `pv2'=1-`pv1'
		
			scalar `mv2' = _b[mv2:_cons]
			scalar `mv1' = -`pv2'*`mv2'/`pv1'
		
			scalar `sv2' = _b[sv2:_cons]
			scalar `sv2' = `sv2'^2
			scalar `sv1' = (1-`pv2'*(`sv2'+(`mv2'^2)))/`pv1'  -  (`mv1'^2)
			if `sv1'>=0 scalar `sv1'=sqrt(`sv1')	
			scalar `sv2'=sqrt(`sv2')
		}
  

//	construct marginals and/or joints
		qui capture gen double `tmpv'=.
		qui capture gen double `tmpv1'=.
		qui capture gen double `bott1'=.
		qui capture gen double `bott2'=.
		qui capture gen double `topp1'=.
		qui capture gen double `topp2'=.
		qui capture gen double `res1'=.
	
// run it twice if option pcond to get the joint and the marginal
		if ("`type'" == "pcond1") | ("`type'" == "pcond2") {
			local indx = 2	
			qui capture gen double `res2'=.
		}
		else {
			local indx = 1		
		}
		forvalues j = 1(1)`indx' {
			if `j' == 2 {
				if ("`type'" == "pcond1") {
					local bot1 = minfloat() 
					local top1 = maxfloat()
				
				}
				
				else if ("`type'" == "pcond2") {
					local bot2 = minfloat() 
					local top2 = maxfloat()				
				}
			
			}
		
			if "`mixture'"=="mix1"|"`mixture'"=="both"|"`mixture'"=="equal" {
			  qui{
				replace `tmpv'=normal((`bot1'-`xbee1'-`mu1')/`su1')
				replace `tmpv'=0 if `bot1'==minfloat() 
				replace `tmpv'=1 if `bot1'==maxfloat()
				replace `tmpv1'=normal((`bot1'-`xbee1'-`mu2')/`su2')
				replace `tmpv1'=0 if `bot1'==minfloat() 
				replace `tmpv1'=1 if `bot1'==maxfloat()
				replace `bott1'=`pu1'*`tmpv' + `pu2'*`tmpv1'
		
			replace `tmpv'=normal((`top1'-`xbee1'-`mu1')/`su1')
				replace `tmpv'=0 if `top1'==minfloat() 
				replace `tmpv'=1 if `top1'==maxfloat()
				replace `tmpv1'=normal((`top1'-`xbee1'-`mu2')/`su2')
				replace `tmpv1'=0 if `top1'==minfloat() 
				replace `tmpv1'=1 if `top1'==maxfloat()
				replace `topp1'=`pu1'*`tmpv'+ `pu2'*`tmpv1'
			  }
			}
			if "`mixture'"=="mix2"|"`mixture'"=="both"|"`mixture'"=="equal" {
			
			  qui{
				replace `tmpv'=normal((`bot2'-`xbee2'-`mv1')/`sv1')
				replace `tmpv'=0 if `bot2'==minfloat() 
				replace `tmpv'=1 if `bot2'==maxfloat()
				replace `tmpv1'=normal((`bot2'-`xbee2'-`mv2')/`sv2')
				replace `tmpv1'=0 if `bot2'==minfloat() 
				replace `tmpv1'=1 if `bot2'==maxfloat()
				replace `bott2'=`pv1'*`tmpv' + `pv2'*`tmpv1'
		
				replace `tmpv'=normal((`top2'-`xbee2'-`mv1')/`sv1')
				replace `tmpv'=0 if `top2'==minfloat() 
				replace `tmpv'=1 if `top2'==maxfloat()
				replace `tmpv1'=normal((`top2'-`xbee2'-`mv2')/`sv2')
				replace `tmpv1'=0 if `top2'==minfloat() 
				replace `tmpv1'=1 if `top2'==maxfloat()
				replace `topp2'=`pv1'*`tmpv' + `pv2'*`tmpv1'
				
				
			  }
			}
			if "`mixture'"=="none" | "`mixture'"=="mix2" {
			  qui{
				replace `tmpv'=normal(`bot1'-`xbee1')
				replace `tmpv'=0 if `bot1'==minfloat() 
				replace `tmpv'=1 if `bot1'==maxfloat()
				replace `bott1'=`tmpv'
		
				replace `tmpv'=normal(`top1'-`xbee1')
				replace `tmpv'=0 if `top1'==minfloat() 
				replace `tmpv'=1 if `top1'==maxfloat()
				replace `topp1'=`tmpv'
			  }
			}	
			if "`mixture'"=="none" | "`mixture'"=="mix1" {
			  qui{
				replace `tmpv'=normal(`bot2'-`xbee2')
				replace `tmpv'=0 if `bot2'==minfloat() 
				replace `tmpv'=1 if `bot2'==maxfloat()
				replace `bott2'=`tmpv'
		
				replace `tmpv'=normal(`top2'-`xbee2')
				replace `tmpv'=0 if `top2'==minfloat() 
				replace `tmpv'=1 if `top2'==maxfloat()
				replace `topp2'=`tmpv'
			  }
			}
			
	
//	evaluate copula likelihood
			qui replace `tmpv'=.

			copu `topp1' `topp2' `theta' `tmpv'
			qui replace  `res`j''=`tmpv'

			copu `bott1' `topp2' `theta' `tmpv'
			qui replace `res`j''=`res`j''-`tmpv'

			copu `topp1' `bott2' `theta' `tmpv'
			qui replace `res`j''=`res`j''-`tmpv'

			copu `bott1' `bott2' `theta' `tmpv'
			qui replace `res`j''=`res`j''+`tmpv'
			if ("`type'" == "pr") {
	
				generate `vtyp' `varn' = `res`j''	
				if (`c1' != -99) &  (`c2' != -99) {
					qui tab `dep1', matrow(matc1)
					qui tab `dep2', matrow(matc2)
					local c1 = el(matc1,`c1'+1,1)
					local c2 = el(matc2,`c2'+1,1)
					label variable `varn' "Pr(`dep1'=`c1', `dep2'=`c2')" 	
					
				}
				else if `c1' == -99 {
					qui tab `dep2', matrow(matc2)
					local c2 = el(matc2,`c2'+1,1)
					label variable `varn' "Pr(`dep2'=`c2')" 	
									
				}
				
				else if `c2' == -99 {
					qui tab `dep1', matrow(matc1)
					local c1 = el(matc1,`c1'+1,1)
					label variable `varn' "Pr(`dep1'=`c1')" 	
								
				}
				exit
			}
			
			else if (("`type'"== "pcond1") | ("`type'"== "pcond2")) & (`j' == 2) {
				generate `vtyp' `varn' = `res1' / `res2'
				if ("`type'"== "pcond1") {
					qui tab `dep1', matrow(matc1)
					local c1 = el(matc1,`c1'+1,1)
					label variable `varn' "Pr(`dep1'=`c1' | `dep2')"
				}					
				else {
					qui tab `dep2', matrow(matc2)
					local c2 = el(matc2,`c2'+1,1)
					label variable `varn' "Pr(`dep2'=`c2' | `dep1')"
				}
				exit
			}
	

			
			}  
	
	} // outcome
	
	
	
	
	    /* xb1  */
	if ("`type'" == "xb1") {
		generate `vtyp' `varn' = `xbee1'
		label variable `varn' "linear prediction of `dep1'"
		exit
	}
	
		 /* xb2  */
	if ("`type'" == "xb2") {
		generate `vtyp' `varn' = `xbee2'
		label variable `varn' "linear prediction of `dep2'"
		exit
	}
	  /* stdp1 */
	if ("`type'" == "stdp1") {
		qui _predict `vtyp' `varn', stdp eq(#1) `offset', if `touse'
		label variable `varn' "S.E. of prediction of `dep1'"
		exit
	}
	
	  /* stdp2 */
	if ("`type'" == "stdp2") {
		qui _predict `vtyp' `varn', stdp eq(#2) `offset', if `touse'
		label variable `varn' "S.E. of prediction of `dep2'"
		exit
	}
	
	
	
	error 198  /** error in the options - not picked up by syntax, i.e. combinations of options**/
	
end


//stolen from ologit_p.ado 
program define Chk
	capture `0'
	if _rc {
		noisily display as error "outcome() must be pair of either values of `e(depvar)',"  _n "or #1, #2, ..."
		exit 111
	}
end

program define copu
version 13
args u v theta cop
local copul = e(copula)
quietly {
	if "$bc_copu"=="indep" {
		 replace `cop'=`u'*`v'
	}
	if "`copul'"=="gaussian" {
		replace `cop'=binorm(invnormal(`u'),invnormal(`v'),`theta') ///
										if `u'>0&`u'<1&`v'>0&`v'<1
	}
	if "`copul'"=="frank" {
		if `theta'==0 replace `cop'=`u'*`v'
		else {
			tempname eta
			scalar `eta'=1-exp(-`theta')
			replace `cop'=`eta'-(1-exp(-`theta'*`u'))*(1-exp(-`theta'*`v')) ///
										if `u'>0&`u'<1&`v'>0&`v'<1
			replace `cop'=-(ln(`cop'/`eta'))/`theta' if `u'>0&`u'<1&`v'>0&`v'<1
		}
	}
	if "`copul'"=="clayton" {
		if `theta'==0 replace `cop'=`u'*`v'
		else {
			replace `cop'=(`u'^(-`theta')) + (`v'^(-`theta')) - 1 /// 
										if `u'>0&`u'<1&`v'>0&`v'<1
			replace `cop'=(max(0,`cop'))^(-1/`theta') if `u'>0&`u'<1&`v'>0&`v'<1
		}
	}
	if "`copul'"=="gumbel" {
		replace `cop'=exp(  - ((-ln(`u'))^`theta'+(-ln(`v'))^`theta') /// 
						^(1/`theta')  ) if `u'>0&`u'<1&`v'>0&`v'<1
	}
	if "`copul'"=="joe" {
		replace `cop'=1 - (((1-`u')^`theta') + ((1-`v')^`theta') ///
					-((1-`u')^`theta')*((1-`v')^`theta') )^(1/`theta') /// 
										if `u'>0&`u'<1&`v'>0&`v'<1
	}
	replace `cop'=0 if `u'==0|`v'==0
	replace `cop'=1 if `u'==1&`v'==1
	replace `cop'=`u' if `u'>0&`u'<1&`v'==1
	replace `cop'=`v' if `u'==1&`v'>0&`v'<1
}
end
