*! version 1.02 16 December 2015
*! version 1.01 18 September 2015
*! version 1.0 6 May 2015
*! Authors: Monica Hernandez and Steve Pudney

/***********************************************************/
/* likelihood function  for the generalized bivariate 
            ordinal regression model  */
/***********************************************************/
program define bicop_lf1
	version 13.1
	
	if ($bc_kx1 != 0) & ($bc_kx2 != 0){
		local nn = 2
	}
	else if ($bc_kx1 == 0) & ($bc_kx2 == 0){
		local nn = 0
	}
	else {
		local nn = 1
	}
	
	if "$bc_mix"=="none" {
		local n_grad   = `nn' + $bc_Nthr1 + $bc_Nthr2 +1
	}	
	if "$bc_mix"=="equal" | "$bc_mix"=="mix1" | "$bc_mix"=="mix2" {
		local n_grad   = `nn' + $bc_Nthr1 + $bc_Nthr2 +4
	}
	if "$bc_mix"=="both" {
		local n_grad   = `nn' + $bc_Nthr1 + $bc_Nthr2 +7
	}	
	if "$bc_copu"=="indep" {
		local n_grad = `n_grad' -1
	}
	
 	forvalues i = 1/`n_grad' { 
		local grad "`grad' g`i'" 
	}

	args todo b lnfj `grad'
	
	
	
	
// equations (y1=xb1) (y2=xb2) Thres for 1 thresh for 2 depend p1 mu1 su2 p2 mv1 sv2

	tempvar xb1 xb2 tmpv tmpv1 bot1 top1 bot2 top2 lik bbot1 ttop1 bbot2 ttop2 d1 d2 d3 d4 ///
			ztop11 ztop12 zbot11 zbot12 ztop21 ztop22 zbot21 zbot22 
			
	tempname depend pu1 pu2 pv1 pv2 mu1 mu2 mv1 mv2 su1 su2 sv1 sv2 thr1 /// 
		thr2 theta tmp tmp1 tmp2 tmp3 tmp4 i nxb b31 b32



quietly{

//	construct x*beta1 and x*beta2	
	
	
	if ($bc_kx1 != 0) & ($bc_kx2 != 0) {
		mleval `xb1' = `b', eq(1)
		mleval `xb2' = `b', eq(2)
		local nxb = 2
	}
	else if ($bc_kx1 == 0) & ($bc_kx2 == 0){
		capture confirm variable bc_offo1
		if !_rc {
			gen double `xb1' = bc_offo1
		}
		else {
			gen double `xb1' = 0
		}
		capture confirm variable bc_offo2
		if !_rc {
			gen double `xb2' = bc_offo2
		}
		else {
			gen double `xb2' = 0
		}			
		local nxb = 0
	}
	
	else if ($bc_kx1 != 0) & ($bc_kx2 == 0) {
		mleval `xb1' = `b', eq(1)
		capture confirm variable bc_offo2
		if !_rc {
			gen double `xb2' = bc_offo2
		}
		else {
			gen double `xb2' = 0
		}
		local nxb = 1
	}	
	else if ($bc_kx1 == 0) & ($bc_kx2 != 0)  {
		capture confirm variable bc_offo1
		if !_rc {
			gen double `xb1' = bc_offo1
		}
		else {
			gen double `xb1' = 0
		}
		mleval `xb2' = `b', eq(1)		
		local nxb = 1	
	}
		
	
//	vectors of threshold parameters 
	mat `thr1' = J($bc_Nthr1+ 2, 1, minfloat())
	forvalues i=1/$bc_Nthr1 {
		local tmp1 = `nxb'+`i'
		if (`i' == 1) & ($bc_kx1 == 0) {
			mleval `tmp3' = `b', eq(`tmp1')
			capture confirm variable bc_offo1
			if !_rc {
				replace `tmp3' = `tmp3' - bc_offo1
			}	
			summ `tmp3' 
			mat `thr1'[`i'+1, 1] = r(mean)	
		}
		else {
			mleval `tmp' = `b', eq(`tmp1') scalar
			mat `thr1'[`i'+1, 1] = `tmp'
		}	
		
	}	
	mat `thr1'[$bc_Nthr1+2, 1] = maxfloat()
	
	mat `thr2' = J($bc_Nthr2+2, 1, minfloat())
	forvalues i=1/$bc_Nthr2 {
		local tmp1 = `nxb'+ $bc_Nthr1 +`i'
		if (`i' == 1) & ($bc_kx2 == 0) {
			mleval `tmp4' = `b', eq(`tmp1')
			capture confirm variable bc_offo2
			if !_rc {
				replace `tmp4' = `tmp4' - bc_offo2
			}	
			summ `tmp4' 
			mat `thr2'[`i'+1, 1] = r(mean)	
		}
		else {
			mleval `tmp' = `b', eq(`tmp1') scalar
			mat `thr2'[`i'+1, 1] = `tmp'
		}	
	}
	mat `thr2'[$bc_Nthr2+2, 1] = maxfloat()
	
	if ($bc_kx1 == 0) & ($bc_kx2 != 0)  {
		// for this case we need to switch the order of the dependent variables
		gen double `bot1'=`thr1'[$ML_y2 ,1]		//	construct variables containing 
		gen double `top1'=`thr1'[$ML_y2+1 ,1]	//	relevant thresholds
		gen double `bot2'=`thr2'[$ML_y1 ,1]
		gen double `top2'=`thr2'[$ML_y1+1 ,1]
		
		
		
		gen double `bbot1'=`thr1'[$ML_y2 ,1]		//	Same as before for the derivatives
		gen double `ttop1'=`thr1'[$ML_y2+1 ,1]	
		gen double `bbot2'=`thr2'[$ML_y1 ,1]
		gen double `ttop2'=`thr2'[$ML_y1+1 ,1]
	
	}
	else {
		gen double `bot1'=`thr1'[$ML_y1 ,1]		//	construct variables containing 
		gen double `top1'=`thr1'[$ML_y1+1 ,1]	//	relevant thresholds
		gen double `bot2'=`thr2'[$ML_y2 ,1]
		gen double `top2'=`thr2'[$ML_y2+1 ,1]
		
	
		
		gen double `bbot1'=`thr1'[$ML_y1 ,1]		//	Same as before for the derivatives
		gen double `ttop1'=`thr1'[$ML_y1+1 ,1]	
		gen double `bbot2'=`thr2'[$ML_y2 ,1]
		gen double `ttop2'=`thr2'[$ML_y2+1 ,1]
		
	}
//	dependency parameter theta
	if "$bc_copu"=="indep" {
		local adjind = 1 //no dependency parameter if independent copula
		scalar `theta'= 0
	}
	else {
		local adjind = 0
	}
	
	if "$bc_copu"!="indep" {
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 +1
	
		mleval `depend' = `b', eq(`tmp1') scalar
	
		if "$bc_copu"=="gaussian" {
			scalar `theta'=tanh(`depend')
		}
		if "$bc_copu"=="frank" {
			scalar `theta'=`depend'
		} 
		if "$bc_copu"=="clayton" {
			scalar `theta'=exp(`depend')
		} 
		if "$bc_copu"=="joe"|"$bc_copu"=="gumbel" {
			scalar `theta'=1+exp(`depend')
		} 
	
	}

//	mixing distribution
	if "$bc_mix"=="both"|"$bc_mix"=="equal"|"$bc_mix"=="mix1" {
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 2 -`adjind'
		mleval `pu1' = `b', eq(`tmp1') scalar
		scalar `pu1'= exp(`pu1')/(1+exp(`pu1'))
		scalar `pu2'=1-`pu1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 3 -`adjind'
		mleval `mu2' = `b', eq(`tmp1') scalar
		
		scalar `mu1' = -`pu2'*`mu2'/`pu1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 4 -`adjind'
		mleval `su2' = `b', eq(`tmp1') scalar
		scalar `b31' = `su2'
		
		scalar `su2' = `su2'^2
		scalar `su1' = (1-`pu2'*(`su2'+(`mu2'^2)))/`pu1'  -  (`mu1'^2)
		if `su1'>=0 scalar `su1'=sqrt(`su1')	// penalize -ve variances later
		scalar `su2'=sqrt(`su2')
		
		if "$bc_mix"=="equal" {
			scalar `b32' = `b31'
			scalar `pv1'=`pu1'
			scalar `pv2'=`pu2'
			scalar `mv1'=`mu1'
			scalar `mv2'=`mu2'
			scalar `sv1'=`su1'
			scalar `sv2'=`su2'
		}
	}
	if "$bc_mix"=="mix2" {
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 2 -`adjind'
		mleval `pv1' = `b', eq(`tmp1') scalar
		scalar `pv1'= exp(`pv1')/(1+exp(`pv1'))
		scalar `pv2'=1-`pv1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 3 -`adjind'
		mleval `mv2' = `b', eq(`tmp1') scalar
		scalar `mv1' = -`pv2'*`mv2'/`pv1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 4 -`adjind'
		mleval `sv2' = `b', eq(`tmp1') scalar
		scalar `b32' = `sv2'
		scalar `sv2' = `sv2'^2
		scalar `sv1' = (1-`pv2'*(`sv2'+(`mv2'^2)))/`pv1'  -  (`mv1'^2)
		if `sv1'>=0 scalar `sv1'=sqrt(`sv1')	// penalize -ve variances later
		scalar `sv2'=sqrt(`sv2')
	}
	if "$bc_mix"=="both" {
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 5 -`adjind'
		mleval `pv1' = `b', eq(`tmp1') scalar
		scalar `pv1'= exp(`pv1')/(1+exp(`pv1'))
		scalar `pv2'=1-`pv1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 6 -`adjind'
		mleval `mv2' = `b', eq(`tmp1') scalar
		scalar `mv1' = -`pv2'*`mv2'/`pv1'
		local tmp1 = `nxb'+ $bc_Nthr1 + $bc_Nthr2 + 7 -`adjind'
		mleval `sv2' = `b', eq(`tmp1') scalar
		scalar `b32' = `sv2'
		scalar `sv2' = `sv2'^2
		scalar `sv1' = (1-`pv2'*(`sv2'+(`mv2'^2)))/`pv1'  -  (`mv1'^2)
		if `sv1'>=0 scalar `sv1'=sqrt(`sv1')	// penalize -ve variances later
		scalar `sv2'=sqrt(`sv2')
	}
	
//	construct marginals
	capture gen double `tmpv'=.
	capture gen double `tmpv1'=.
	if "$bc_mix"=="mix1"|"$bc_mix"=="both"|"$bc_mix"=="equal" {
		
		replace `tmpv'=normal((`bot1'-`xb1'-`mu1')/`su1')
		replace `tmpv'=0 if `bot1'==minfloat() 
		replace `tmpv'=1 if `bot1'==maxfloat()
		gen double `zbot11' = (`bot1'-`xb1'-`mu1')/`su1'
		replace `tmpv1'=normal((`bot1'-`xb1'-`mu2')/`su2')
		replace `tmpv1'=0 if `bot1'==minfloat() 
		replace `tmpv1'=1 if `bot1'==maxfloat()
		gen double `zbot12' = (`bot1'-`xb1'-`mu2')/`su2'
		replace `bot1'=`pu1'*`tmpv' + `pu2'*`tmpv1'
		
		replace `tmpv'=normal((`top1'-`xb1'-`mu1')/`su1')
		replace `tmpv'=0 if `top1'==minfloat() 
		replace `tmpv'=1 if `top1'==maxfloat()
		gen double `ztop11' = (`top1'-`xb1'-`mu1')/`su1'
		replace `tmpv1'=normal((`top1'-`xb1'-`mu2')/`su2')
		replace `tmpv1'=0 if `top1'==minfloat() 
		replace `tmpv1'=1 if `top1'==maxfloat()
		gen double `ztop12' = (`top1'-`xb1'-`mu2')/`su2'
		replace `top1'=`pu1'*`tmpv'+ `pu2'*`tmpv1'
	}
	if "$bc_mix"=="mix2"|"$bc_mix"=="both"|"$bc_mix"=="equal" {
		replace `tmpv'=normal((`bot2'-`xb2'-`mv1')/`sv1')
		replace `tmpv'=0 if `bot2'==minfloat() 
		replace `tmpv'=1 if `bot2'==maxfloat()
		gen double `zbot21' = (`bot2'-`xb2'-`mv1')/`sv1'
		replace `tmpv1'=normal((`bot2'-`xb2'-`mv2')/`sv2')
		replace `tmpv1'=0 if `bot2'==minfloat() 
		replace `tmpv1'=1 if `bot2'==maxfloat()
		gen double `zbot22' = (`bot2'-`xb2'-`mv2')/`sv2'
		replace `bot2'=`pv1'*`tmpv' + `pv2'*`tmpv1'
		
		replace `tmpv'=normal((`top2'-`xb2'-`mv1')/`sv1')
		replace `tmpv'=0 if `top2'==minfloat() 
		replace `tmpv'=1 if `top2'==maxfloat()
		gen double `ztop21' = (`top2'-`xb2'-`mv1')/`sv1'
		replace `tmpv1'=normal((`top2'-`xb2'-`mv2')/`sv2')
		replace `tmpv1'=0 if `top2'==minfloat() 
		replace `tmpv1'=1 if `top2'==maxfloat()
		gen double `ztop22' = (`top2'-`xb2'-`mv2')/`sv2'
		replace `top2'=`pv1'*`tmpv' + `pv2'*`tmpv1'
	}
	if "$bc_mix"=="none" | "$bc_mix"=="mix2" {
		replace `tmpv'=normal(`bot1'-`xb1')
		replace `tmpv'=0 if `bot1'==minfloat() 
		replace `tmpv'=1 if `bot1'==maxfloat()
		replace `bot1'=`tmpv'
		
		replace `tmpv'=normal(`top1'-`xb1')
		replace `tmpv'=0 if `top1'==minfloat() 
		replace `tmpv'=1 if `top1'==maxfloat()
		replace `top1'=`tmpv'
	}	
	if "$bc_mix"=="none" | "$bc_mix"=="mix1" {
		replace `tmpv'=normal(`bot2'-`xb2')
		replace `tmpv'=0 if `bot2'==minfloat() 
		replace `tmpv'=1 if `bot2'==maxfloat()
		replace `bot2'=`tmpv'
		
		replace `tmpv'=normal(`top2'-`xb2')
		replace `tmpv'=0 if `top2'==minfloat() 
		replace `tmpv'=1 if `top2'==maxfloat()
		replace `top2'=`tmpv'
	}
	
//	evaluate copula likelihood
	replace `tmpv'=.

	copu `top1' `top2' `theta' `tmpv'
	replace  `lnfj'=`tmpv'
	
	//noi di "lnf1" `lnfj'

	copu `bot1' `top2' `theta' `tmpv'
	replace `lnfj'=`lnfj'-`tmpv'
	//noi di  `lnfj'

	copu `top1' `bot2' `theta' `tmpv'
	replace `lnfj'=`lnfj'-`tmpv'
	//noi di  `lnfj'
	copu `bot1' `bot2' `theta' `tmpv'
	replace `lnfj'=`lnfj'+`tmpv'
	//noi di  `lnfj'

//	sort out invalid likelihoods
	replace `lnfj'=0.1D-60 if `lnfj'<=0		//	underflowing likelihood
	if "$bc_mix"=="both" | "$bc_mix"=="equal" {
		replace `lnfj'=0.1D-60 if `su1'<0|`sv1'<0  // invalid parameter value
	}
	if "$bc_mix"=="mix1" {
		replace `lnfj'=0.1D-60 if `su1'<0  // invalid parameter value
	}
	if "$bc_mix"=="mix2" {
		replace `lnfj'=0.1D-60 if `sv1'<0  // invalid parameter value
	}
//	log-likelihood
	capture gen double `lik' = `lnfj' //keep value for derivatives
	
	replace `lnfj'=ln(`lnfj') 
 
	

	if (`todo'==0 | `lnfj'>=.)  {
		exit // end of loglikelihood evaluation
	}
	
	
/***************************************************************/	
// gradients
/****************************************************************/

	forvalues i = 1/`n_grad' { 
		tempvar gg`i'
		gen double `gg`i'' = .
	}
	
	forvalues i = 1/2 {
		tempvar d`i'a d`i'b
		gen double `d`i'a' = .
		gen double `d`i'b' = .
	}
	
		
	
	if "$bc_mix"=="none" | "$bc_mix"=="mix2"   {
		//du/dxb1
		replace `d1a' = - normalden(`ttop1'-`xb1')
		replace `d1a' = 0 if (`ttop1' == minfloat())  | (`ttop1' == maxfloat())
		replace `d1b' = - normalden(`bbot1'-`xb1')
		replace `d1b' = 0 if (`bbot1' == minfloat()) | (`bbot1' == maxfloat()) 
	}
	if "$bc_mix"=="none" | "$bc_mix"=="mix1"   {
	
		// dv/dxb2
		replace `d2a' = - normalden(`ttop2'-`xb2')
		replace `d2a' = 0 if (`ttop2' == minfloat())  | (`ttop2' == maxfloat())
		replace `d2b' = - normalden(`bbot2'-`xb2')
		replace `d2b' = 0 if (`bbot2' == minfloat())  | (`bbot2' == maxfloat())
	}
	
	if "$bc_mix"=="mix1" | "$bc_mix"=="equal"  | "$bc_mix"=="both"  {
	
		//du/dxb1
		replace `d1a' = -`pu1' * normalden(`ztop11') /`su1' - `pu2' * normalden(`ztop12') /`su2'
		replace `d1a' = 0 if (`ttop1' == minfloat())  | (`ttop1' == maxfloat())
		replace `d1b' = -`pu1' * normalden(`zbot11') /`su1' - `pu2' * normalden(`zbot12') /`su2'
		replace `d1b' = 0 if (`bbot1' == minfloat()) | (`bbot1' == maxfloat()) 
	}
	
	if "$bc_mix"=="mix2" | "$bc_mix"=="equal"  | "$bc_mix"=="both"  {
	
		// dv/dxb2
		replace `d2a' = -`pv1' * normalden(`ztop21') /`sv1' - `pv2' * normalden(`ztop22') /`sv2'
		replace `d2a' = 0 if (`ttop2' == minfloat())  | (`ttop2' == maxfloat())
		replace `d2b' = -`pv1' * normalden(`zbot21') /`sv1' - `pv2' * normalden(`zbot22') /`sv2'
		replace `d2b' = 0 if (`bbot2' == minfloat())  | (`bbot2' == maxfloat())
	}
 

	// first derivative with respect to u and v	
	forvalues i = 11/14 { 
		tempvar tmp`i'
		gen double `tmp`i'' = . 
	}
		
	forvalues i = 21/24 { 
		tempvar tmp`i'
		gen double `tmp`i'' = . 
	}
	

	if `theta'==0 & (("$bc_copu"=="frank") | ("$bc_copu"=="clayton") | ("$bc_copu"=="indep")) {	
		replace `tmp11' =`top2'
		replace `tmp12' =`top2'
		replace `tmp13' =`bot2'
		replace `tmp14' =`bot2'
		
		replace `tmp21' =`top1'
		replace `tmp22' =`bot1'
		replace `tmp23' =`top1'
		replace `tmp24' =`bot1'		
	}
	else {
		dcopu `top1' `top2' `theta' "1" `tmp11'
		dcopu `bot1' `top2' `theta' "1" `tmp12'
		dcopu `top1' `bot2' `theta' "1" `tmp13'
		dcopu `bot1' `bot2' `theta' "1" `tmp14'
		
		dcopu `top1' `top2' `theta' "2" `tmp21'
		dcopu `bot1' `top2' `theta' "2" `tmp22'
		dcopu `top1' `bot2' `theta' "2" `tmp23'
		dcopu `bot1' `bot2' `theta' "2" `tmp24'
	}		
	if ($bc_kx1 != 0) & ($bc_kx2 != 0){
		replace `gg1' = `tmp11' *`d1a' - `tmp12' *`d1b' - `tmp13' *`d1a'+ `tmp14' *`d1b'
		replace `gg2' = `tmp21' *`d2a' - `tmp22' *`d2a' - `tmp23' *`d2b'+ `tmp24' *`d2b'
	}
	
	else if ($bc_kx1 == 0) & ($bc_kx2 != 0) {
		replace `gg1' = `tmp21' *`d2a' - `tmp22' *`d2a' - `tmp23' *`d2b'+ `tmp24' *`d2b'
	}
	else if ($bc_kx1 != 0) & ($bc_kx2 == 0) {
		replace `gg1' = `tmp11' *`d1a' - `tmp12' *`d1b' - `tmp13' *`d1a'+ `tmp14' *`d1b'
	}
	
	
					
	local n1 = `nn' +1
	local n2 = `nn' + $bc_Nthr1
		
	forvalues i = `n1'/`n2' {
			
		replace `d1a' = .
		replace `d1b' = .
		replace `tmpv' = 0
		replace `tmpv' = 1 if (`ttop1' == `thr1'[`i'-`nn'+1,1])
				
				
		replace `tmpv1' = 0
		replace `tmpv1' = 1 if (`bbot1' == `thr1'[`i'-`nn'+1,1])
				
				
		if "$bc_mix"=="none" | "$bc_mix"=="mix2"   {
			//du/dthr1
			replace `d1a' =  normalden(`ttop1'-`xb1') if `tmpv' == 1
			replace `d1a' = 0 if `tmpv' == 0
			replace `d1b' = normalden(`bbot1'-`xb1')  if `tmpv1' == 1
			replace `d1b' = 0  if `tmpv1' == 0
		}
	
		if "$bc_mix"=="mix1" | "$bc_mix"=="equal"  | "$bc_mix"=="both"  {
	
			//du/dthr1
			replace `d1a' = `pu1' * normalden(`ztop11') /`su1' + `pu2' * normalden(`ztop12') /`su2' if `tmpv' == 1
			replace `d1a' = 0 if `tmpv' == 0
			replace `d1b' = `pu1' * normalden(`zbot11') /`su1' + `pu2' * normalden(`zbot12') /`su2' if `tmpv1' == 1
			replace `d1b' = 0  if `tmpv1' == 0
		}
		
		replace `gg`i'' = `tmp11' *`d1a'  - `tmp12' *`d1b'  - `tmp13' *`d1a' + `tmp14' *`d1b'
	}
	
	
	local n1 = `nn'+ $bc_Nthr1 + 1
	local n2 = `nn'+ $bc_Nthr1 + $bc_Nthr2
		
	forvalues i = `n1'/`n2' {
			
		replace `d1a' = .
		replace `d1b' = .
		replace `tmpv' = 0
		replace `tmpv' = 1 if (`ttop2' == `thr2'[`i'-$bc_Nthr1 -`nn'+1 ,1])
			
				
		replace `tmpv1' = 0
		replace `tmpv1' = 1 if (`bbot2' == `thr2'[`i'-$bc_Nthr1 -`nn'+1 ,1])
				
				
		if "$bc_mix"=="none" | "$bc_mix"=="mix1"   {
			//du/dthr2
			replace `d1a' =  normalden(`ttop2'-`xb2') if `tmpv' == 1
			replace `d1a' = 0 if `tmpv' == 0
			replace `d1b' = normalden(`bbot2'-`xb2')  if `tmpv1' == 1
			replace `d1b' = 0  if `tmpv1' == 0
		}
	
		if "$bc_mix"=="mix2" | "$bc_mix"=="equal"  | "$bc_mix"=="both"  {
	
			//du/dthr2
			replace `d1a' = `pv1' * normalden(`ztop21') /`sv1' + `pv2' * normalden(`ztop22') /`sv2' if `tmpv' == 1
			replace `d1a' = 0 if `tmpv' == 0
			replace `d1b' = `pv1' * normalden(`zbot21') /`sv1' + `pv2' * normalden(`zbot22') /`sv2' if `tmpv1' == 1
			replace `d1b' = 0  if `tmpv1' == 0
		}
			
		replace `gg`i'' = `tmp21' *`d1a'  - `tmp22' *`d1a'  - `tmp23' *`d1b' + `tmp24' *`d1b'
			
	}
	
	 //dependency parameter
	local ndep = `nn'+$bc_Nthr1 + $bc_Nthr2 +1
	if `theta'==0 & (("$bc_copu"=="frank") | ("$bc_copu"=="clayton")) {	
		replace `gg`ndep'' = 0
	}
	
	else if ("$bc_copu"!="indep"){		
		forvalues i = 31/34 { 
			tempvar tmp`i' 
			gen double `tmp`i'' = . 
		}
			
		ddcopu `top1' `top2' `theta' `tmp31'
		ddcopu `bot1' `top2' `theta' `tmp32'
		ddcopu `top1' `bot2' `theta' `tmp33'
		ddcopu `bot1' `bot2' `theta' `tmp34'
			
					
		if ("$bc_copu"=="gaussian") {
		
			replace `gg`ndep'' = (`tmp31' - `tmp32' - `tmp33' +`tmp34') * (1-`theta'^2)
		}
		if ("$bc_copu"=="frank") {
		
			replace `gg`ndep'' = `tmp31' - `tmp32' - `tmp33' +`tmp34'
		}
		if  ("$bc_copu"=="clayton") {
			replace `gg`ndep'' = (`tmp31' - `tmp32' - `tmp33' +`tmp34' ) * (`theta') //exp(/depend) = theta
		}
			
		if ("$bc_copu"=="gumbel") | ("$bc_copu"=="joe") {
			replace `gg`ndep'' = (`tmp31' - `tmp32' - `tmp33' +`tmp34' ) * (`theta' - 1)  //exp(/depend) = theta -1
		}
			
	}
	
	// mixing parameters

	
	local n2 = `nn'+$bc_Nthr1 + $bc_Nthr2 +2 - `adjind'
	
	
	
			
	if "$bc_mix"=="mix1" | "$bc_mix"=="both" | "$bc_mix"=="equal" {
		//du/db1 (b1 - first mixing parameter)
		replace `d1a' = .
		replace `d1b' = .
		tempname der1 der2 der3 der4
		scalar `der1' = `pu1'* (1- `pu1')  // dpu1/db1
		scalar `der2' = `mu2' * `pu2' /`pu1'    //dmu1/db1
		scalar `der4' = - `der1'    //dpu2/db1
		scalar `der3' = (-1/(2*`su1' * `pu1'^2))*((`su2'^2 + `mu2'^2) * `der4' * `pu1' + (1- `pu2' * (`su2'^2 + `mu2'^2)) * `der1'  + 2* `mu1' * `der2' * `pu1'^2)   //dsu1/db1
			
				
		replace `d1a' = `der1' * normal(`ztop11') - (`pu1'/`su1') * normalden(`ztop11') * (`der2' +`ztop11' * `der3') + `der4' * normal(`ztop12')
		replace `d1b' = `der1' * normal(`zbot11') - (`pu1'/`su1') * normalden(`zbot11') * (`der2' +`zbot11' * `der3') + `der4' * normal(`zbot12')
			
				
		replace `gg`n2'' = `tmp11' *`d1a'  - `tmp12' *`d1b'  - `tmp13' *`d1a' + `tmp14' *`d1b'
				
		//du/db2 (b2 - second mixing parameter)
		local i1 = `nn'+$bc_Nthr1 + $bc_Nthr2 +3 - `adjind'
		tempname der5 der6 der7
		replace `d1a' = .
		replace `d1b' = .
				
		scalar `der5' = -`pu2'/`pu1'  // dmu1/db2
		scalar `der6' = (`pu2' * (`mu1'-`mu2'))/(`su1' * `pu1') // dsu1/db2
		scalar `der7' = 1 // dmu2/db2
				
		replace `d1a' = -`pu1' * normalden(`ztop11') *( `der5' +`ztop11' * `der6' ) / `su1' - `pu2' *  normalden(`ztop12') * `der7' / `su2'
		replace `d1b' = -`pu1' * normalden(`zbot11') *( `der5' +`zbot11' * `der6' ) / `su1' - `pu2' *  normalden(`zbot12') * `der7' / `su2'
			
		replace `gg`i1'' = `tmp11' *`d1a'  - `tmp12' *`d1b'  - `tmp13' *`d1a' + `tmp14' *`d1b'
				
				
	
		//du/db3 (b3 - third mixing parameter)
		
		local i2 = `nn'+$bc_Nthr1 + $bc_Nthr2 +4 - `adjind'
		tempname der8
		replace `d1a' = .
		replace `d1b' = .
			
		scalar `der8' = - (`pu2' * `su2') / (`pu1' * `su1')  // dsu1/db3
		replace `d1a' = -`pu1' * normalden(`ztop11') * `ztop11' * `der8' / `su1' - `pu2' *  normalden(`ztop12') * `ztop12' / `su2'
		replace `d1b' = -`pu1' * normalden(`zbot11') * `zbot11' * `der8' / `su1' - `pu2' *  normalden(`zbot12') * `zbot12' / `su2'
		replace `gg`i2'' =sign(`b31') * ( `tmp11' *`d1a'  - `tmp12' *`d1b'  ///
		   - `tmp13' *`d1a' + `tmp14' *`d1b') // this is du/dsu2 * dsu2/db3
			                                      // and dsu2/db3=b3/su2 = sign(b3)                                                      
	}
			
	if "$bc_mix"=="mix2" | "$bc_mix"=="both" | "$bc_mix"=="equal" {
		//dv/db1 (b1 - first mixing parameter)
		replace `d1a' = .
		replace `d1b' = .
		tempname der1 der2 der3 der4
		scalar `der1' = `pv1' - `pv1'^2  // dpv1/db1
		scalar `der2' = `mv2' * `pv2' / `pv1'  //dmv1/db1
		scalar `der4' = - `der1'    //dpv2/db1
		scalar `der3' = (-1/(2*`sv1' * `pv1'^2))*((`sv2'^2 + `mv2'^2) * `der4' * `pv1' + (1- `pv2' * (`sv2'^2 + `mv2'^2)) * `der1'  + 2* `mv1' * `der2' * `pv1'^2)    //dsv1/db1
				
				
		replace `d1a' = `der1' * normal(`ztop21') - (`pv1'/`sv1') * normalden(`ztop21') * (`der2' +`ztop21' * `der3') + `der4' * normal(`ztop22')
		replace `d1b' = `der1' * normal(`zbot21') - (`pv1'/`sv1') * normalden(`zbot21') * (`der2' +`zbot21' * `der3') + `der4' * normal(`zbot22')
			
		tempvar tempg1 tempg2 tempg3
		gen double `tempg1' = `tmp21' *`d1a'  - `tmp22' *`d1a'  - `tmp23' *`d1b' + `tmp24' *`d1b'
				
		//dv/db2 (b2 - second mixing parameter)
		local i1 = `nn'+$bc_Nthr1 + $bc_Nthr2 +3 - `adjind'
		tempname der5 der6 der7
		replace `d1a' = .
		replace `d1b' = .
				
		scalar `der5' = -`pv2'/`pv1'  // dmv1/db2
		scalar `der6' = (`pv2' * (`mv1'-`mv2'))/(`sv1' * `pv1') // dsv1/db2
		scalar `der7' = 1 // dmv2/db2
				
		replace `d1a' = -`pv1' * normalden(`ztop21') *( `der5' +`ztop21' * `der6' ) / `sv1' - `pv2' *  normalden(`ztop22') * `der7' / `sv2'
		replace `d1b' = -`pv1' * normalden(`zbot21') *( `der5' +`zbot21' * `der6' ) / `sv1' - `pv2' *  normalden(`zbot22') * `der7' / `sv2'
				
		gen double `tempg2' = `tmp21' *`d1a'  - `tmp22' *`d1a'  - `tmp23' *`d1b' + `tmp24' *`d1b'
				
				
		//dv/db3 (b3 - third mixing parameter)
			
		local i2 = `nn'+$bc_Nthr1 + $bc_Nthr2 +4 - `adjind'
		tempname der8
		replace `d1a' = .
		replace `d1b' = .
			
		scalar `der8' = - (`pv2' * `sv2') / (`pv1' * `sv1')  // dsv1/db3
		replace `d1a' = -`pv1' * normalden(`ztop21') * `ztop21' * `der8' / `sv1' - `pv2' *  normalden(`ztop22') * `ztop22' / `sv2'
		replace `d1b' = -`pv1' * normalden(`zbot21') * `zbot21' * `der8' / `sv1' - `pv2' *  normalden(`zbot22') * `zbot22' / `sv2'

		gen double `tempg3' =sign(`b32') * ( `tmp21' *`d1a'  - `tmp22' *`d1a'  ///
		   - `tmp23' *`d1b' + `tmp24' *`d1b') // this is dv/dsv2 * dsv2/db3
				                                      // and dsv2/db3=b3/sv2 = sign(b3)   
		if "$bc_mix"=="mix2" {
			replace `gg`n2'' = `tempg1'
			replace `gg`i1'' = `tempg2'
			replace `gg`i2'' = `tempg3'
		}
		else if  "$bc_mix"=="both"  {
			local j1 = `nn'+$bc_Nthr1 + $bc_Nthr2 +5 - `adjind'
			local j2 = `nn'+$bc_Nthr1 + $bc_Nthr2 +6 - `adjind'
			local j3 = `nn'+$bc_Nthr1 + $bc_Nthr2 +7 - `adjind'
			replace `gg`j1'' = `tempg1'
			replace `gg`j2'' = `tempg2'
			replace `gg`j3'' = `tempg3'	
		}
			
		else if "$bc_mix"=="equal" {
			replace `gg`n2'' = `gg`n2'' + `tempg1'
			replace `gg`i1'' = `gg`i1'' + `tempg2'
			replace `gg`i2'' = `gg`i2'' + `tempg3'
		}
			
	}
	
	forvalues i = 1/`n_grad' {
		replace `g`i''=`gg`i'' / `lik'
	}

}
end



program define copu
version 13
args u v theta cop
quietly {
	if "$bc_copu"=="indep" {
		 replace `cop'=`u'*`v'
	}
	if "$bc_copu"=="gaussian" {
		replace `cop'=binorm(invnormal(`u'),invnormal(`v'),`theta') ///
										if `u'>0&`u'<1&`v'>0&`v'<1
	}
	if "$bc_copu"=="frank" {
		if `theta'==0 replace `cop'=`u'*`v'
		else {
			tempname eta
			scalar `eta'=1-exp(-`theta')
			replace `cop'=`eta'-(1-exp(-`theta'*`u'))*(1-exp(-`theta'*`v')) ///
										if `u'>0&`u'<1&`v'>0&`v'<1
			replace `cop'=-(ln(`cop'/`eta'))/`theta' if `u'>0&`u'<1&`v'>0&`v'<1
		}
	}
	if "$bc_copu"=="clayton" {
		if `theta'==0 replace `cop'=`u'*`v'
		else {
			replace `cop'=(`u'^(-`theta')) + (`v'^(-`theta')) - 1 /// 
										if `u'>0&`u'<1&`v'>0&`v'<1
			replace `cop'=(max(0,`cop'))^(-1/`theta') if `u'>0&`u'<1&`v'>0&`v'<1
		}
	}
	if "$bc_copu"=="gumbel" {
		replace `cop'=exp(  - ((-ln(`u'))^`theta'+(-ln(`v'))^`theta') /// 
						^(1/`theta')  ) if `u'>0&`u'<1&`v'>0&`v'<1
	}
	if "$bc_copu"=="joe" {
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



program define dcopu
/*derivative of likelihood parts wrt u (k=1) and v (k=2)*/
version 13
args u v theta k dcop
tempvar eta
quietly {
	if "$bc_copu"=="indep" {
		if `k' == 1 {
			replace `dcop' = `v'
		}
		else {
			replace `dcop' = `u'
		}
	
	}
	if "$bc_copu"=="gaussian" {
		
		if `k' == 1 {
			replace `dcop' = normal((invnormal(`v')- `theta' * invnormal(`u'))/ sqrt(1-`theta'^2))
			replace `dcop' = 0 if `u' == 1 &`v'>0&`v'<1
			replace `dcop' = 1 if `v' == 1 &`u'>0&`u'<1
		}
		else {
			replace `dcop' = normal((invnormal(`u')- `theta' * invnormal(`v'))/ sqrt(1-`theta'^2))
			replace `dcop' = 0 if `v' == 1 &`u'>0&`u'<1
			replace `dcop' = 1 if `u' == 1 &`v'>0&`v'<1
		}
		
	}
	if "$bc_copu"=="frank" {
		if `theta'==0 {
			if `k' == 1 {
				replace `dcop'=`v'
			}
			else {
				replace `dcop'=`u'
			}
		}
		else {
			
			tempname eta
			scalar `eta'=-1+exp(-`theta')
			replace `dcop'=1+ (((-1+exp(-`theta'*`u'))*(-1+exp(-`theta'*`v')))/`eta')										
			replace `dcop'=`eta'*`dcop' 
			if `k' == 1 {
				replace `dcop' = ((exp(-`theta'*`u'))*(-1+exp(-`theta'*`v'))) / `dcop'
				replace `dcop' = 0 if `u' == 1 &`v'>0&`v'<1
				replace `dcop' = 1 if `v' == 1 &`u'>0&`u'<1
			}
			else {
				replace `dcop' = ((exp(-`theta'*`v'))*(-1+exp(-`theta'*`u'))) / `dcop'
				replace `dcop' = 0 if `v' == 1 &`u'>0&`u'<1
				replace `dcop' = 1 if `u' == 1 &`v'>0&`v'<1
			}
			
		}
	}
	if "$bc_copu"=="clayton" {
		if `theta'==0 {
			if `k' == 1 {
				replace `dcop'=`v'
			}
			else {
				replace `dcop'=`u'
			}
		}
		else {
			replace `dcop'=(`u'^(-`theta')) + (`v'^(-`theta')) - 1 
			replace `dcop' = (`dcop')^(-(`theta'+1)/`theta')
			if `k' == 1 {
				replace `dcop'=`dcop' * (`u')^(-`theta'-1)
				replace `dcop' = 0 if `u' == 1 &`v'>0&`v'<1
				replace `dcop' = 1 if `v' == 1 &`u'>0&`u'<1
			}
			else {
				replace `dcop'=`dcop' * (`v')^(-`theta'-1)
				replace `dcop' = 0 if `v' == 1 &`u'>0&`u'<1
				replace `dcop' = 1 if `u' == 1 &`v'>0&`v'<1
			}
		}
	}
	if "$bc_copu"=="gumbel" {
		tempvar uu vv
		gen double `uu' = (-ln(`u'))^`theta'
		gen double `vv' = (-ln(`v'))^`theta'
		
		if `k' == 1 {
			replace `dcop'= - ((exp(-(`uu'+`vv')^(1/`theta'))) * ((`uu'+`vv')^(1/`theta')) * `uu')/ (`u' * (`uu' + `vv') * ln(`u'))
			replace `dcop' = 0 if `u' == 1 &`v'>0&`v'<1
			replace `dcop' = 1 if `v' == 1 &`u'>0&`u'<1
		}
		else {
			replace `dcop'= - ((exp(-(`uu'+`vv')^(1/`theta'))) * ((`uu'+`vv')^(1/`theta')) * `vv')/ (`v' * (`uu' + `vv') * ln(`v'))
			replace `dcop' = 0 if `v' == 1 &`u'>0&`u'<1
			replace `dcop' = 1 if `u' == 1 &`v'>0&`v'<1
		}
	}
	if "$bc_copu"=="joe" {
		replace `dcop' =-( (((1-`u')^`theta') + ((1-`v')^`theta') ///
					-((1-`u')^`theta')*((1-`v')^`theta') )^((1-`theta')/`theta'))
		if `k' == 1 {
			replace `dcop'=`dcop' * ((1-`u')^(`theta'-1)) * (-1+(1-`v')^(`theta'))
			replace `dcop' = 0 if `u' == 1 &`v'>0&`v'<1
			replace `dcop' = 1 if `v' == 1 &`u'>0&`u'<1
		}
		else {
			replace `dcop'=`dcop' * ((1-`v')^(`theta'-1)) * (-1+(1-`u')^(`theta'))
			replace `dcop' = 0 if `v' == 1 &`u'>0&`u'<1
			replace `dcop' = 1 if `u' == 1 &`v'>0&`v'<1
		}			
	}
	replace `dcop' = 0 if `u' == 1 & `v' == 1
	replace `dcop' = 0 if `u' == 0 | `v' == 0 
}
end

program define ddcopu
/*derivative of likelihood parts wrt theta*/
version 13
args u v theta dcop

quietly {
	if "$bc_copu"=="gaussian" {
		tempvar uu vv
		gen double `uu' =invnormal(`u')
		gen double `vv' =invnormal(`v')
		replace `dcop' = (exp(-(`uu'^2 + `vv'^2 - 2*`theta'*`uu' * `vv')/(2*(1-`theta'^2))))/(2 * c(pi) *sqrt(1-`theta'^2))
	}
	if "$bc_copu"=="frank" {
		if `theta'==0 {
			replace `dcop' = 0 //the derivative wrt theta is zero as theta doesnt exist. 
		}
		else {
			tempname eta
			scalar `eta'=-1+exp(-`theta')
			tempvar uu vv insi
			gen double `uu' = -1 + exp(-`theta'*`u')
			gen double `vv' = -1 + exp(-`theta'*`v')
			gen double `insi' = 1+ `uu' * `vv'/`eta'
			
			replace `dcop'=(ln(`insi'))/ (`theta'^2)
			tempvar a b c
			gen double `a' = - (`u' * (exp(-`theta'*`u')) * `vv') / `eta'
			gen double `b' = - (`v' * (exp(-`theta'*`v')) * `uu') / `eta'
			gen double `c' = (`uu' * `vv' * exp(-`theta'))/ (`eta'^2)
			
			
			replace `dcop' = `dcop' - ((`a')+(`b')+(`c')) / (`theta' * `insi')
			
		}
	}
	if "$bc_copu"=="clayton" {
		if `theta'==0 {
			replace `dcop' = 0
		}
		else {
			tempvar f 
			gen double `f' = (`u'^(-`theta')) + (`v'^(-`theta')) -1
			replace `dcop' = (`f'^(-1/`theta')) * (((ln(`f'))/(`theta'^2))+(((`u'^(-`theta'))* ln(`u') + (`v'^(-`theta'))*ln(`v') )/(`theta' * `f')))
			
		}
	
	}
	if "$bc_copu"=="gumbel" {
		tempvar uu vv
		gen double `uu' = (-ln(`u'))^`theta'
		gen double `vv' = (-ln(`v'))^`theta'
		
		replace `dcop' = -(exp(-(`uu'+`vv')^(1/`theta'))) * ((`uu'+`vv')^(1/`theta')) * ///
		     ( - ((ln(`uu'+`vv'))/ (`theta'^2)) + (( ((-ln(`u'))^`theta') * ln(-ln(`u')) + ///
			 ((-ln(`v'))^`theta') * ln(-ln(`v'))  )/ (`theta' * (`uu'+`vv'))))
		
	}
	if "$bc_copu"=="joe" {
		
				
		tempvar uu vv a
		gen double `uu' = ( 1-`u' )^( `theta' )
		gen double `vv' = ( 1-`v' )^( `theta' )
		gen double `a' = `uu' + `vv' - `uu' * `vv' 
		replace `dcop' = - (( `a' )^(1/ `theta' ) ) * ///
		  (- (ln( `a' ) )/ ( `theta' ^2) + ((`uu' * ln(1-`u')+ `vv' * ln(1-`v') ///
		  -`uu' * `vv' * ln(1-`u')-`uu' * `vv' * ln(1-`v'))/(`theta' * `a')) )
		  
		replace `dcop' = -(( `vv' )^(1/`theta')) * (- (ln(`vv')/(`theta'^2))+( (ln(1-`v') )/`theta' ))   if `u' ==1 & `v' != 1
		replace `dcop' = -(( `uu' )^(1/`theta')) * (- (ln(`uu')/(`theta'^2))+( (ln(1-`u') )/`theta' ))   if `v' ==1 & `u' != 1

	}
	replace `dcop' = 0 if `u' == 1 | `v' == 1 
	replace `dcop' = 0 if `u' == 0 | `v' == 0 
}
end

