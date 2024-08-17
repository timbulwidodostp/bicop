{smcl}
{* documented: 16 December 2015}{...}
{hi:help bicop}{right: ({browse "http://www.stata-journal.com/article.html?article=st0429":SJ16-1: st0429})}
{right:also see:  {help bicop postestimation}}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col: {cmd:bicop} {hline 2}}Generalized bivariate ordinal regression
model{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{phang}Syntax 1: Same independent variables in both equations

{p 8 17 2}
{cmd:bicop} 
{it:{help depvar:depvar1}} {it:{help depvar:depvar2}}
[{indepvars}] {ifin} [{it:{help bicop##weight:weight}}]
[{cmd:,}
{it:{help bicop##syntax1_options:syntax1_options}}]


{phang}Syntax 2: Different independent variables in each equation

{p 8 17 2}
{cmd:bicop}
{cmd:(}{it:equation1}{cmd:)} {cmd:(}{it:equation2}{cmd:)}
{ifin}
[{it:{help bicop##weight:weight}}]
[{cmd:,} {it:{help bicop##syntax2_options:syntax2_options}}]


{phang}
{it:equation1} and {it:equation2} are specified as 

{p 12 12 2}{cmd:(}[{it:eqname}{cmd::}] {depvar} [{cmd:=}] [{indepvars}]
		[{cmd:,} {opth off:set(varname)}]{cmd:)}


{synoptset 20}{...}
{marker syntax1_options}{...}
{synopthdr:syntax1_options}
{synoptline}
{synopt :{opt mix:ture(mixturetype)}}specify marginal distribution of each
residual; {it:mixturetype} may be {opt none}, {opt mix1}, {opt mix2}, or
{opt both}; default is {cmd:mixture(none)}{p_end}
{synopt :{opt cop:ula(copulatype)}}specify copula function to control the
pattern of stochastic dependence; {it:copulatype} may be {opt indep},
{opt gaussian}, {opt frank}, {opt clayton}, {opt gumbel}, or {opt joe};
default is {cmd:copula(gaussian)}{p_end}
{synopt :{opth c:onstraints(numlist)}}apply specified linear constraints{p_end}
{synopt:{opt col:linear}}keep collinear variables{p_end}
{synopt :{opth offset1(varname)}}specify offset variable for first equation{p_end}
{synopt :{opth offset2(varname)}}specify offset variable for second equation{p_end}
{synopt :{opth vce(vcetype)}}specify how to estimate variance-covariance
matrix; {it:vcetype} may be {opt oim}, {opt r:obust}, {opt cl:uster}
{it:clustvar}, or {opt opg}{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt from(init_specs)}}specify initial values for the coefficients{p_end}
{synopt :{opt search(spec)}}specify whether to use {helpb ml}'s initial search
algorithm{p_end}
{synopt :{opt repeat(#)}}specify the number of random attempts to find a
better initial-value vector; must specify with {cmd:search()}{p_end}
{synopt :{it:{help bicop##bicop_maximize:maximize_options}}}control the
maximization process; some options may be especially useful{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 20}{...}
{marker syntax2_options}{...}
{synopthdr :syntax2_options}
{synoptline}
{synopt :{opt mix:ture(mixturetype)}}specify marginal distribution of each
residual; {it:mixturetype} may be {opt none}, {opt mix1}, {opt mix2}, or
{opt both}; default is {cmd:mixture(none)}{p_end}
{synopt :{opt cop:ula(copulatype)}}specify copula function to control the
pattern of stochastic dependence; {it:copulatype} may be {opt indep},
{opt gaussian}, {opt frank}, {opt clayton}, {opt gumbel}, or {opt joe};
default is {cmd:copula(gaussian)}{p_end}
{synopt :{opth c:onstraints(numlist)}}apply specified linear constraints{p_end}
{synopt:{opt col:linear}}keep collinear variables{p_end}
{synopt :{opth vce(vcetype)}}specify how to estimate variance-covariance
matrix; {it:vcetype} may be {opt oim}, {opt r:obust}, {opt cl:uster}
{it:clustvar}, or {opt opg}{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt from(init_specs)}}specify initial values for the coefficients{p_end}
{synopt :{opt search(spec)}}specify whether to use {helpb ml}'s initial search
algorithm{p_end}
{synopt :{opt repeat(#)}}specify the number of random attempts to find a
better initial-value vector; must specify with {cmd:search()}{p_end}
{synopt :{it:{help bicop##bicop_maximize:maximize_options}}}control the
maximization process; some options may be especially useful{p_end}
{synoptline}
{p2colreset}{...}

INCLUDE help fvvarlist
{p 4 6 2}{it:depvar1}, {it:depvar2}, {it:indepvars}, and {it:depvar} may
contain time-series operators; see {help tsvarlist}.{p_end}
{marker weight}{...}
{p 4 6 2}{opt pweight}s, {opt fweight}s, and {opt iweight}s are allowed; see 
{help weight}.{p_end}
{p 4 6 2} See {help bicop postestimation} for features available after
estimation.{p_end}
{p 4 6 2}{cmd:bicop} typed without arguments redisplays previous results.
{p_end}


{title:Description}

{pstd}
{cmd:bicop} is a user-written command that fits a generalized bivariate
ordinal regression model using maximum likelihood estimation.  It is
implemented as an {cmd:lf1 ml} evaluator.  The model involves a pair of latent
regression equations, each with a standard threshold-crossing condition to
generate ordinal observed dependent variables.  The bivariate residual
distribution is specified to have marginals, each with the form of a two-part
normal mixture, and a choice of copula functions to represent the pattern of
dependence between the two residuals.


{title:Options}

{phang}
{opt mixture(mixturetype)} specifies the marginal distribution of each
residual.  There are five choices for {it:mixturetype}: {opt none} specifies
that each marginal distribution be N(0,1); {opt mix1} specifies that the
residual from equation 1 has a two-part normal mixture distribution but that
the residual from equation 2 be N(0,1); {opt mix2} specifies N(0,1) for
equation 1 and a normal mixture for equation 2; {opt both} allows each
residual to have a different normal mixture distribution; and {opt equal}
specifies that both residuals have the same normal mixture distribution. The
default is {cmd:mixture(none)}.

{phang}
{opt copula(copulatype)} specifies the copula function to be
used to control the pattern of stochastic dependence of the two residuals.
There are six choices for {it:copulatype}: {opt indep}, which specifies the
special form c(u,v)=uv, {opt gaussian}, {opt clayton}, {opt frank},
{opt gumbel}, and {opt joe}.  The default is {cmd:copula(gaussian)}.  Note
that if both {cmd:mixture()} and {cmd:copula()} are omitted, the {cmd:bicop}
command produces the same results as the existing {cmd:bioprobit} (Sajaia
2008) and (if both dependent variables are binary) {cmd:biprobit} commands.

{phang}
{opt c:onstraints(numlist)} applies specified linear constraints; see
{manhelp constraints R}.

{phang}
{opt collinear} retains collinear variables.  Usually, there is no
reason to leave collinear variables in place, and doing so would cause the
estimation to fail because of matrix singularity.  However, in some
constrained cases, the model may be fully identified despite the collinearity.
The {opt collinear} option then allows estimation to occur, leaving the
equations with collinear variables intact.  This option is seldom used.

{phang}
{opth offset1(varname)} and {opth offset2(varname)} specify an offset
variable for the first equation and the second equation, respectively.

{phang}
{opt vce(vcetype)} specifies how to estimate the
variance-covariance matrix corresponding to the parameter estimates.
The supported options are {opt oim}, {opt opg}, {opt robust}, and
{opt cluster}.  The current version of the command does not allow
{opt bootstrap} or {opt jackknife} estimators.  See
{manhelpi vce_option R}.

{phang}
{opt level(#)} sets the significance level to be used for
confidence intervals; see {manhelp level R}.

{marker bicop_maximize}{...}
{phang}
{opt from(init_specs)}, where {it:init_specs} is either
{it:matname}, the name of a matrix containing the starting values, or
{it:matname}{cmd:,} {cmd:copy}|{cmd:skip}.  The {cmd:copy}
suboption specifies that the initialization vector be copied into the
initial-value vector by position rather than by name, and the {cmd:skip}
suboption specifies that any irrelevant parameters found in the specified
initialization vector be ignored.  Poor values in {cmd:from()} may lead to
convergence problems.

{phang}
{opt search(spec)} specifies whether {opt ml}'s ({manhelp ml R})
initial search algorithm is used.  {it:spec} may be {opt on} or {opt off}.

{phang}
{opt repeat(#)} specifies the number of random attempts to be made to find a
better initial-value vector.  This option should be used in conjunction with
{opt search()}.

{phang}
{it:maximize_options}: {opt dif:ficult}, {opt tech:nique(algorithm_spec)},
{opt iter:ate(#)}, [{cmdab:no:}]{opt lo:g}, {opt tr:ace},
{opt grad:ient}, {opt showstep}, {opt hess:ian},
{opt show:tolerance}, {opt tol:erance(#)}, {opt ltol:erance(#)},
{opt gtol:erance(#)}, {opt nrtol:erance(#)}, and {opt nonrtol:erance}; see
{manhelp maximize R}.


{title:Remarks}

{pstd}
The likelihood functions of mixture models are known to have multiple optima
(McLachlan and Peel 2000).  We recommend that several different
starting values be used to ensure convergence to the global maximum.

{pstd}
Copulas can represent different types of dependence.  Nonconvergence of
a particular copula might be due to discrepancies between the pattern of
dependence in the data and implied by the copula (Trivedi and Zimmer 2005).


{title:Stored results}

{pstd}
{cmd:bicop} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_aux)}}number of auxiliary parameters{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(chi2_ind)}}chi-squared of independence test{p_end}
{synopt:{cmd:e(df_ind)}}degrees of freedom of independence test{p_end}
{synopt:{cmd:e(p_ind)}}significance of independence test{p_end}

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:bicop}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}names of dependent variable{p_end}
{synopt:{cmd:e(indepvars1)}}names of independent variables in equation 1{p_end}
{synopt:{cmd:e(indepvars2)}}names of independent variables in equation 2{p_end}
{synopt:{cmd:e(copula)}}name of copula{p_end}
{synopt:{cmd:e(mixture)}}name of mixture{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(offset1)}}offset for first equation{p_end}
{synopt:{cmd:e(offset2)}}offset for second equation{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared
	test{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform
                         maximization or minimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(Cns)}}constraints matrix{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(extpar)}}matrix of the transformed parameters corresponding to the dependency and mixtures parameters{p_end}

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. use ukhlsfwb}{p_end}

{pstd}Generalized bivariate ordinal regression model with N(0,1) marginals
using a Clayton copula{p_end}
{phang2}{cmd:. bicop finnow finfut female homeowner unempsick, copula(clayton)}{p_end}

{pstd}Generalized bivariate ordinal regression model with N(0,1) marginals
using a Clayton copula and different variables in each equation {p_end}
{phang2}{cmd:. bicop (finnow =  female homeowner unempsick) (finfut =  homeowner unempsick), copula(clayton) }{p_end}

{pstd}Generalized bivariate ordinal regression model using a Clayton copula
and a two-component normal mixture distribution for U and a pure normal
distribution for V; user-supplied starting parameter values are also
provided{p_end}
{phang2}{cmd:. matrix x =(-.1651802, .5308612, -.7150698, -.0501328, -.2089726, -.1409637, -1.644831, -.9164829, .1108187, 1.05824, -1.054552, .4758982, 0773722, -2.246348, .1525841, .902352)}{p_end}
{phang2}{cmd:. bicop finnow finfut female homeowner unempsick, copula(clayton) mixture(mix1) from(x,copy) search(off)}{p_end}


{title:References}

{phang}
McLachlan, G., and D. Peel. 2000. {it:Finite Mixture Models}.
New York: Wiley.

{phang}
Sajaia, Z. 2008. bioprobit: Stata module for bivariate ordered probit
regression. Statistical Software Components S456920, Department of
Economics, Boston College.
{browse "https://ideas.repec.org/c/boc/bocode/s456920.html"}.

{phang}
Trivedi, P. K., and D. M. Zimmer. 2005. Copula modeling: An introduction for
practitioners. {it:Foundations and Trends in Econometrics} 1: 1-111.


{title:Authors}

{phang}Monica Hernandez-Alava, HEDS, ScHARR, University of Sheffield, UK{p_end}
{phang}monica.hernandez@sheffield.ac.uk{p_end}
{phang}{browse "http://www.shef.ac.uk/scharr/sections/heds/staff/hernandez_m"}

{phang}Stephen Pudney, ISER, University of Essex, UK{p_end}
{phang}spudney@essex.ac.uk{p_end}
{phang}{browse "https://www.iser.essex.ac.uk/people/spudney"}


{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 16, number 1: {browse "http://www.stata-journal.com/article.html?article=st0429":st0429}

{p 7 14 2}
Help:  {help bicop postestimation}
{p_end}
