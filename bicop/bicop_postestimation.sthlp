{smcl}
{* documented: 13 December 2015}{...}
{hi:help bicop postestimation}{right: ({browse "http://www.stata-journal.com/article.html?article=st0429":SJ16-1: st0429})}
{right:also see:  {helpb bicop}}
{hline}

{title:Title}

{p2colset 5 29 31 2}{...}
{p2col: {hi:bicop postestimation} {hline 2}}Postestimation tools for bicop
{p2colreset}{...}


{title:Postestimation commands}

{pstd}
The following postestimation commands are available after {cmd:bicop}:

{synoptset 14 tabbed}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb estat}}postestimation statistics{p_end}
INCLUDE help post_estimates
INCLUDE help post_lincom
INCLUDE help post_lrtest_star
INCLUDE help post_nlcom
{synopt :{helpb biprobit postestimation##predict:predict}}predictions, residuals, influence statistics, and other diagnostic measures{p_end}
INCLUDE help post_test
INCLUDE help post_testnl
{synoptline}
{p2colreset}{...}
{phang}* {cmd:lrtest} is not appropriate with {cmd:svy} estimation results.
{p_end}


{marker syntax_predict}{...}
{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} 
{dtype}
{newvar} 
{ifin}
[{cmd:,} {it:statistic} {opt o:utcome(outcome_pair)} {opt nooff:set}]

{p 8 16 2}
{cmd:predict}
{dtype}
{c -(}{it:stub}{cmd:*}{c |}{it:{help newvar:newvar_eq1}}
                     {it:{help newvar:newvar_eq2}} ...
                     {it:{help newvar:newvar_eqk}}{c )-}
{ifin}
{cmd:,}
{opt sc:ores}

{synoptset 14}{...}
{synopthdr :statistic}
{synoptline}
{synopt :{opt pr}}predicted probabilities; the default{p_end}
{synopt :{opt xb1}}linear prediction for equation 1{p_end}
{synopt :{opt xb2}}linear prediction for equation 2{p_end}
{synopt :{opt stdp1}}standard error of the linear prediction for equation 1{p_end}
{synopt :{opt stdp2}}standard error of the linear prediction for equation 2{p_end}
{synopt :{opt pcond1}}conditional (on equation 2) probability; {opt o:utcome(outcome_pair)} is required {p_end}
{synopt :{opt pcond2}}conditional (on equation 1) probability; {opt o:utcome(outcome_pair)} is required {p_end}
{synoptline}
{p2colreset}{...}
INCLUDE help esample


{title:Description for predict}

{pstd}
Following {cmd:bicop}, the {cmd:predict} command can be used to
construct several alternative predictions.  The predictions include the linear
indices X_{i1}B_1 and X_{i2}B_2 and corresponding standard errors;
probabilities of the form Pr(Y_{ij}=r|X_{ij}) or
Pr(Y_{i1}=r,Y_{i2}=s|X_{i1},X_{i2}); and conditional probabilities of the
form Pr(Y_{ij}=r|Y_{ik}=s,X_{i1},X_{i2}).


{marker options_predict}{...}
{title:Options for predict}

{phang}
{opt pr} calculates the predicted probabilities.

{phang}
{opt xb1} calculates the linear prediction for equation 1.

{phang}
{opt xb2} calculates the linear prediction for equation 2.

{phang}
{opt stdp1} calculates the standard error of the linear prediction of
equation 1.

{phang}
{opt stdp2} calculates the standard error of the linear prediction of
equation 2.

{phang}
{opt pcond1} calculates the conditional probability Pr(depvar1 = A | depvar2 =
B) using {opt outcome(A,B)}.

{phang}
{opt pcond2} calculates the conditional probability Pr(depvar2 = A | depvar1 =
B) using {opt outcome(B,A)}.

{phang}
{opt outcome(outcome_pair)} specifies the outcome pair for which the predicted
probabilities are calculated.  {opt outcome()} should contain a pair of
either values of the dependent variables or one of {opt #1}, {opt #2},
..., with {opt #1} meaning the first category of a dependent variable,
{opt #2} the second category, etc.  If {opt outcome(outcome_pair)} is not
specified, {cmd:predict} will compute {opt outcome(#1,#1)}.  If one of the
arguments is missing, then {cmd:predict} will return the marginal probability;
that is, {opt outcome(., k)} will return Pr(y2=k).

{phang}
{opt nooffset} is relevant only if you specified {opth offset1(varname)} or 
{opt offset2(varname)} for
{cmd:bicop}.  It modifies the calculations made by {opt predict} so that
they ignore the offset variables; the linear predictions are treated as 
xb1 rather than as xb1 + offset1 and as xb2 rather than as xb2 + offset2.

{phang}
{opt scores} calculates equation-level score variables.

{pmore}
The first new variable will contain the derivative of the log likelihood with
respect to the first regression equation.

{pmore}
The second new variable will contain the derivative of the log likelihood with
respect to the second regression equation.

{pmore}
The following variables will contain the derivatives of the log likelihood with
respect to each of the threshold parameters, the dependency parameter and the mixture parameters (if using mixtures).


{marker examples}{...}
{title:Examples}

{pstd}Predicting the probability of y1=1 and y2=2 and storing it in the
variable {opt pr1}

{phang}{cmd:. bicop y1 y2 x1 x2 x3, copula(clayton)}{p_end}
{phang}{cmd:. predict pr1,pr outcome(1,2)}{p_end}

{pstd}Predicting the probability of y1=2  and storing it in the variable
{opt pr2}

{phang}{cmd:. bicop y1 y2 x1 x2 x3, copula(clayton)}{p_end}
{phang}{cmd:. predict pr2,pr outcome(2,.)}{p_end}

{pstd}Predicting the probability of y1=2 conditional on y2=3  and storing it
in the variable {opt pr3}

{phang}{cmd:. bicop y1 y2 x1 x2 x3, copula(clayton)}{p_end}
{phang}{cmd:. predict pr3,pcpnd1 outcome(2,3)}{p_end}


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
Help:  {helpb bicop}
{p_end}
