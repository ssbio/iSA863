
$INLINECOM /*  */
$onlisting
$offdigit

* GAMS script for running Flux variability Analysis 
* on Genome-scale metabolic models.
* The models, from Systems-Biology Markup Language, 
* needs to be converted to GAMS-specific input files.
* - A list of reactions
* - A list of metabolites
* - A list of reaction types for each of the reactions 
*     - Irreversible reactions: 0
*     - Reversible reactions: 1
*     - Uptake reactions: 3
*     - Export reaction: 4
* - The stoichiometric matrix
* - The list of media components
* - List(s) of regulated reactions in differen conditions (optional)


* The following options can be overriden using the optfile
options limrow = 1000
	optCR = 1E-9
	optCA = 1E-9
	iterlim = 1000000
	decimals = 8
	reslim = 1000000
	work = 50000000
    lp = cplex;
        

SETS
        jSA   Reactions of SA 
/
$include "reactions.txt"
/

        iSA   Metabolites of SA
/
$include "metabolites.txt"
/

	medium_SA(jSA) Exchange rxns corresponding to compounds in the growth medium (SA)
$include "medium.txt"

* Additional Amino acids 
    aa_additional(jSA) 
/
$include "aa_additional.txt"
/


* The following sections are for model curation and 
* can be omitted at the final stage of model reconstruction
$ontext
	modelfix_off(jSA)	Reactions turned off for model curation
*contains model fixes for cycles as well as GNG fixes by Growmatch
/
$include "modelfix_off.txt"
$include "rxn_GrowmatchGNGresults.txt"
/

	modelfix_irrev(jSA)	Reactions made irreversible for model curation
/
$include "modelfix_irrev.txt"
$include "rxn_irrev_GrowmatchGNGresults.txt"
/
	modelfix_reverse(jSA)	Reactions with reversed directionality for model curation
*contains model fixes for cycles as well as GNG and NGG fixes by Growmatch
/
$include "modelfix_reverse.txt"
$include "rxn_reverse_GrowmatchGNGresults.txt"
$include "rxn_reverse_GrowmatchNGGresults.txt"
/

    regulation(jSA)     Reactions to be turned off for regulation effect
/
$include "downregulated_reactions.txt"
/
$offtext
;


PARAMETERS
  UBSA(jSA) UB on rxns of SA 
  LBSA(jSA) LB on rxns of SA


  rxntype_SA(jSA) Reaction type for rxns of SA 
/
$include "rxntype.txt"
/


******************  Set the values of LB and UB ******************
SCALAR Vmax /1000/;

*** Set the experimental conditions

* irreversible reactions (rxntype = 0)
UBSA(jSA)$(rxntype_SA(jSA) = 0) = Vmax;
LBSA(jSA)$(rxntype_SA(jSA) = 0) = 0;

* For reversible reactions (rxntype = 1) 
UBSA(jSA)$(rxntype_SA(jSA) = 1) = Vmax;
LBSA(jSA)$(rxntype_SA(jSA) = 1) = -Vmax;

* uptake reactions (rxntype = 3)
UBSA(jSA)$(rxntype_SA(jSA) = 3) = 0;
LBSA(jSA)$(rxntype_SA(jSA) = 3) = 0;

* export reactions (rxntype = 4)
UBSA(jSA)$(rxntype_SA(jSA) = 4) = Vmax;
LBSA(jSA)$(rxntype_SA(jSA) = 4) = 0;

* growth medium 
UBSA(jSA)$medium_SA(jSA) = 0;

* ATPM
* 7.51 is the minimum in CDMG for a feasible solution. (so that the atpA mutant can grow)
* 5.00 is the minimum in CDM for a feasible solution.
*LBSA('R_ATPM') = 7.51;
*UBSA('R_ATPM') = 7.51;

*$ontext
*Glucose uptake
*UBSA('EX_u_glc_e') = 10;
UBSA('EX_u_glc_e') = 0;

*Gunanine uptake
UBSA('EX_u_gua_e') = 2;

*Nicotinate uptake
UBSA('EX_u_nac_e') = 1;

*oxygen uptake
UBSA('EX_u_o2_e') = 10;

*L-leucine uptake
UBSA('EX_u_leu__L_e') = 2;

*Thiamin uptake
UBSA('EX_u_thm_e') = 1;

* Additional amino acid uptake
UBSA(jSA)$aa_additional(jSA) = 5;

*$offtext


$ontext
* Reactions turned off for model curation
UBSA(jSA)$modelfix_off(jSA) = 0;
LBSA(jSA)$modelfix_off(jSA) = 0;

* Reactions made irreversible for model curation
LBSA(jSA)$modelfix_irrev(jSA) = 0;

* Reactions reveresed for model curation
UBSA(jSA)$modelfix_reverse(jSA) = 0;
LBSA(jSA)$modelfix_reverse(jSA) = -1000;
$offtext

**************** mutant test Section**********************************
* Set the upper and lower bounds of the reactions catalyzed by 
* the gene(s) of interest to 0.0

************ Example *****************
*pfkA
*UBSA('PFK') = 0;
*LBSA('PFK') = 0;

****************stoichiometry*****************
PARAMETERS
  SSA(iSA,jSA) Stoichiometric matrix for SA 
/
$include "sij.txt"
/
;

PARAMETERS
C(jSA)
max(jSA)
min(jSA);
alias (j1SA, jSA);

VARIABLES
	zprimal_SA     primal objective function (biomass_WT) for SA
    vSA(jSA)       Flux of SA rxns
;
	
vSA.lo(jSA)=LBSA(jSA);
vSA.up(jSA)=UBSA(jSA);


*********************EQUATIONS**************************
EQUATIONS

	primalobj_SA			primal objective function (biomass_WT) for SA

	massbalance_SA(iSA)  	Mass balance for SA
;


****************************** SA *******************************

primalobj_SA..		  zprimal_SA =e= sum(jSA$(C(jSA) eq 1),vSA(jSA));

massbalance_SA(iSA)..     sum(jSA,SSA(iSA,jSA)*vSA(jSA)) =e= 0;


******************************************************************************
model primalSA
/
primalobj_SA
massbalance_SA
/
;

file result /FVA_results.txt/;
put result;

put 'FVA results for SA when no Carbon source is present'/;
put 'Rxn            Min             Max'/;
put '------------------------------------------------------------------------'/;

scalar count /0/;


LOOP(j1SA,
	C(jSA) = 0;
	C(j1SA) = 1;
	
primalSA.optfile=1;
Solve primalSA using lp maximizing zprimal_SA;
max(j1SA) = zprimal_SA.l;

Solve primalSA using lp minimizing zprimal_SA;
min(j1SA) = zprimal_SA.l;
put j1SA.tl:25, system.tab min(j1SA):0:8, system.tab max(j1SA):0:8/;
);
