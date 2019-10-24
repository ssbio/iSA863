
$INLINECOM /*  */
$onlisting
$offdigit

* GAMS script for running Essentiality Analysis 
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

* Set the biomass reaction ID
$SET biomass R_Sa_biomass_universal

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

* If regulation needs to be imposed for a mutant, 
* save the WT flux distributions in "WT_flux.txt" file
$ontext
  WT_flux(jSA)   WT aerobic flux distribution
/
$include "WT_flux.txt"
/
$offtext


  rxntype_SA(jSA) Reaction type for rxns of SA 
/
$include "rxntype.txt"
/



  S(iSA,jSA) contains the Stoichiometric matrix of the metabolic model
/
$include "sij.txt"
/


totalrxn total number of reactions
count
essanaero
nessaero
c
;

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
UBSA(jSA)$medium_SA(jSA) = 1000;

* ATPM
* In CDMG media, use 7.51 
* In CDM media, use 5.00
LBSA('R_ATPM') = 5.0;
UBSA('R_ATPM') = 5.0;

*Glucose uptake
* In CDM media, set to 0.0 
*UBSA('EX_u_glc_e') = 10;
UBSA('EX_u_glc_e') = 0.0;

*Gunanine uptake
UBSA('EX_u_gua_e') = 2;

*Nicotinate uptake
UBSA('EX_u_nac_e') = 1;

* Oxygen uptake
* For anaerobic regulations, turn oxygen off.
UBSA('EX_u_o2_e') = 10;

*L-leucine uptake
UBSA('EX_u_leu__L_e') = 2;

*Thiamin uptake
UBSA('EX_u_thm_e') = 1;

* Additional amino acid uptake
UBSA(jSA)$aa_additional(jSA) = 5;

$ontext

* Reactions turned off for model curation
UBSA(jSA)$modelfix_off(jSA) = 0;
LBSA(jSA)$modelfix_off(jSA) = 0;

* Reactions made irreversible for model curation
LBSA(jSA)$modelfix_irrevSA(jSA) = 0;

* Reactions reveresed for model curation
UBSA(jSA)$modelfix_reverse(jSA) = 0;
LBSA(jSA)$modelfix_reverse(jSA) = -1000;

* Transport reactions turned off for model curation
UBSA(jSA)$modelfix_transport_off(jSA) = 0;
LBSA(jSA)$modelfix_transport_off(jSA) = 0;

* Transport reactions made irreversible for model curation
LBSA(jSA)$modelfix_transport_irrevSA(jSA) = 0;

* Transport reactions reveresed for model curation
UBSA(jSA)$modelfix_transport_reverse(jSA) = 0;
LBSA(jSA)$modelfix_transport_reverse(jSA) = -1000;

*turning off rxns to fix GNG inconsistencies identified by Growmatch
UBSA(jSA)$GrowmatchGNGfixes(jSA) = 0;
LBSA(jSA)$GrowmatchGNGfixes(jSA) = 0;
$offtext
**********************

essanaero=0;
nessaero=0;

VARIABLES
        vSA(jSA)      Flux
        z         Objective function
;

vSA.lo(jSA)=LBSA(jSA);
vSA.up(jSA)=UBSA(jSA);

*set at a percentage of the max biomass at the given uptake condition
* to check for different percentages, change the value of vmaxbiom 
Parameter vmaxbiom/0.0015/;

Scalar n;
n=card(jSA);
Scalar x;

EQUATIONS
        massbal(iSA)   Mass balance equations for each metabolite i
        objSA              Objective function
        reactionflux     constraint for reaction flux to be set ot zero
        OxFlux       Fixing the oxygen uptake rate
;


objSA..		  z =e= vSA('%biomass%');
reactionflux..      sum(jSA$(ord(jSA) eq count),vSA(jSA))=e=0;
massbal(iSA)..  sum(jSA,S(iSA,jSA)*vSA(jSA) )=e=0;

Model Essentiality
/objSA
massbal
reactionflux
/;

FILE D /Essentiality.txt/;
PUT D;

for(x= 1 to n by 1,
         count=x;
         SOLVE Essentiality USING LP MAXIMIZING z;

if (z.l<vmaxbiom,
         essanaero=essanaero+1;
         loop(jSA$(ord(jSA) eq count),
         put jSA.tl:25/; 
         );
);
);


PUT  essanaero//  ;
