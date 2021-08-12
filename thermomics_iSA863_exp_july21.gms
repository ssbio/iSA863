
$INLINECOM /*  */
$onlisting
$offdigit

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



        exp(jSA) rxns with expression data available
/
$include "rxn_with_exp.txt"
/

        thermo(jSA) rxns with thermodynamic data available
/
$include "rxn_with_dG.txt"
/

        finitebounds(jSA)   RXNS WITH FINITE UPPER OR LOWER BOUNDS
/
$include "finitebounds.txt"
/

        iSA   Metabolites of SA
/
$include "metabolites.txt"
/

	medium_SA(jSA) Exchange rxns corresponding to compounds in the growth medium (SA)
$include "medium.txt"

* Allowing Amino acids to capture Paul fey's results
* https://mbio.asm.org/content/8/1/e01434-16
* "Amino Acid Catabolism in Staphylococcus aureus and the Function of Carbon Catabolite Repression"
* Cortney R. Halsey, Shulei Lei, Jacqueline K. Wax, Mckenzie K. Lehman, Austin S. Nuxoll, Laurey Steinke, Marat Sadykov, Robert Powers, Paul D. Fey
* mBio, 2017
* media information from  http://www.microbiologyresearch.org/docserver/fulltext/jmm/34/3/medmicro-34-3-143.pdf?expires=1544116213&id=id&accname=sgid026404&checksum=3FCAED918F485B2948ED5BF590E399F7

    AA_additional(jSA) 
/
$include "AA_additional.txt"
/

* The following sections are for model curation and 
* can be omitted at the final stage of model reconstruction
$ontext
	modelfix_off(jSA)	Reactions turned off for model curation
$include "modelfix_off.txt"

	modelfix_irrev(jSA)	Reactions made irreversible for model curation
*contains model fixes for cycles as well as GNG and NGG fixes by Growmatch
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

	modelfix_transport_off(jSA)	 Transport reactions turned off for model curation
$include "modelfix_transport_off.txt"

	modelfix_transport_irrev(jSA)	Transport reactions made irreversible for model curation
$include "modelfix_transport_irrev.txt"

	modelfix_transport_reverse(jSA)	Transport reactions with reveresed directionality for model curation
$include "modelfix_transport_reverse.txt"


    GrowmatchGNGfixes(jSA)    Reaction to be turned off  to fix GNGs (identified by Growmatch) 
/
*direction changes are included in the modelfix section
$include "rxn_GrowmatchGNGresults.txt"
/
$offtext

    regulation(jSA)     Reactions to be turned off for regulation effect
* For anaerobic regulations, turn oxygen off.
/
$include "downregulated_reactions.txt"
/
;


PARAMETERS
  UBSA(jSA) UB on rxns of SA 
  LBSA(jSA) LB on rxns of SA

  WT_flux(jSA)   WT aerobic flux distribution
/
$include "WT_flux.txt"
/


  rxntype_SA(jSA) Reaction type for rxns of SA 
/
$include "rxntype.txt"
/
;


******************  Set the values of LB and UB ******************
SCALAR Vmax /1000/;
SCALAR M /10000/;
SCALAR lambda /1/;

* Change the following parameters based on what condition and/or mutant is being simulated
SCALAR deltaG_avg /-18.2/;

* range of positive and negative deltaG value within which we allow for both directions
SCALAR deltaG_threshold /4.9/;

SCALAR flux_total /3500000/; 

* total expression value
* standard
SCALAR exp_total /35000/; 
* acetate
* SCALAR exp_total /3700/; 


SCALAR exp_avg /4920/;

SCALAR exp_max /9975/;

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

* constraint biomass, current value 0.164943  
LBSA('Sa_biomass_CDMG') = 0.16;
*UBSA('Sa_biomass_CDMG') = 0.165;



* ATPM
* 7.51 is the minimum in CDMG for a feasible solution. (so that the atpA mutant can grow)
* 5.00 is the minimum in CDM for a feasible solution.
LBSA('ATPM') = 7.51;
UBSA('ATPM') = 7.51;

*Glucose uptake
UBSA('EX_u_glc(e)') = 10;

*Gunanine uptake
UBSA('EX_u_gua(e)') = 2;

*Nicotinate uptake
UBSA('EX_u_nac(e)') = 1;

*oxygen uptake
LBSA('EX_u_o2(e)') = 10;

*L-leucine uptake
UBSA('EX_u_leu-L(e)') = 2;

*Thiamin uptake
UBSA('EX_u_thm(e)') = 1;


* Allowing Amino acids to capture Paul fey's results
* https://mbio.asm.org/content/8/1/e01434-16
* "Amino Acid Catabolism in Staphylococcus aureus and the Function of Carbon Catabolite Repression"
* Cortney R. Halsey, Shulei Lei, Jacqueline K. Wax, Mckenzie K. Lehman, Austin S. Nuxoll, Laurey Steinke, Marat Sadykov, Robert Powers, Paul D. Fey
* mBio, 2017
* media information from  http://www.microbiologyresearch.org/docserver/fulltext/jmm/34/3/medmicro-34-3-143.pdf?expires=1544116213&id=id&accname=sgid026404&checksum=3FCAED918F485B2948ED5BF590E399F7
UBSA(jSA)$AA_additional(jSA) = 5;



* regulation(s)
$ontext
Loop(jSA$regulation(jSA),
    UBSA(jSA) = 0.5*WT_aerobic_flux(jSA);
    );
$offtext


$ontext
* Reactions turned off for model curation
UBSA(jSA)$modelfix_off(jSA) = 0;
LBSA(jSA)$modelfix_off(jSA) = 0;

* Reactions made irreversible for model curation
LBSA(jSA)$modelfix_irrev(jSA) = 0;


* Reactions reveresed for model curation
UBSA(jSA)$modelfix_reverse(jSA) = 0;
LBSA(jSA)$modelfix_reverse(jSA) = -1000;

* Transport reactions turned off for model curation
UBSA(jSA)$modelfix_transport_off(jSA) = 0;
LBSA(jSA)$modelfix_transport_off(jSA) = 0;

* Transport reactions made irreversible for model curation
LBSA(jSA)$modelfix_transport_irrev(jSA) = 0;

* Transport reactions reveresed for model curation
UBSA(jSA)$modelfix_transport_reverse(jSA) = 0;
LBSA(jSA)$modelfix_transport_reverse(jSA) = -1000;

*turning of rxns to fix GNG inconsistencies identified by Growmatch
UBSA(jSA)$GrowmatchGNGfixes(jSA) = 0;
LBSA(jSA)$GrowmatchGNGfixes(jSA) = 0;
$offtext

****************test Section**********************************
$ONTEXT
* TCA cycle repression, in presence of Glucose, may not be needed for this formulation
UBSA('ACONTa') = 0.5;
UBSA('SUCOAS') = 0.5;
UBSA('AKGDH') = 0.5;
UBSA('CS') = 0.5;
UBSA('FUM') = 0.5;
UBSA('MDH3') = 0.5;
UBSA('CITL') = 0.5;
*UBSA('FRD2') = 0.5;
$OFFTEXT

* Forcing minimum flux through pta/ackA pathway
*UBSA('PTAr') = 0;
*LBSA('PTAr') = 0.5;

*** SSI grant mutants *****************
*pfkA
*UBSA('PFK') = 0;
*LBSA('PFK') = 0;

*pyc
*UBSA('PC') = 0;
*LBSA('PC') = 0;

*citZ
*UBSA('CS') = 0;
*LBSA('CS') = 0;

*sucC
*UBSA('SUCOAS') = 0;
*LBSA('SUCOAS') = 0;

*sucA
*UBSA('AKGDH') = 0;
*LBSA('AKGDH') = 0;

*ackA
*UBSA('ACKr') = 0;
*LBSA('ACKr') = 0;

*gudB
*UBSA('GLUDx') = 0;
*LBSA('GLUDx') = 0;

*ldhA
*UBSA('LDH_D') = 0;
*LBSA('LDH_D') = 0;

*UBSA('LDH_L') = 0;
*LBSA('LDH_L') = 0;

*ndhA
*UBSA('NADH10') = 0;
*LBSA('NADH10') = 0;

*UBSA('NADH7') = 0;
*LBSA('NADH7') = 0;

* menD
*UBSA('2S5EPAC') = 0;
*LBSA('2S5EPAC') = 0;

*atpA
*UBSA('ATPS24') = 0;
*LBSA('ATPS24') = 0;

*UBSA('ALAR') = 5;
*LBSA('ALAR') = 5;

****************stoichiometry*****************
PARAMETERS
  SSA(iSA,jSA) Stoichiometric matrix for SA 
/
$include "sij.txt"
/


    deltaG_min(jSA) min delta G values calculated for each reaction
/
$include "deltaG_min_iSA863.txt"
/

    rxnExp(jSA)  expression values calcualted for each reaction
/
$include "rxn_expression_WT.A.txt"
/
;


VARIABLES
	z     primal objective function (biomass_WT) for SA
	
	
    vSA(jSA)       Flux of SA rxns
;
	
vSA.lo(jSA)=LBSA(jSA);
vSA.up(jSA)=UBSA(jSA);
	
	
POSITIVE VARIABLES        
    epsp(jSA)  positive deviation variable for normalization
    epsn(jSA)  negative deviation variable for normalization
    
    xp(jSA)    positive deviation variable for flux
    xn(jSA)    negative deviation variable for flux
;

BINARY VARIABLES
    y(jSA)    Binary variables for reaction knockout or repression
;


*********************EQUATIONS**************************
EQUATIONS
	obj			    primal objective function (biomass_WT) 
    thermodynamics_linearization(jSA)   linearizing the absolute value in the objective function (thermodynamic normalization)
	flux_linearization(jSA)     linearizing the absolute value in the objective function (flux)
	massbalance(iSA)  	    Mass balance 
	ubound1(jSA)
	lbound1(jSA)
	ubound2(jSA)
	lbound2(jSA)
	ubound3(jSA)
	lbound3(jSA)
	
	universal_bound_exponential1(jSA)
	universal_bound_exponential2(jSA)
;


****************************** SA *******************************
* Primal 
*obj..      z =e= sum(j, (epsp(jSA) + epsn(jSA))) + lambda*sum(j, (xp(jSA) + xn(jSA)));


obj..       z =e= vSA('Sa_biomass_CDMG') - sum(jSA, (xp(jSA) + xn(jSA)));

*obj..       z =e= vSA('Sa_biomass_CDMG');

*thermodynamics_linearization(jSA)$thermo(jSA)..     epsp(jSA) - epsn(jSA) - vSA(jSA)/sum(j,vSA(jSA)) + (rxnExp(jSA)/exp_total)*(deltaG_min(jSA)/deltaG_avg) =e= 0;

thermodynamics_linearization(jSA)$(thermo(jSA) and exp(jSA))..     epsp(jSA) - epsn(jSA) - vSA(jSA)/flux_total + (rxnExp(jSA)/exp_total)*(deltaG_min(jSA)/deltaG_avg) =e= 0;

flux_linearization(jSA).. xp(jSA) - xn(jSA) - vSA(jSA) =e= 0;

massbalance(iSA)..     sum(jSA,SSA(iSA,jSA)*vSA(jSA)) =e= 0;

ubound1(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) < -deltaG_threshold)..   vSA(jSA) =l= UBSA(jSA)*(deltaG_min(jSA)/deltaG_avg)*(rxnExp(jSA)/exp_avg);

ubound2(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) > deltaG_threshold)..  vSA(jSA) =l= 0;


lbound1(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) < -deltaG_threshold)..  vSA(jSA) =g= 0;

lbound2(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) > deltaG_threshold)..  vSA(jSA) =g= LBSA(jSA)*(- deltaG_min(jSA)/deltaG_avg)*(rxnExp(jSA)/exp_avg);


ubound3(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) <= deltaG_threshold and deltaG_min(jSA) >= -deltaG_threshold)..   vSA(jSA) =l= UBSA(jSA)*(-abs(deltaG_min(jSA))/deltaG_avg)*(rxnExp(jSA)/exp_avg);


lbound3(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) <= deltaG_threshold and deltaG_min(jSA) >= -deltaG_threshold)..  vSA(jSA) =g= LBSA(jSA)*(-abs(deltaG_min(jSA))/deltaG_avg)*(rxnExp(jSA)/exp_avg);

*universal_bound_exponential1(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) < -deltaG_threshold)..  vSA(jSA) =l= UBSA(jSA)*(rxnExp(jSA)/exp_avg)*(1- system.exp(1.6774*(deltaG_min(jSA))));

*universal_bound_exponential2(jSA)$(thermo(jSA) and exp(jSA) and deltaG_min(jSA) > deltaG_threshold)..  vSA(jSA) =g= -LBSA(jSA)*(rxnExp(jSA)/exp_avg)*(1- system.exp(1.6774*(deltaG_min(jSA))));



universal_bound_exponential1(jSA)$(thermo(jSA) and (not finitebounds(jSA)) and exp(jSA) and (rxnExp(jSA) ne 0) and deltaG_min(jSA) < -deltaG_threshold)..  vSA(jSA) =l= UBSA(jSA)*(rxnExp(jSA)/exp_max)*(1- system.exp(1.6774*(deltaG_min(jSA))));

universal_bound_exponential2(jSA)$(thermo(jSA) and (not finitebounds(jSA)) and exp(jSA) and (rxnExp(jSA) ne 0) and deltaG_min(jSA) > deltaG_threshold)..  vSA(jSA) =g= -LBSA(jSA)*(rxnExp(jSA)/exp_max)*(1- system.exp(1.6774*(deltaG_min(jSA))));

******************************************************************************

file result /thermomics_results_iSA863_july21.txt/;
put result;

model thermomics
/
obj
flux_linearization

*ubound1
*lbound1
*ubound2
*lbound2
*ubound3
*lbound3

massbalance

universal_bound_exponential1
universal_bound_exponential2
/
;

thermomics.optfile=1;
Solve thermomics using lp maximizing z;



put 'primal model status:', thermomics.modelstat/;
put 'primal objective funtion:	', z.l:0:8 //;


loop (jSA,
	PUT ' 	', jSA.tl:25, "=", vSA.l(jSA):0:8/;
*		PUT vSA.l(jSA):0:8/;
	);
