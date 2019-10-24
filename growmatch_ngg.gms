****************************************************************************************
*                                     GrowMatch for NGGs                               *
*                                     ------------------                               *
* Related publications: PMID: 19282964 & PMID:21190580                                 *
*                                                                                      *                                        
*  - This code works directly with the genes, i.e., there is no need to manually       *
*    delete the rxns that correspond to a specific gene knockout (NGG). To this end,   
*    the reaciton fluxes are correlated with genes through using binary                *
*    variables/parameters.                          
*  - For this sample code, the inconsistencies are fixed by addings rxns from a        
*    set called 'database' which constains only the irreversible reactions in the  
*    reverse direction and uptake and transport reactions that are missing in the     
*    original model.                                                                   *
*                                                                                                 *
****************************************************************************************

$INLINECOM /*  */
$onlisting
$offdigit

$SET biomass R_Sa_biomass_universal

OPTIONS
        limrow = 1000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        lp = cplex
;


SETS
        i The set of metabolites in the model and those involved in the rxns of the Database
$include "metabolites.txt"
$include "metabolites_database.txt"
/* include the name of the file containing your own data */

        j The set of reactions in the original model and those in the Database
$include "reactions.txt"
$include "database_reactions.txt"
/* include the name of the file containing your own data */

        database(j) The set of reactions in database
$include "database_reactions.txt"
/* include the name of the file containing your own data */

        medium_rxns(j) Exchange fluxes for metabolites in the growth medium
$include "medium.txt"
/* or, you can include the name of the file containing your own data */

        g genes included in model
$include "gene_names.txt"
/* include the name of the file containing your own data */

        ngg(g)  The set of NGGs
$include "NGGs.txt"

* Dynamic sets used for program control and results output
        currentflux(j) flux that is currently being deleted
        blocked(j) blocked variables set
        essential(j) essential variables set

        iter /1*10000/
        eqncounter(iter) Dynamic set containing the current number of runs of loop
;

* At first there is nothing in the set eqncounter
eqncounter(iter)=no;

PARAMETERS
  UB(j) UB on reactions
  LB(j) LB on reactions

  maxbm Max biomass flux when no genes or rxns is knocked out

  S(i,j) Stoichiomteric matrix for rxns in the origianl model and those in the Database except those in the reverse direction of irreversible rxns
$include "Sij.txt"
$include "Sij_database"

 prev_w(database,iter) Stores the binary variables which are 1 or 0 in each iteration

 n_min  The minimum number of added rxns

 counter

 done  Represents when we are done with finding the rxns to resolve a NGG
 
 		rxntype(j) Rxn Types
$include "rxntype.txt"
$include "rxntype_database.txt"
;

* The stoichiomteric coefficients for the reverse direction of irreversible rxns
* This file can also incorporated into the file for S(i,j) if writtein the that format
*$INCLUDE "S_matrix_irrev.txt"

* Import the file containing the lowerbound and upperbound on rxns
*$INCLUDE "bounds_on_reactions.txt"



* Adjust the Viability threshold
SCALAR viabilityThr /0.01/;


** Set up the experimental conditions, e.g., for the Minimal medium: aerobic grwoth on glucose **
*starts with all upperbounds are a thousand so I don't need to mainly add UB for the newly added reactions.
*for all the backward reactions
* irreversible reactions (rxntype = 0)
UB(j)$(rxntype(j) = 0) = 1000;
LB(j)$(rxntype(j) = 0) = 0;

* For reversible reactions (rxntype = 1) 
UB(j)$(rxntype(j) = 1) = 1000;
LB(j)$(rxntype(j) = 1) = -1000;

* uptake reactions (rxntype = 3)
UB(j)$(rxntype(j) = 3) = 0;
LB(j)$(rxntype(j) = 3) = 0;

* export reactions (rxntype = 4)
UB(j)$(rxntype(j) = 4) = 1000;
LB(j)$(rxntype(j) = 4) = 0;


* growth medium 
UB(j)$medium_rxns(j) = 1000;


***** Set the conditions for the growth medium *****
UB(j)$medium_rxns(j) = 1000;

*Glucose uptake
UB('EX_u_glc_e') = 10;

*Gunanine uptake
UB('EX_u_gua_e') = 2;

*Nicotinate uptake
UB('EX_u_nac_e') = 1;

*oxygen uptake
UB('EX_u_o2_e') = 2;

*L-leucine uptake
UB('EX_u_leu__L_e') = 2;

*Thiamin uptake
UB('EX_u_thm_e') = 1;

* ATPM
UB('R_ATPM') = 7.51;
LB('R_ATPM') = 7.51;



*******************  Define the variables ************************

VARIABLES
        v(j)          Fluxes
        z_biomass     objective function
        z_growmatch
;

BINARY VARIABLES
        y(g)          zero if gene g is knocked out
        w(database)   One if reaction rxn 'database' is added to the model

* NOTE: If one wants to work with reactions instead of genes, then there is
* no need to define the binary variable y(g)

* We always fix the value of y(g) in this code to zero or one and also use the 'holdfixed'
* option of the GAMS. This implies that y(g) in fact  acts as a parameter

;

***************** Define the equations  ****************************
EQUATIONS
        massbalance(i)        mass balance equations for each metabolite i
        biomassobj            define objective function
        growtheq              Maintain growth
        LBeq(database)        Lowerbound on rxns
        UBeq(database)        Upperbound on rxns
        growmatchobj          Objective function for GrowMatch

* Equations for deleting the previous solutions
       integercut(iter)  This eqn is defined over a dynamic set containing the current number of iterations

* Define the name of equations for gene to rxn maps
$include "gene_rxn_map_eqn_names.txt"
;

massbalance(i)..           sum(j,S(i,j)*v(j)) =e= 0;
biomassobj..               z_biomass =e= v('%biomass%');
growtheq..                 v('%biomass%') =g= viabilityThr*maxbm;
LBeq(database)..           v(database) =g= LB(database)*w(database);
UBeq(database)..           v(database) =l= UB(database)*w(database);
growmatchobj..             z_growmatch =e= sum(database,w(database));

* The eqns defining gene to rxn maps
$include "gene_rxn_map_eqns.txt";


* Integer cuts to preclude the previously found solutions
integercut(eqncounter)..   sum(database$(prev_w(database,eqncounter) eq 1),w(database)) =l= n_min-1;

*************** Define the Models ********************
Model FBA_gene
/
 massbalance
 biomassobj
$include "gene_rxn_map_eqn_names.txt"
/;

Model GrowMatch
/
 massbalance
 LBeq
 UBeq
 growtheq
 growmatchobj
 integercut
$include "gene_rxn_map_eqn_names.txt"
/;

* Include the option files (stored in coplex.opt)
FBA_gene.optfile=1;
GrowMatch.optfile=1;

* Do not consider the fixed varibles when optimization model is created
FBA_gene.holdfixed=1;
GrowMatch.holdfixed=1;

* turn all the genes on
y.fx(g) = 1;

v.lo(j)=LB(j);
v.up(j)=UB(j);


*************** Solve the models  ************
** First find out the maximum biomass flux of the wild-type network
v.fx(database)=0;

SOLVE FBA_gene USING mip MAXIMIZING z_biomass;
maxbm = z_biomass.l;

* Now reset the bounds for the rxns in the set 'database'
v.lo(database)=LB(database);
v.up(database)=UB(database);

* Define a file to store the results
FILE resultfile  /resolved_nggs.txt/;
PUT resultfile;
PUT "";
resultfile.ap=1;

* Initializing prev_w for the first run of the loop
prev_w(database,iter)=no;
n_min=0;
done=0;
counter=0;


alias(g,g1);

* The model GrowMatch is now solved in a loop to find out alternative solutions

LOOP(g1$(ngg(g1)),
   y.fx(g1)=0;
   PUT resultfile;
   PUT /"********** ",g1.tl:0," ***************"/;
   eqncounter(iter)=no;
   SOLVE GrowMatch using mip minimizing z_growmatch;
   if(GrowMatch.modelstat = 10,       /*If the model is infeasible */
     PUT resultfile;
     PUT " Reactions in the set Dataset cannot fix this NGG! "/;
   elseif (GrowMatch.modelstat ne 1),  /* If a non-optimum solution is found */
     PUT resultfile;
     PUT "ERROR: modelstat= "GrowMatch.modelstat/; /* Refer to GAMS manual for the modelstat*/
   elseif (z_growmatch.l=0),   /* If it turns out that no rxns is needed */
     PUT resultfile;
     PUT "Objective = 0, No rxns need to be added! This is not a NGG"/;
   else     /* If at least one optimal solution exists */
      counter=0;
      n_min=z_growmatch.l;   /* store the minimal number of rxns needed */
      PUT "This NGG can be solved by adding  a minimum of ",n_min:0:1," rxns to the model"//;
      done=0;
      LOOP(iter$(not done),  /* Now find all alternative solutions */
         counter=counter+1;
         SOLVE GrowMatch using mip minimizing z_growmatch;
         if(GrowMatch.modelstat ne 1,
            done=1;
            PUT resultfile;
            PUT "ERROR: modelstat= "GrowMatch.modelstat,"  No optimal solution!"/;
         elseif (z_growmatch.l > n_min),   /* If number of rxns needed is greater than n_min */
            done=1;           /* Do not continue */
            PUT "obj > n_min = ",n_min:0:1/;
         else                /* Continue finding other alternative solutions */
            eqncounter(iter)=yes;  /* Add one more element to our dynamic set */
            prev_w(database,iter)=0;  /* Initiate prev_w for this iteration */
            PUT resultfile;
            PUT counter:0:0,") ";
            LOOP(database$(w.l(database) = 1),
                PUT " '",database.tl:0,"'";
                prev_w(database,iter)=1;  /* Store the soln at this iteration in prev_w */
           );
           PUT /;
         );
         if(counter = 10,      /* Find a maximum of 10 alternative solutions */
             done=1;
         );
      );
   );
   y.fx(g1)=1;
);

