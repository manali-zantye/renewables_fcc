$set DIM %gams.user1%
$set num_scen %gams.user2%
$offdigit

Set
   n             'nodes'            / 1*%DIM% /
   t             'time periods'     / t1*t%gams.user3% /
   stage(n,t)    'stage mapping'
   ancestor(n,n) 'ancestor matrix'
   leaf(n)       'leaf nodes';

stage('1','t1') = yes;
stage( n  ,'t2') = ord(n) > 1 and ord(n) <= %num_scen%+1;
stage( n  ,'t3') = ord(n) > %num_scen%+1 and ord(n) <= 2*%num_scen%+1;
stage( n  ,'t4') = ord(n) > 2*%num_scen%+1 and ord(n) <= 3*%num_scen%+1;

leaf(n) = stage(n,'t4');

ancestor(n,'1')$stage(n,'t2') = yes;
ancestor(n,n-%num_scen%)$stage(n,'t3') = yes;
ancestor(n,n-%num_scen%)$stage(n,'t4') = yes;

display stage
display leaf
display ancestor



* Random variables (demand) and probabilities

Parameter price(n) node demands
/
$include nodes_demands.gms
/
;

Parameter prob(n) node probabilities
/
$include nodes_probabilities.gms
/
;

* First stage
price('1') = 1;
prob('1')  = 1;


* Initialize ScenRed2
$set sr2prefix my_SR_problem
$libInclude scenred2

File fopts 'Scenred option file' / 'sr2%sr2prefix%.opt' /;
putClose fopts 'order           1'
             / 'section   epsilon'
             / ' 2           0.05'
             / ' 3           0.05'
             / 'end';

* Scenred2 method choice
ScenRedParms('construction_method') = 2;
ScenRedParms('reduction_method'   ) = 2;

$include SR_category.gms

ScenRedParms('out_scen') = 1;

ScenRedParms('sroption'           ) = 1;


$libInclude runscenred2 %sr2prefix% scen_red n ancestor prob ancestor prob price

display prob
display ancestor

Set srn(n) 'subset of nodes after tree construction';
srn(n) = prob(n);

display srn

File results /results.csv/;
put results;

loop(n$(srn(n)$(ord(n)>1)), put ord(n):<12:6/);
loop(n$(srn(n)$(ord(n)>1)), put prob(n):<12:6/);
loop(n$(srn(n)$(ord(n)>1)), put price(n):<12:6/);

*****Post processing******
*execute_unload "results_scenred.gdx";
*execute "gdx2sqlite -i results.gdx -o results_scenred.db";
