option minlp = baron;

option optcr = 1e-4;
option optca = 1e-4;
$offdigit

*master sets
Set m "set of all energy sources/sinks" /cl, w, st, sp, g, a, d, c, b/;

Alias(m,mm);

*user input
Sets 
     t   "total time steps in optimization"    /1*%gams.user1%/,
     sor(m) "major energy sources"    /cl, w, st, sp/,
     sks(m) "major energy sinks" /g, a, d, c/
     ren(m) "renewables" /w, st, sp/
     ;
     
Sets
     JB(m,mm) "set of blocks mm receiving energy streams from block m" 
     /cl.(g,a,d,c),
      w.(g,a,b,c),
      b.d,
      st.(g,a,d,b,c),
      sp.(g,a,b,c)
     /;
     
Sets IB(m,mm) "set of blocks mm from which there is energy input to block m"
    /b.(w,st,sp),
     g.(cl, w, st, sp),
     a.(cl, w, st, sp),
     c.(cl, w, st, sp),
     d.(cl, b, st)
    /;


Positive variables 
            	CC_total "Total Capital cost of project ($)"
            	
            	CC(m) "cap cost in $"
            	
            CC_t "capital cost of solvent tanks ($)"
            CC_b "capital cost of electric boiler ($)"
            
            sz(m) "max cap in MW"
            nc_b "max capacity of electric boiler in MW"
            A_st "collector area of solar thermal technology in km2"
            
            Sa_max "rich/lean solvent tank volume in m3"

            g1(m,mm,t) "electrical energy flow in MW from block m to block mm at hour t"
            g_in(m,t) "total electrical energy flow in block m at hour t"
            g_out(m,t) "total electrical energy flow from block m at hour t"
            q_st(t) "solar thermal output in MW"

            OM_fixed
            om_fix_one(m)
            OM_var(t)
            
            Cg1(t) "Generation cost at hour t in dollars"
            out_tot
            E_gtot "we need this variable to ensure that total net co2 emissions are positive"
            ramp_cost(t) "ramping cost in hour t in dollars"
            
            	
	           ra(t) "CO2 adsorption rate in scrubber at hour t"
		   rd(t) "CO2 desorption rate in stripper at hour t"
		   rc(t) "CO2 compression rate at hour t"
		   Cg(t) "Generation cost per unit of gross output at hour t in dollars per MWh"
		   eeta_g(t) "Gross generation efficiency of coal plant at hour t"

		   Sa_final(t)	"Volume in rich solvent tank at time t in m3"
		   Sd_final(t)	"Volume in lean solvent tank at time t in m3"
		   eg(t) 	"CO2 emission intensity per unit of gross power output in tons/MWh"
		   Gen_cost(t)     "Electricity gen costs in $"
		   Storage_cost(t)  "CO2 storage and transport cost in $"
		   eor_sale(t) "revenue from selling captured co2 for enhanced oil recovery in $"
		   
	           S_tot(t)
		   DailyCO2int
		   
		   E_gN(t) "Net CO2 emissions at hour t in tons"
           E_gV(t) "CO2 emissions vented out before going to capture system at hour t in tons"
;


Parameter

Price_spot(t) "Hourly spot price in $/MWh" /
$include price.txt
/,

cf_renew(m,t)/
$include cap_fac.txt
/,


G_s(t) "Hourly direct normal irradiation of solar in W/m2"/
$include solar_rad.txt
/,

time_period(t) "time for which scenario t is valid"
/
$include time_periods.txt
/,

OM_fixed_one(m) "Unit operating and main cost of energy sources in $/MW.year"/
    w 0
    st 0
    sp 0
/,

OM_var_one(m) "variable o&m expenses in $/MWh"/
    st 3.5
/

CO(m) "unit renew cap cost ($/MW, $/m2, $/MW)"
/
w 
$include price_wt.gms
st 154
sp 
$include price_pv.gms
/
;


Scalar
    R_tax "Tax rate" /0.15/,
    r_disc "annual discount rate" /0.1/,
    t_dp "useful project life in years" /20/,
    t_lf "project lifetime in yrs" /25/,
    

    co_t "capital cost of one tank in $/m3" /300/,
    co_b "capital cost of one boiler in $/MW" /[88*1000]/,
    
    G_des "design point useful DNI for CSP in W/m2" /850/,
    eeta_rep "thermal to electrical conv efficiency for repowering case of CSP" /0.39/,
    eeta_te "thermal to electricity conversion factor" /0.183/,
    eeta_opt "CSP collector optical efficiency" /0.65/,
    
    c_1s "solar thermal heat loss coeff in W/m2K" /0/,
    c_2s "solar thermal heat loss coeff in W/m2K2" /0.0004/,
    T_rep_mean "csp mean T for rep case in C" /237.5/,
    T_boil_mean "csp mean T for reboiler case in C" /87.5/,
    T_rep_in "csp inlet T for rep case in C" /200/,
    T_boil_in "csp inlet T for reboiler case in C" /45/,
    delta_TS "T_ms - T_is: constant delta T for hot oil temp from CSP" /80/,
    L_st "land cost in $/m2" /[2554/4046.86]/,
    U_st "land utilization of solar thermal techn CHECK!! in %" /0.8/,
    
    L_pv "area occupied by Pv panels in m2 per MW of capacity" /[8*4046.86]/,
    
*    eor_price "price of CO2 in eor market ($/ton)" /35/
    
    Price_rampg "unit ramping cost in $/deltaMW" /2/,
    
    
	Price_genL "Average long term contract price in $/MWh" /51.7/

    Price_TS "Transport and storage costs of CO2 in $/ton: Mantripragada et al. (2019)" /4.6/

*Base case data
	mu_A0 "Base case efficiency penalty due to adsorption" /0.02/
	mu_B0 "Base case basic efficiency penalty" /0.01/
	mu_D0 "Base case efficiency penalty due to desorption" /0.04/
	mu_C0 "Base case efficiency penalty due to compression" /0.02/
	mu_G0 "Base case gross efficiency" /0.44/

	E0 	"Base case total power output in MW"
	alpha_A "Constant in net power output eqn for adsorption"
	alpha_D "Constant in net power output eqn for desorption"
	Cg0	"Base case electricity gen costs in $/MWh" /31/

	eg0_ref /0.76/


*Scrubber and Stripper data
	gamma_A "CO2 removal rate of scrubber" /0.90/
	gmin	 "Min power output in MW"
    lfmin "Min load fac for coal power plant" /0.2/
	ra_max	"Max rate of CO2 absorption" /1/
	rd_max	"Max rate of CO2 desorption" /1.25/
	delta_g_ramp	"Max ramping rate (MW/hr)" 
	delta_ra_max	 "Max ramping rate of scrubber (per hr)" /1/
	delta_rd_max	"Max ramping rate of stripper (per hr)" /1.25/
	e_gmax	"Baseline CO2 emission intensity in tons/MWh" /[300*0.001]/
	g_contract "Hourly Power output at hour t based on long term contract schedule in MW" 
	Rev_contract "Hourly Revenue from long term contracts in $"

*Wind power
	eeta_EB "Electricity to heat conversion efficiency of electric boiler" /0.96/
	
*	z_st /1/
    orig_scen "original number of total scenarios or time steps" /8760/
;

Scalar

g0 "Base case gross power output of coal plant in MW" /
$include namepcap.gms
/
	
e_g0	"Base case CO2 emissions in tons/MWh" /
$include co2eminten.gms
/

eor_price "price of CO2 in eor market ($/ton)" /
$include eor_price.gms
/
	
Price_CO2 /
$include tax_carbon.gms
/

cc_capt "total capital costs of capture in $"

co_capt "$/kW cost of capture in 2002 units(Rao,Rubin 07)" /810/

cepi_02 /395.6/
cepi_17 /567.5/
;

*Scalar assignments
E0 = g0/mu_G0;
alpha_A = E0*mu_A0;
alpha_D = E0*(mu_C0 + mu_D0);

delta_g_ramp = 0.3*g0;

g_contract = (2/3)*g0;
Rev_contract = g_contract*Price_genL;
gmin = lfmin*g0; 

cc_capt = co_capt*(cepi_17/cepi_02)*g0*(1e3)*(e_g0/eg0_ref);

Variables

	NPV "Net present value ($)"
	PF_net "Annual after tax operational profit ($)"
	PF_gro "annual gross operational profit ($)"
	PF(t) "Profit at scenario t"
	Rev_spot(t)     "Revenue from spot electricity price in $"
    CO2cost(t)
	profit
	Result_gross
	Result_net
	E_gN(t) "Net CO2 emissions at hour t in tons"
	E_gE(t) "Gross CO2 emissions from generation in tons"
	E_gS(t) "Amount of CO2 treated by scrubber in tons"
;

Binary variables
     z_st
     
     y_capt "for selection of capture system"
    ;

equations 
    objeq
    grossprof
	scenprof
	profiteq
	
	Spoteq
	Spoteq1(t)
	co2costeq
	Gencosteq
	Storageeq
	rampcosteq1(t) 
    rampcosteq2(t)
	omcostfixed
    omcostvar
    eorsaleeq(t)
	

    cap_tanks
	cap_boiler
	cap_boiler_max
	cap_st
	cap_total
    varvolbal1(t)
    varvolbal2(t)
	
	

	gineq(m,t)
	gouteq(m,t)
	boilerbal(t)
	boilermaxcap(t)
    maxtrans(t)


	st1(t)
	st5(t)
	gouteq_st(t)
	energybal(t)
    Totalvolads(t)
    Totalvoldes(t)
    CO2result1
    CO2result2
    sztoarea



	
	
	CO2rateeq(t)
	Costeq(t)
	Costeq1(t)
	Netemissioneq(t)
	Rampcons1(t)
	Rampcons2(t)
	Rampcons3(t)
	Rampcons4(t)
	Rampcons5(t)
	Rampcons6(t)

    emminteq(t)
	Grossemeq(t)
	Scrubemeq(t)
	Absbal(t)
	Compbal(t)
	Desbal(t)
	
	ventemeq(t)
	captsel1(t)
	captsel2(t)

    om_fixone(m)
    cap_ren(m)
    ren_pout(m,t)
    ;


objeq.. NPV =e= (-CC_total + (CC_total*R_tax/(r_disc*t_dp))*(1 -  1/(power((1+r_disc),t_dp))) + PF_net*((1/r_disc) - (1/(r_disc*(power((1+r_disc),t_lf))))))/(1e6);

scenprof(t)$(ord(t) ne card(t)).. PF(t) =e= Rev_spot(t) + Rev_contract + eor_sale(t) - CO2cost(t) - Gen_cost(t) - Storage_cost(t)  - ramp_cost(t)/time_period(t) - OM_var(t);

grossprof.. PF_gro =e= sum(t$(ord(t) ne card(t)),time_period(t)*PF(t)) - OM_fixed;

profiteq.. PF_net =e= (1 - R_tax)*PF_gro;

**********Operational profit eqs******************
Spoteq(t)$(ord(t) ne card(t)).. Rev_spot(t)  =e= (g_in('g',t) - g_contract)*Price_spot(t);

Spoteq1(t)$(ord(t) ne card(t)).. g_in('g',t) =g= g_contract;

CO2costeq(t)$(ord(t) ne card(t)).. co2cost(t) =e= E_gN(t)*Price_CO2;

Gencosteq(t)$(ord(t) ne card(t)).. Gen_cost(t) =e= Cg1(t);

Storageeq(t)$(ord(t) ne card(t)).. Storage_cost(t) =e=  E_gS(t)*Price_TS;

eorsaleeq(t)$(ord(t) ne card(t)).. eor_sale(t) =e=  E_gS(t)*eor_price;


rampcosteq1(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. ramp_cost(t) =g= Price_rampg*(g_out('cl',t+1) - g_out('cl',t));
rampcosteq2(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. ramp_cost(t) =g= Price_rampg*(g_out('cl',t) - g_out('cl',t+1));

omcostvar(t)$(ord(t) ne card(t)).. OM_var(t) =e= sum(m$ren(m),g_out(m,t)*OM_var_one(m));

om_fixone(m)$(ren(m)).. OM_fix_one(m) =e= OM_fixed_one(m)*sz(m);

omcostfixed.. OM_fixed =e= sum(m$ren(m),OM_fix_one(m));

**********Capital cost equations**************

cap_ren(m)$(ren(m) and (ord(m) ne 3) ).. CC(m) =e= CO(m)*sz(m);

cap_tanks.. CC_t =e= 2*Sa_max*CO_t;

cap_boiler.. CC_b =e= CO_b*nc_b;

cap_boiler_max.. nc_b =l= sum(m$ren(m), sz(m));

cap_st.. CC('st') =e= A_st*(1e6)*CO('st');

cap_total.. CC_total =e=  sum(m$ren(m), CC(m)) + cc_b + (cc_t + cc_capt)*y_capt;


*******Energy balance eqs***********************
gineq(m,t)$(ord(t) ne card(t)).. g_in(m,t) =e= sum(mm$IB(m,mm), g1(mm,m,t));

gouteq(m,t)$((ord(m) ne 3) and (ord(t) ne card(t))).. g_out(m,t) =e= sum(mm$JB(m,mm), g1(m,mm,t));

gouteq_st(t)$(ord(t) ne card(t)).. g_out('st',t) =e= z_st*(g1('st','g',t) + g1('st','a',t) + g1('st','c',t) + g1('st','b',t)) + (1 - z_st)*g1('st','d',t);

energybal(t)$(ord(t) ne card(t)).. sum(m$sor(m),g_out(m,t)) =e= sum(m$sks(m),g_in(m,t)) + g_in('b',t)*(1 - (eeta_EB*eeta_te));

*********Coal system eqs********************************************

Costeq1(t)$(ord(t) ne card(t)).. Cg(t)*eeta_g(t) =e= Cg0*mu_G0;

Costeq(t)$(ord(t) ne card(t)).. Cg1(t) =e= Cg(t)*g_out('cl',t);

Rampcons1(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. g_out('cl',t+1) - g_out('cl',t) =l= delta_g_ramp*time_period(t);
Rampcons2(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. g_out('cl',t+1) - g_out('cl',t) =g= -delta_g_ramp*time_period(t);

maxtrans(t)$(ord(t) ne card(t)).. g_in('g',t) =l= g0;

*******Boiler******
boilerbal(t)$(ord(t) ne card(t)).. g_out('b',t) =e= g_in('b',t)*(eeta_EB*eeta_te);

boilermaxcap(t)$(ord(t) ne card(t)).. g_in('b',t) =l= nc_b;

*********Renewables pout eqs*******************************
ren_pout(m,t)$(ren(m) and (ord(m) ne 3) and (ord(t) ne card(t))).. g_out(m,t) =e= cf_renew(m,t)*sz(m);

*****Solar thermal eqs*************
st1(t)$(ord(t) ne card(t)).. q_st(t) =e= A_st*(eeta_opt*G_s(t));

sztoarea.. sz('st') =e= A_st*G_des*eeta_opt*(eeta_rep*z_st + eeta_te*(1 - z_st));

st5(t)$(ord(t) ne card(t)).. g_out('st',t) =e= q_st(t)*eeta_rep*z_st + q_st(t)*eeta_te*(1 - z_st);


*****Capture system eqs***********

emminteq(t)$(ord(t) ne card(t)).. eg(t)*eeta_g(t) =e= e_g0*mu_G0;

Grossemeq(t)$(ord(t) ne card(t)).. E_gE(t) =e= g_out('cl',t)*eg(t);

Scrubemeq(t)$(ord(t) ne card(t)).. E_gS(t) =e= rc(t)*g0*e_g0*gamma_A;

Netemissioneq(t)$(ord(t) ne card(t)).. E_gN(t) =e= E_gE(t) - ra(t)*g0*e_g0*gamma_A;

CO2rateeq(t)$(ord(t) ne card(t)).. rd(t) =e= rc(t);

Rampcons3(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. ra(t + 1) - ra(t) =l= delta_ra_max*time_period(t);
Rampcons4(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. ra(t + 1) - ra(t) =g= -delta_ra_max*time_period(t);
Rampcons5(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. rd(t + 1) - rd(t) =l= delta_rd_max*time_period(t);
Rampcons6(t)$((ord(t) ne card(t)) and (ord(t) ne (card(t)-1))).. rd(t + 1) - rd(t) =g= -delta_rd_max*time_period(t);

Absbal(t)$(ord(t) ne card(t)).. E0*mu_A0*ra(t) =e= g_in('a',t);
Desbal(t)$(ord(t) ne card(t)).. E0*mu_C0*rc(t) =e= g_in('c',t);
Compbal(t)$(ord(t) ne card(t)).. E0*mu_D0*rd(t) =e= g_in('d',t);


******Solvent storage bal********
Totalvolads(t)$( ord(t) ne card(t) ).. Sa_final(t+1) =e= Sa_final(t) + (g0*7300*e_g0/(600*eg0_ref))*(ra(t) - rd(t))*time_period(t);
Totalvoldes(t)$( ord(t) ne card(t) ).. Sd_final(t+1) =e= Sd_final(t) + (g0*7300*e_g0/(600*eg0_ref))*(rd(t) - ra(t))*time_period(t);
varvolbal1(t).. Sd_final(t) =l= Sa_max;
varvolbal2(t).. Sa_final(t) =l= Sa_max;
Sa_final.fx('1') = (g0*7300*e_g0/(600*eg0_ref));
Sd_final.fx('1') = (g0*7300*e_g0/(600*eg0_ref));
Sa_final.fx(t)$(ord(t) eq card(t)) = (g0*7300*e_g0/(600*eg0_ref));
Sd_final.fx(t)$(ord(t) eq card(t)) = (g0*7300*e_g0/(600*eg0_ref));


*******Constraint on CO2 emissions*******

CO2result1.. E_gtot =e= sum(t$(ord(t) ne card(t)),E_gN(t)*time_period(t));

CO2result2.. out_tot =e= sum(t,g_in('g',t)*time_period(t));


*Upper and lower bounds

g_out.lo('cl',t) = gmin;
g_out.up('cl',t) = g0;

ra.lo(t) = 0;

rd.lo(t) = 0;


captsel1(t)$( ord(t) ne card(t) ).. ra(t) =l= ra_max*y_capt;
captsel2(t)$( ord(t) ne card(t) ).. rd(t) =l= rd_max*y_capt;


Sa_max.fx = 2*(g0*7300*e_g0/(600*eg0_ref));

*bounds on design vars
A_st.up = 3;
A_st.fx = 0;
z_st.fx = 1;

eeta_g.fx(t) = 0.44;

*New set of capture constraints
ventemeq(t)$(ord(t) ne card(t)).. E_gV(t) =e= E_gE(t) - ra(t)*g0*e_g0;



*fixing design decisions
$include "fix_design.gms";


option reslim = 172800;
option limrow = 60;

model flexible /all/;


solve flexible using MINLP maximizing NPV;


*****Post processing******
execute_unload "results.gdx";
execute "gdx2sqlite -i results.gdx -o results.db";
