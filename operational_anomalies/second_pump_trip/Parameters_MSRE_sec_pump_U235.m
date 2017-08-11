% MSRE 
format longe
%clear;
close all; clc;
% Perturbations
% SOURCE INSERTION
% No source insertion
sourcedata = [0 0 0];
sourcetime = [0 50 100];
% % 1 (n/no)/s for 10 seconds
% sourcedata = [0 10 0];
% sourcetime = [0 10 20];
source = timeseries(sourcedata,sourcetime);

% REACTIVITY INSERTION
% No reactivity insertion
simtime = 100;
reactdata = [0 0E-4];
reacttime = [0 1500];
% Periodic 60 PCM for 50 seconds
% simtime = 500;
% periodic = [0, 0; 50, 6e-4; 100, 0; 150, -6e-4; 200, 0; 250, 6e-4; 300, 0; 350, -6e-4; 400, 0]; 
% reactdata = periodic(:,2);
% reacttime = periodic(:,1);
% Step up 60 pcm 
% simtime = 1000;
% reactdata = [0 6e-3];
% reacttime = [0 300];
% % Step down -60 pcm for 10 sec
% simtime = 100;
% reactdata = [0 -6e-4];
% reacttime = [0 50];
% % Pulse 600 pcm for 0.1 sec
% simtime = 30;
% reactdata = [0 6e-3 0];
% reacttime = [0 10 10.1];




react = timeseries(reactdata,reacttime);

ts_max = 1e-1; % maximum timestep (s) *enter 'ts_max' in field 'max step size' in 'solver' pane of model configuration parameters (ctrl+E) within simulink to use this parameter

%% Pump trip/New flow rate fraction

% Time constant for exopnial deacy of flow (coastdown)
tau_p = 2;
therm_circ = .01;

%% NEUTRONICS DATA FOR U-235 VALUES
tau_l  = 16.73; % ORNL-TM-0728 %16.44; % (s)
tau_c  = 8.46; % ORNL-TM-0728 %8.460; % (s)
P      = 8;%8; % Thermal Power in MW ORNL-TM-1070, p.2
n_frac0 = 1; % initial fractional nuetron density n/n0 (n/cm^3/s)
Lam  = 2.400E-04;  % mean generation time ORNL-TM-1070 p.15 U235
% Lam    = 4.0E-04;  % mean generation time ORNL-TM-1070 p.15 U233
lam    = [1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00];
beta = [0.000223,0.001457,0.001307,0.002628,0.000766,0.00023]; % U235
% beta   = [0.00023,0.00079,0.00067,0.00073,0.00013,0.00009]; % U233
beta_t = sum(beta); % total delayed neutron fraction MSRE
rho_0  = beta_t - bigterm(beta,lam,tau_l,tau_c); % reactivity change in going from stationary to circulating fuel
C0(1)  = ((beta(1))/Lam)*(1.0/(lam(1) - (exp(-lam(1)*tau_l) - 1.0)/tau_c));
C0(2)  = ((beta(2))/Lam)*(1.0/(lam(2) - (exp(-lam(2)*tau_l) - 1.0)/tau_c));
C0(3)  = ((beta(3))/Lam)*(1.0/(lam(3) - (exp(-lam(3)*tau_l) - 1.0)/tau_c));
C0(4)  = ((beta(4))/Lam)*(1.0/(lam(4) - (exp(-lam(4)*tau_l) - 1.0)/tau_c));
C0(5)  = ((beta(5))/Lam)*(1.0/(lam(5) - (exp(-lam(5)*tau_l) - 1.0)/tau_c));
C0(6)  = ((beta(6))/Lam)*(1.0/(lam(6) - (exp(-lam(6)*tau_l) - 1.0)/tau_c));
% Feedback co-efficients
a_f    = -8.71E-05; % U235 (drho/°C) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
a_g    = -6.66E-05; % U235 (drho/°C) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101

%% CORE HEAT TRANSFER PARAMETERS
% FUEL PARAMETERS - DONE
vdot_f   = 7.5708E-02; % ORNL-TM-0728 % 7.571e-2; % vol. flow rate (m^3/s) ORNL-TM-1647 p.3, ORNL-TM-0728 p.12
rho_f    = 2.14647E+03; % (partially enriched U-235)ORNL-TM-0728 p.8 2.243E+03; % (Th-U) density of fuel salt (kg/m^3) ORNL-TM-0728 p.8
W_f      = 1.83085e+02;%vdot_f*rho_f; % 182.78; % calcd from m_dot*cp*delT=P; vdot_f*rho_f; % fuel flow rate (kg/s)
% tau_f_c  = tau_c; % ORNL-TM-0728 % 8.45; % transit time of fuel in core (s) ORNL-TM-1070 p.15, TDAMSRE p.5
m_f      = W_f*tau_c; % fuel mass in core (kg)
nn_f     = 2; % number of fuel nodes in core model
mn_f     = m_f/nn_f; % fuel mass per node (kg)
% cp_f     = 4.2*9/5; % (MJ/deg-C) total fuel heat capacity TDAMSRE p.5
scp_f    = 1.9665E-3;% specific heat capacity of fuel salt (MJ/kg-C) ORNL-TM-0728 p.8


% Core Upflow - DONE
v_g      = 1.95386; % graphite volume(m^3) ORNL-TM-0728 p. 101
rho_g    = 1.860E3; % graphite density (kg/m^3) ORNL-3812 p.77, ORNL-TM-0728 p.87
m_g      = v_g*rho_g; % graphite mass (kg)
cp_g     = 3.6*9/5; % TDAMSRE p.5 graphite total heat capacity (MW-s/C) ORNL-TM-1647 p.3
scp_g    = cp_g/m_g; % graphite specific heat capacity (MW-s/kg-C) ORNL-TM-1647 p.3
mcp_g1   = m_g*scp_g; % (mass of material x heat capacity of material) of graphite per lump (MW-s/°C)
mcp_f1   = mn_f*scp_f; % (mass of material x heat capacity of material) of fuel salt per lump (MW-s/°C)
mcp_f2   = mn_f*scp_f; % (mass of material x heat capacity of material) of fuel salt per lump (MW-s/°C)
hA_fg    = 0.02*9/5; % (fuel to graphite heat transfer coeff x heat transfer area) (MW/°C) ORNL-TM-1647 p.3, TDAMSRE p.5
k_g      = 0.07; % fraction of total power generated in the graphite  ORNL-TM-0728 p.9
k_1      = 0.5; % fraction of heat transferred from graphite which goes to first fuel lump
k_2      = 0.5; % fraction of heat transferred from graphite which goes to second fuel lump
k_f      = 0.93;     % fraction of heat generated in fuel - that generated in external loop ORNL-TM-0728 p.9
k_f1     = k_f/nn_f; % fraction of total power generated in lump f1
k_f2     = k_f/nn_f; % fraction of total power generated in lump f2


% % New node for power deposited in fuel outside the core
% k_out  = 1 - (k_g + k_f); % fraction of power generated in fuel in external loop ORNL-TM-0728 p.9
% m_out  = W_f; % (kg) Mass of node such that resident time is 1 sec

% Initial conditions - DONE
Tf_in  = 6.3222E+02; % in °C ORNL-TM-1647 p.2
T0_f2  = 6.5444E+02; % in °C 6.461904761904777e+02; ORNL-TM-1647 p.2
T0_f1  = Tf_in + (T0_f2-Tf_in)/2; % 6.405952380952389e+02; in °C
T0_g1  = T0_f1+(k_g*P/hA_fg); % 6.589285714285924e+02; in °C 
%T0_out = k_out*P/m_out/scp_f + T0_f2; % in °C

%% Heat Exchanger - DONE
% Geometry
d_he = 16; %(in) he diameter ORNL-TM-0728 p. 164
h_he = 72; % (in) active height % 96; %(in) he height ORNL-TM-0728 p. 164
od_tube = 0.5; %(in) coolant tube OD ORNL-TM-0728 p. 164
id_tube = od_tube - 2*0.042; % (in) coolant tube ID ORNL-TM-0728 p. 164
n_tube = 159; % number of coolant tubes ORNL-TM-0728 p. 164
a_tube = 254*144; % (in^2) total area of tubes ORNL-TM-0728 p. 164
l_tube = a_tube/n_tube/(pi*od_tube); % (in) tube length
v_tube = n_tube*pi*(od_tube/2)^2*l_tube; % (in^3) hx shell volume occupied by tubes
v_cool = n_tube*pi*(id_tube/2)^2*l_tube; % (in^3) hx volume occupied by coolant
v_he = (d_he/2)^2*pi*h_he; % (in^3) volume of heat exchanger shell
v_he_fuel = v_he-v_tube; % (in^3) volume available to fuel in shell

% Unit conversions
in_m = 1.63871e-5; % 1 cubic inch = 1.63871e-5 cubic meters

% PRIMARY FLOW PARAMETERS - DONE
W_p  = W_f; % fuel flow rate (kg/s)

m_p  = v_he_fuel*in_m*rho_f; % fuel mass in PHE (kg) 
nn_p = 4; % number of fuel nodes in PHE
mn_p = m_p/nn_p; % fuel mass per node (kg)
cp_p = scp_f; % fuel heat capacity (MJ/(kg-C))

% SECONDARY FLOW PARAMETERS - DONE
vdot_s   = 5.36265E-02; % ORNL-TM-0728 p. 164 % 5.236E-02; % coolant volume flow rate (m^3/s) ORNL-TM-1647 p.3
rho_s    = 1.922e3; % coolant salt density (kg/m^3) ORNL-TM-0728 p.8
W_s      = 1.0827e+02;%vdot_s*rho_s; % calcd from mdot*cp*delT; vdot_s*rho_s; % coolant flow rate (kg/s) ORNL-TM-1647 p.3 

m_s      = v_cool*in_m*rho_s; % coolant mass in PHE (kg) 
nn_s     = 4; % number of coolant nodes in PHE
mn_s     = m_s/nn_s; % coolant mass per node (kg)
scp_s    = 2.217E-3;%cp_s/m_s; %2.219E-03; % specific heat capacity of coolant (MJ/(kg-C) ORNL-TM-0728 p.8

A_phe    = 2.359E+01; % effective area for heat transfer (primary and secondary, m^2) ORNL-TM-0728 p.164

ha_p     = 6.480E-01; % heat transfer*area coefficient from primary to tubes (MW/C) ORNL-TM-1647 p.3
ha_s     = 3.060E-01; % heat transfer*area coefficient from tubes to secondary (MW/C) ORNL-TM-1647 p.3

% Primary Side 
mcp_pn   = mn_p*cp_p; % (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
hA_pn    = ha_p/nn_s; %3.030; % (primary to tube heat transfer coeff x heat transfer area) in MW/°C

% Tubes - DONE
nn_t     = 2; % number of nodes of tubes in model
rho_tube = 8.7745E+03; % (kg/m^3) density of INOR-8 ORNL-TM-0728 p.20
m_tn     = (v_tube - v_cool)*in_m*rho_tube/nn_t; % mass of tubes (kg)                                 
scp_t    = 5.778E-04; % specific heat capacity of tubes (MJ/(kg-C)) ORNL-TM-0728 p.20      
mcp_tn   = m_tn*scp_t; % from ratio of (A_phe/mcp_t)msbr = (A_phe/mcp_t)msre m_tn*cp_tn; % mass*(heat capacity) of tubes per lump in MW-s/°C      

% Secondary Side - DONE
mcp_sn   = mn_s*scp_s; % (mass of material x heat capacity of material) of coolant salt per lump in MW-s/°C
hA_sn    = ha_s/nn_s; % (tube to secondary heat transfer coeff x heat transfer area) in MW/°C

% Initial conditions -  DONE
% Primary nodes
Tp_in    = T0_f2; % 6.628E+02; % in °C ORNL-TM-1647 p.2
T0_p4    = Tf_in; % in °C ORNL-TM-1647 p.2
T0_p1    = Tp_in - (Tp_in-T0_p4)/4; % in °C
T0_p2    = Tp_in - 2*(Tp_in-T0_p4)/4; % in °C
T0_p3    = Tp_in - 3*(Tp_in-T0_p4)/4; % in °C

% Secondary nodes
Ts_in    = 5.4611E+02; % in °C ORNL-TM-1647 p.2
T0_s4    = 5.7944E+02; % in °C ORNL-TM-1647 p.2
T0_s1    = Ts_in + (T0_s4-Ts_in)/nn_s; % in °C
T0_s2    = Ts_in + 2*(T0_s4-Ts_in)/nn_s; % in °C
T0_s3    = Ts_in + 3*(T0_s4-Ts_in)/nn_s; % in °C
% Tube nodes
T0_t1    = (T0_p1*hA_pn+T0_s3*hA_sn)/(hA_pn+hA_sn); % in °C
T0_t2    = (T0_p3*hA_pn+T0_s1*hA_sn)/(hA_pn+hA_sn); % in °C
%% Radiator Parameters - NOT DONE

% Initial conditions -  DONE
% Primary nodes
Trp_in   = T0_s4; %5.933E+02; % in °C ORNL-TM-1647 p.2
T0_rp    = Ts_in; % in °C ORNL-TM-1647 p.2

% Secondary nodes -  DONE
Trs_in   = 37.78; % (C) air inlet temperature ORNL-TM-1647 p.2
T0_rs    = 148.9; % (C) air exit temperature ORNL-TM-1647 p.2

% Radiator Geometry
od_rad   = 0.01905; % (m) outer diameter of tubes in radiator ORNL-TM-0728 p.296
tube_wall_thick = 0.0018288; % (m) thickness of tubes in radiator ORNL-TM-0728 p.296
id_rad   = od_rad-2*tube_wall_thick;
n_rtubes = 120; %  number of tubes in radiator (rows times tubes per row) ORNL-TM-0728 p.296
l_rtube  = 9.144; % (m) length of tubes in radiator ORNL-TM-0728 p.296
v_rp     = pi*(id_rad/2)^2*l_rtube*n_rtubes; % volume available to salt in radiator
%v_rtube = pi*(od_rad/2)^2*l_rtube*n_rtubes-v_rp; % volume of metal in radiator tubes *TUBES NOT MODELED

n_tpr = 12; % number of tubes per row in radiator matrix
n_row = 10; % number rows in radiator matrix
tube_space = 0.0381; % (m) spacing between tubes and rows of matrix
v_rs = (n_row*od_rad+(n_row-1)*tube_space)*(n_tpr*od_rad+(n_tpr-1)*tube_space)*l_rtube; % volume of air inside radiator

% PRIMARY FLOW PARAMETERS - DONE
W_rp  = W_s; % coolant salt flow rate (kg/s)
m_rp  = v_rp*rho_s; % coolant salt mass in rad (kg)                     
nn_rp = 1; % number of coolant salt nodes in radiator
mn_rp = m_rp/nn_rp; % coolant mass per node (kg)
cp_rp = scp_s; % coolant specific heat capacity (MJ/(kg-C)) ORNL-TM-0728 p.8

% SECONDARY FLOW PARAMETERS - DONE
vdot_rs    = 94.389; % ORNL-TM-0728 p. 296; 78.82; % air volume flow rate (m^3/s) ORNL-TM-1647 p.2
rho_rs     = 1.1237; % air density (kg/m^3) REFPROP (310K and 0.1MPa)
W_rs       = vdot_rs*rho_rs; % air flow rate (kg/s) 

m_rs       = v_rs*rho_rs; % coolant air mass in rad (kg) 
nn_rs      = 1; % number of coolant nodes in rad
mn_rs      = m_rs/nn_rs; % coolant mass per node (kg)
scp_rs     = 1.0085E-3; % (MJ/kg-C) specific heat capacity of air at (air_out+air_in)/2 REFPROP

A_rad    = 6.503E1; % (m^2) surface area of radiator ORNL-TM-0728 p.14
h_roverall = P/A_rad/((T0_rp+Trp_in)/2-(T0_rs+Trs_in)/2); % cald as: P/A_rad/((T0_rp+Trp_in)/2-(T0_rs+Trs_in)/2)  3.168E-4; % (MW/m^2-C) polimi thesis

% Primary Side 
mcp_rpn   = mn_rp*cp_rp; % (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
hA_rpn    = h_roverall*A_rad/nn_rs; %3.030; % (primary to secondary heat transfer coeff x heat transfer area) in MW/°C

% Secondary Side - DONE
mcp_rsn   = mn_rs*scp_rs; % (mass of material x heat capacity of material) of coolant salt per lump in MW-s/°C
hA_rsn = h_roverall*A_rad/nn_rs; % (tube to secondary heat transfer coeff x heat transfer area) in MW/°C

%% Pure time delays between components - DONE

tau_hx_c = 8.67; % (sec) delay from hx to core TDAMSRE p.6
tau_c_hx = 3.77; % (sec) subtracted 1 sec for external loop power generation node resident time; delay from core to fuel hx TDAMSRE p.6
tau_hx_r = 4.71; % (sec) fertile hx to core TDAMSRE p.6
tau_r_hx = 8.24; % (sec) core to fertile hx TDAMSRE p.6


%% heat transfer coeffeint lookup table

xdat=[0.505782322224352
3.95959209339203
7.09516909533313
9.91333351559723
12.1007737102551
14.5990649861935
16.4674522240752
18.0200672553791
19.8851737430626
22.6910353500834
25.4968969571042
28.3052191267736
31.1143614839927
34.2359952975914
37.3584492987396
41.7325095005058
46.4198813462012
49.8589277414769
52.9846624928234
56.4237088880991
59.8627552833748
64.2417366104382
67.9949148371928
72.0597643327774
75.8121223719824
79.5644804111874
83.9442819258004
87.3849686961752
91.4522787544085
94.8921453372337
98.0195204636793
100];

ydat=[3.55277906881373
3.87511277578805
4.97990540503595
6.71050113459278
8.90942395494437
11.735790250704
14.8738278152938
18.3243568362632
22.0890176887115
26.1694507477376
30.2498838067638
33.8603493998961
37.3141591710638
41.0821007737103
44.6933865543921
49.2478880170599
53.9598654892419
57.1020039915794
60.0866664844028
63.2288049867403
66.3709434890778
69.9855100199579
73.128468709845
76.7422150531755
80.0418295650272
83.3414440768789
86.7993547857943
89.6281816442026
92.7719605216393
95.7574432020122
98.4287940509063
100];

xdat=0.01*xdat;
ydat=0.01*ydat;