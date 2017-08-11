% MSRE 
% format longe
% clear; close all; clc;
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
reactdata = [0 2E-4];
reacttime = [0 2500];
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
%% NEUTRONICS DATA 
tau_l  = 16.73; % ORNL-TM-0728 %16.44; % (s)
tau_c  = 8.46; % ORNL-TM-0728 %8.460; % (s)
P      = .1; % Thermal Power in MW ORNL-TM-1070, p.2
n_frac0 = 1; % initial fractional nuetron density n/n0 (n/cm^3/s)
% Lam  = 2.400E-04;  % mean generation time ORNL-TM-1070 p.15 U235
Lam    = 4.0E-04;  % mean generation time ORNL-TM-1070 p.15 U233
lam    = [1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00];
% beta = [0.000223,0.001457,0.001307,0.002628,0.000766,0.00023]; % U235
beta   = [0.00023,0.00079,0.00067,0.00073,0.00013,0.00009]; % U233
beta_t = sum(beta); % total delayed neutron fraction MSRE
rho_0  = beta_t - bigterm(beta,lam,tau_l,tau_c); % reactivity change in going from stationary to circulating fuel
C0(1)  = ((beta(1))/Lam)*(1.0/(lam(1) - (exp(-lam(1)*tau_l) - 1.0)/tau_c));
C0(2)  = ((beta(2))/Lam)*(1.0/(lam(2) - (exp(-lam(2)*tau_l) - 1.0)/tau_c));
C0(3)  = ((beta(3))/Lam)*(1.0/(lam(3) - (exp(-lam(3)*tau_l) - 1.0)/tau_c));
C0(4)  = ((beta(4))/Lam)*(1.0/(lam(4) - (exp(-lam(4)*tau_l) - 1.0)/tau_c));
C0(5)  = ((beta(5))/Lam)*(1.0/(lam(5) - (exp(-lam(5)*tau_l) - 1.0)/tau_c));
C0(6)  = ((beta(6))/Lam)*(1.0/(lam(6) - (exp(-lam(6)*tau_l) - 1.0)/tau_c));
% Feedback co-efficients
a_f    = -11.034E-5; % U233 (drho/°C) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
a_g    = -05.814E-5; % U233 (drho/°C) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101

%% CORE HEAT TRANSFER PARAMETERS
% FUEL PARAMETERS - DONE
vdot_f   = 7.5708E-02; % ORNL-TM-0728 % 7.571e-2; % vol. flow rate (m^3/s) ORNL-TM-1647 p.3, ORNL-TM-0728 p.12
rho_f    = 2.1465E+03; % (partially enriched U-235)ORNL-TM-0728 p.8 2.243E+03; % (Th-U) density of fuel salt (kg/m^3) ORNL-TM-0728 p.8
W_f      = 1.623879934566580e+02; %vdot_f*rho_f; % 182.78; % calcd from m_dot*cp*delT=P; vdot_f*rho_f; % fuel flow rate (kg/s)
% tau_f_c  = tau_c; % ORNL-TM-0728 % 8.45; % transit time of fuel in core (s) ORNL-TM-1070 p.15, TDAMSRE p.5
m_f      = W_f*tau_c; % fuel mass in core (kg)
nn_f     = 2; % number of fuel nodes in core model
mn_f     = m_f/nn_f; % fuel mass per node (kg)
cp_f     = 4.19*9/5; % (MJ/deg-C) total fuel heat capacity TDAMSRE p.5
scp_f    = 1.9665E-3;%cp_f/m_f;% specific heat capacity of fuel salt (MJ/kg-C) ORNL-TM-0728 p.8

m_fcp = cp_f/scp_f; 

% Core Upflow - DONE
v_g      = 1.954; % graphite volume(m^3) ORNL-TM-0728 p. 101
rho_g    = 1.860E3; % graphite density (kg/m^3) ORNL-3812 p.77, ORNL-TM-0728 p.87
m_g      = v_g*rho_g; % graphite mass (kg)
cp_g     = 3.6*9/5; % TDAMSRE p.5 graphite total heat capacity (MW-s/C) ORNL-TM-1647 p.3
scp_g    = 1.773E-3; % http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/c.html cp_g/m_g; % graphite specific heat capacity (MW-s/kg-C) ORNL-TM-1647 p.3
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


%% multiregion core 

tau_fa   = [1.386, 2.083, 1.139, 1.424, 2.084, 1.139, 1.424, 2.371, 1.610]; % tau for fuel lump 'a' in each region 'i' in sec. ORNL-TM-1070 pg 67
tau_fb   = [1.454, 1.424, 1.139, 2.772, 1.424, 1.139, 2.774, 1.380, 2.700]; % tau for fuel lump 'b' in each region 'i' in sec. ORNL-TM-1070 pg 67


ka       = [0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333]; % power generation fraction for fuel lumps 'a' ORNL-TM-1070 pg 67
kb       = [0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311]; % power generation fraction for fuel lumps 'b' ORNL-TM-1070 pg 67
kga      = [0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277]; % heat transferred from for graphite lump to fuel lump 'a' ORNL-TM-1070 pg 67
kgb      = [0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305]; % heat transferred from for graphite lump to fuel lump 'b' ORNL-TM-1070 pg 67

mcp_fa   = [0.0151, 0.0512, 0.0280, 0.0350, 0.0866, 0.0473, 0.0592, 0.2380, 0.1615] * (9/5);
m_fa     = mcp_fa/scp_f;
mcp_fb   = [0.0158, 0.0349, 0.0280, 0.0682, 0.0592, 0.0473, 0.1152, 0.1384, 0.2710] * (9/5);
m_fb     = mcp_fb/scp_f;
mcp_g    = [0.0700, 0.2114, 0.1606, 0.2056, 0.3576, 0.2718, 0.3478, 0.9612, 0.9421] * (9/5);
m_g      = mcp_g/scp_g;
hA       = [0.000392, 0.001204, 0.000900, 0.001174, 0.001977, 0.001525, 0.001985, 0.005445, 0.005360] * (9/5); 
Ifa      = [0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936]; 
Ifb      = [0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473];
Ig       = [0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343];

m_f1 = (mcp_fa(1) + mcp_fb(1))/scp_f;
m_f2 = (sum(mcp_fa([2 3 4])) + sum(mcp_fb([2 3 4])))/scp_f;
m_f3 = (sum(mcp_fa([5 6 7])) + sum(mcp_fb([5 6 7])))/scp_f;
m_f4 = (sum(mcp_fa([8 9])) + sum(mcp_fb([8 9])))/scp_f;

W_f1 = (mcp_fa(1)/scp_f)/tau_fa(1); % 
W_f2 = (mcp_fa(2)/scp_f)/tau_fa(2); % 
W_f3 = (mcp_fa(5)/scp_f)/tau_fa(5); % 
W_f4 = (mcp_fa(8)/scp_f)/tau_fa(8); % 

tau_1 = m_f1/W_f1;
tau_2 = m_f2/W_f2;
tau_3 = m_f3/W_f3;
tau_4 = m_f4/W_f4;
tau_m = 2;

v_fm = 0.2793269;
m_fm = W_f * tau_m;
m_fm1 = m_fm*(m_f1/m_f);
m_fm2 = m_fm*(m_f2/m_f);
m_fm3 = m_fm*(m_f3/m_f);
m_fm4 = m_fm*(m_f4/m_f);

tau_m1 = ((2*W_f)-(m_fm-m_fm1))/W_f1;
tau_m2 = ((2*W_f)-(m_fm-m_fm2))/W_f2;
tau_m3 = ((2*W_f)-(m_fm-m_fm3))/W_f3;
tau_m4 = ((2*W_f)-(m_fm-m_fm4))/W_f4;


% Initial conditions - DONE
% fuel nodes
Tf_in  = 6.3222E+02; % in °C ORNL-TM-1647 p.2
Tf_out = 6.6353E+02; % 6.5727E+02; % 6.6355E+02; % in °C ORNL-TM-1647 p.2
T0_fm  = Tf_out;  % in °C 

P1 = P*(Ifa(1)+Ifb(1))/(sum(Ifa)+sum(Ifb));
P2 = P*(sum(Ifa([2 3 4]))+sum(Ifb([2 3 4])))/(sum(Ifa)+sum(Ifb));
P3 = P*(sum(Ifa([5 6 7]))+sum(Ifb([5 6 7])))/(sum(Ifa)+sum(Ifb));
P4 = P*(sum(Ifa([8 9]))+sum(Ifb([8 9])))/(sum(Ifa)+sum(Ifb));

Tf_out1 = Tf_in + (P1/(W_f1*scp_f)); 
Tf_out2 = Tf_in + (P2/(W_f2*scp_f));
Tf_out3 = Tf_in + (P3/(W_f3*scp_f));
Tf_out4 = Tf_in + (P4/(W_f4*scp_f));

T0_f1a = Tf_in + ((Ifa(1)/(Ifa(1)+Ifb(1)))*(Tf_out1-Tf_in)); % in °C 
T0_f1b = T0_f1a + ((Ifb(1)/(Ifa(1)+Ifb(1)))*(Tf_out1-Tf_in)); % in °C 

T0_f2a = Tf_in + ((Ifa(2)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 
T0_f2b = T0_f2a + ((Ifb(2)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 
T0_f3a = T0_f2b + ((Ifa(3)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 
T0_f3b = T0_f3a + ((Ifb(3)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 
T0_f4a = T0_f3b + ((Ifa(4)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 
T0_f4b = T0_f4a + ((Ifb(4)/(Ifa(2)+Ifb(2)+Ifa(3)+Ifb(3)+Ifa(4)+Ifb(4)))*(Tf_out2-Tf_in)); % in °C 

T0_f5a = Tf_in + ((Ifa(5)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C  
T0_f5b = T0_f5a + ((Ifb(5)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C 
T0_f6a = T0_f5b + ((Ifa(6)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C 
T0_f6b = T0_f6a + ((Ifb(6)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C  
T0_f7a = T0_f6b + ((Ifa(7)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C 
T0_f7b = T0_f7a + ((Ifb(7)/(Ifa(5)+Ifb(5)+Ifa(6)+Ifb(6)+Ifa(7)+Ifb(7)))*(Tf_out3-Tf_in)); % in °C 

T0_f8a = Tf_in + ((Ifa(8)/(Ifa(8)+Ifb(8)+Ifa(9)+Ifb(9)))*(Tf_out4-Tf_in)); % in °C
T0_f8b = T0_f8a + ((Ifb(8)/(Ifa(8)+Ifb(8)+Ifa(9)+Ifb(9)))*(Tf_out4-Tf_in)); % in °C
T0_f9a = T0_f8b + ((Ifa(9)/(Ifa(8)+Ifb(8)+Ifa(9)+Ifb(9)))*(Tf_out4-Tf_in)); % in °C
T0_f9b = T0_f9a + ((Ifb(9)/(Ifa(8)+Ifb(8)+Ifa(9)+Ifb(9)))*(Tf_out4-Tf_in)); % in °C


% T0_f1a = Tf_in + ((m_fa(1)/(m_fa(1)+m_fb(1)))*(T0_fm-Tf_in)); % in °C 
% T0_f1b = T0_f1a + ((m_fb(1)/(m_fa(1)+m_fb(1)))*(T0_fm-Tf_in)); % in °C 
% 
% T0_f2a = Tf_in + ((m_fa(2)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% T0_f2b = T0_f2a + ((m_fb(2)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% T0_f3a = T0_f2b + ((m_fa(3)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% T0_f3b = T0_f3a + ((m_fb(3)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% T0_f4a = T0_f3b + ((m_fa(4)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% T0_f4b = T0_f4a + ((m_fb(4)/(m_fa(2)+m_fb(2)+m_fa(3)+m_fb(3)+m_fa(4)+m_fb(4)))*(T0_fm-Tf_in)); % in °C 
% 
% T0_f5a = Tf_in + ((m_fa(5)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% T0_f5b = T0_f5a + ((m_fb(5)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% T0_f6a = T0_f5b + ((m_fa(6)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% T0_f6b = T0_f6a + ((m_fb(6)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% T0_f7a = T0_f6b + ((m_fa(7)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% T0_f7b = T0_f7a + ((m_fb(7)/(m_fa(5)+m_fb(5)+m_fa(6)+m_fb(6)+m_fa(7)+m_fb(7)))*(T0_fm-Tf_in)); % in °C 
% 
% T0_f8a = Tf_in + ((m_fa(8)/(m_fa(8)+m_fb(8)+m_fa(9)+m_fb(9)))*(T0_fm-Tf_in)); % in °C 
% T0_f8b = T0_f8a + ((m_fb(8)/(m_fa(8)+m_fb(8)+m_fa(9)+m_fb(9)))*(T0_fm-Tf_in)); % in °C 
% T0_f9a = T0_f8b + ((m_fa(9)/(m_fa(8)+m_fb(8)+m_fa(9)+m_fb(9)))*(T0_fm-Tf_in)); % in °C 
% T0_f9b = T0_f9a + ((m_fb(9)/(m_fa(8)+m_fb(8)+m_fa(9)+m_fb(9)))*(T0_fm-Tf_in)); % in °C 

% graphite nodes
T0_g1  = T0_f1a+((kga(1)+kgb(1))*P/hA(1)); % in °C 
T0_g2  = T0_f2a+((kga(2)+kgb(2))*P/hA(2)); % in °C 
T0_g3  = T0_f3a+((kga(3)+kgb(3))*P/hA(3)); % in °C 
T0_g4  = T0_f4a+((kga(4)+kgb(4))*P/hA(4)); % in °C 
T0_g5  = T0_f5a+((kga(5)+kgb(5))*P/hA(5)); % in °C 
T0_g6  = T0_f6a+((kga(6)+kgb(6))*P/hA(6)); % in °C 
T0_g7  = T0_f7a+((kga(7)+kgb(7))*P/hA(7)); % in °C 
T0_g8  = T0_f8a+((kga(8)+kgb(8))*P/hA(8)); % in °C 
T0_g9  = T0_f9a+((kga(9)+kgb(9))*P/hA(9)); % in °C 

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
scp_p = scp_f; % fuel heat capacity (MJ/(kg-C))

% SECONDARY FLOW PARAMETERS - DONE
vdot_s   = 5.362E-02; % ORNL-TM-0728 p. 164 % 5.236E-02; % coolant volume flow rate (m^3/s) ORNL-TM-1647 p.3
rho_s    = 1.922e3; % coolant salt density (kg/m^3) ORNL-TM-0728 p.8
W_s      = 1.005793369810108e+02; % calcd from mdot*cp*delT; vdot_s*rho_s; % coolant flow rate (kg/s) ORNL-TM-1647 p.3 

m_s      = v_cool*in_m*rho_s; % coolant mass in PHE (kg) 
nn_s     = 4; % number of coolant nodes in PHE
mn_s     = m_s/nn_s; % coolant mass per node (kg)
scp_s    = 2.39E-3;%cp_s/m_s; %2.219E-03; % specific heat capacity of coolant (MJ/(kg-C) ORNL-TM-0728 p.8

A_phe    = 2.359E+01; % effective area for heat transfer (primary and secondary, m^2) ORNL-TM-0728 p.164

ha_p     = 6.480E-01; % heat transfer*area coefficient from primary to tubes (MW/C) ORNL-TM-1647 p.3
ha_s     = 3.060E-01; % heat transfer*area coefficient from tubes to secondary (MW/C) ORNL-TM-1647 p.3

% Primary Side 
mcp_pn   = mn_p*scp_p; % (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
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
Tp_in    = Tf_out; % 6.628E+02; % in °C ORNL-TM-1647 p.2
T0_p4    = Tf_in; % in °C ORNL-TM-1647 p.2
T0_p1    = Tp_in - (Tp_in-T0_p4)/4; % in °C
T0_p2    = Tp_in - 2*(Tp_in-T0_p4)/4; % in °C
T0_p3    = Tp_in - 3*(Tp_in-T0_p4)/4; % in °C

% Secondary nodes
Ts_in    = 5.4611E+02; % in °C ORNL-TM-1647 p.2
T0_s4    = 5.7939E+02; % in °C ORNL-TM-1647 p.2
T0_s1    = Ts_in + (T0_s4-Ts_in)/nn_s; % in °C
T0_s2    = Ts_in + 2*(T0_s4-Ts_in)/nn_s; % in °C
T0_s3    = Ts_in + 3*(T0_s4-Ts_in)/nn_s; % in °C
% Tube nodes
T0_t1    = (T0_p1*hA_pn+T0_s3*hA_sn)/(hA_pn+hA_sn); % in °C
T0_t2    = (T0_p3*hA_pn+T0_s1*hA_sn)/(hA_pn+hA_sn); % in °C

%% Radiator Parameters - NOT DONE

% Initial conditions -  DONE
% Primary nodes
Trp_in   = T0_s4; % in °C ORNL-TM-1647 p.2
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
W_rs       = 8.923430894987993e+01; %vdot_rs*rho_rs; % air flow rate (kg/s) 

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

%% Secondary-side dynamics as implemented in ORNL-TM-1070

%Heat exchanger

he_n1   = 1.10;   % section length hA/WCp, dimensionless 
he_n2   = 1.366; 
he_ns   = 0.1363;
he_t1   = 2.01;   % transport time, sec
he_t2   = 2.29; 
he_tau1 = 0.569;  % heat transfer time constant mCp/hA, sec
he_tau2 = 0.304;
he_taus = 1.14;

Tp_in   = Tf_out;
T0_p2   = Tf_in;
T0_p1   = Tp_in - (Tp_in-T0_p2)/2;
Ts_in   = 5.4611E+02; %
T0_s2   = 5.7939E+02; % 
T0_s1   = Ts_in + (T0_s2 - Ts_in)/2;
T0_t    = ((he_tau1*he_tau2)/(he_tau1 + he_tau2))*(T0_p1/he_tau2 + T0_s1/he_tau1);
T0_sh   = T0_s1;

% Radiator
ra_n1   = 0.8803; % section length hA/WCp, dimensionless 
ra_n2   = 0.2591;
ra_ns   = 0;
ra_t1   = 6.52;   % transport time, sec
ra_t2   = 0.01;
ra_tau1 = 2.35;   % heat transfer time constant mCp/hA, sec
ra_tau2 = 19.7;

Trp_in  = T0_s4;
T0_rp2  = Ts_in;
T0_rp1  = Trp_in - (Trp_in-T0_rp2)/2;
Trs_in  = 37.78;
T0_rs2  = 148.9;
T0_rs1  = T0_rs2 - (T0_rs2-Trs_in)/2; 
T0_rt   = ((ra_tau1*ra_tau2)/(ra_tau1 + ra_tau2))*(T0_rp1/ra_tau1 + T0_rs1/ra_tau2);

%% Pure time delays between components - DONE

tau_hx_c = 8.67; % (sec) delay from hx to core TDAMSRE p.6
tau_c_hx = 3.77; % (sec) subtracted 1 sec for external loop power generation node resident time; delay from core to fuel hx TDAMSRE p.6
tau_hx_r = 4.71; % (sec) fertile hx to core TDAMSRE p.6
tau_r_hx = 8.24; % (sec) core to fertile hx TDAMSRE p.6