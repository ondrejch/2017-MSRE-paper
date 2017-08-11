clear; close all; clc;
run('Parameters_MSRE_pump_trip_U235.m')
i=1;
for therm_circ = [.01 .02 .03 .04 .05];
sim('MSRE_pump.slx')

indices        = find(tout>980);
new_tout{i}    = tout(indices);
new_fuel_hx{i} = msre_he_mux(indices,6);
new_power{i}   = msre_core_mux(indices,4);
new_pc{i}      = msre_core_mux(indices,[5 7 9 11 13 15]);
new_grph{i}    = msre_core_mux(indices,17);
new_fuel_1{i}  = msre_core_mux(indices,18);
new_fuel_2{i}  = msre_core_mux(indices,19);
new_hx_out{i}  = msre_he_mux(indices,12);
new_hx_in{i}   = msre_he_mux(indices,2);
i=i+1;
end
figure;plot(new_tout{1}-1000,new_power{1},new_tout{2}-1000,new_power{2},new_tout{3}-1000,new_power{3},new_tout{4}-1000,new_power{4},new_tout{5}-1000,new_power{5}); grid on;
figure;plot(new_tout{1}-1000,new_fuel_2{1},new_tout{2}-1000,new_fuel_2{2},new_tout{3}-1000,new_fuel_2{3},new_tout{4}-1000,new_fuel_2{4},new_tout{5}-1000,new_fuel_2{5}); grid on;
figure;plot(new_tout{1}-1000,(new_fuel_1{1}+new_fuel_2{1})/2,new_tout{2}-1000,(new_fuel_1{2}+new_fuel_2{2})/2,new_tout{3}-1000,(new_fuel_1{3}+new_fuel_2{3})/2,new_tout{4}-1000,(new_fuel_1{4}+new_fuel_2{4})/2,new_tout{5}-1000,(new_fuel_1{5}+new_fuel_2{5})/2); grid on;
figure;plot(new_tout{1}-1000,new_grph{1},new_tout{2}-1000,new_grph{2},new_tout{3}-1000,new_grph{3},new_tout{4}-1000,new_grph{4},new_tout{5}-1000,new_grph{5}); grid on;
figure;plot(new_tout{1}-1000,new_hx_out{1},new_tout{2}-1000,new_hx_out{2},new_tout{3}-1000,new_hx_out{3},new_tout{4}-1000,new_hx_out{4},new_tout{5}-1000,new_hx_out{5}); grid on;
figure;plot(new_tout{1}-1000,new_hx_in{1},new_tout{2}-1000,new_hx_in{2},new_tout{3}-1000,new_hx_in{3},new_tout{4}-1000,new_hx_in{4},new_tout{5}-1000,new_hx_in{5}); grid on;
figure;plot(new_tout{1}-1000,new_fuel_hx{1},new_tout{2}-1000,new_fuel_hx{2},new_tout{3}-1000,new_fuel_hx{3},new_tout{4}-1000,new_fuel_hx{4},new_tout{5}-1000,new_fuel_hx{5}); grid on;


figure(1);title('n/n_0');xlabel('Time (s)');ylabel('Fractional Power'); legend ('1%','2%','3%','4%','5%');
figure(2);title('Core Outlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
figure(3);title('Core Average Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
figure(4);title('Graphite Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
figure(5);title('Heat Exchanger Outlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
figure(6);title('Heat Exchanger Inlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
figure(7);title('Fuel salt leaving HE');xlabel('Time (s)');ylabel('Temperature (\circC)'); legend ('1%','2%','3%','4%','5%');
