clear; close all; clc;
run('Parameters_MSRE_cold_slug_U233.m')
i=1;
for temp = [482.22 537.79 593.33];
for volm = [1003.74 1007.48 1011.22];
   
% Cold Slug
% 482.22C = 900F % 537.79C = 1000F % 593.33C = 1100F
colddata = [0 temp 0];
% 3.74sec = 10ft^3 % 7.48sec = 20ft^3 % 11.22sec = 30ft^3
coldtime = [0 1000 volm];

cold = timeseries(colddata,coldtime);

sim('MSRE_final.slx');

indices       = find(tout>950);
new_tout{i}   = tout(indices);
new_power{i}  = msre_core_mux(indices,4);
% new_pc{i}     = msre_core_mux(indices,[5 7 9 11 13 15]);
% new_grph{i}   = msre_core_mux(indices,17);
new_fuel_1{i} = msre_core_mux(indices,18);
new_fuel_2{i} = msre_core_mux(indices,19);
i=i+1;
end
end
figure;semilogy(new_tout{1}-1000,new_power{1},new_tout{2}-1000,new_power{2},new_tout{3}-1000,new_power{3}); grid on;
figure;semilogy(new_tout{4}-1000,new_power{4},new_tout{5}-1000,new_power{5},new_tout{6}-1000,new_power{6}); grid on;
figure;semilogy(new_tout{7}-1000,new_power{7},new_tout{8}-1000,new_power{8},new_tout{9}-1000,new_power{9}); grid on;
figure;plot(new_tout{1}-1000,(new_fuel_1{1}+new_fuel_2{1})/2,new_tout{2}-1000,(new_fuel_1{2}+new_fuel_2{2})/2,new_tout{3}-1000,(new_fuel_1{3}+new_fuel_2{3})/2); grid on;
figure;plot(new_tout{4}-1000,(new_fuel_1{4}+new_fuel_2{4})/2,new_tout{5}-1000,(new_fuel_1{5}+new_fuel_2{5})/2,new_tout{6}-1000,(new_fuel_1{6}+new_fuel_2{6})/2); grid on;
figure;plot(new_tout{7}-1000,(new_fuel_1{7}+new_fuel_2{7})/2,new_tout{8}-1000,(new_fuel_1{8}+new_fuel_2{8})/2,new_tout{9}-1000,(new_fuel_1{9}+new_fuel_2{9})/2); grid on;
% figure;plot(new_tout{1}-1000,new_fuel_2{1},new_tout{2}-1000,new_fuel_2{2},new_tout{3}-1000,new_fuel_2{3},new_tout{4}-1000,new_fuel_2{4},new_tout{5}-1000,new_fuel_2{5}); grid on;
% figure;plot(new_tout{1}-1000,(new_fuel_1{1}+new_fuel_2{1})/2,new_tout{2}-1000,(new_fuel_1{2}+new_fuel_2{2})/2,new_tout{3}-1000,(new_fuel_1{3}+new_fuel_2{3})/2,new_tout{4}-1000,(new_fuel_1{4}+new_fuel_2{4})/2,new_tout{5}-1000,(new_fuel_1{5}+new_fuel_2{5})/2); grid on;
% figure;plot(new_tout{1}-1000,new_grph{1},new_tout{2}-1000,new_grph{2},new_tout{3}-1000,new_grph{3},new_tout{4}-1000,new_grph{4},new_tout{5}-1000,new_grph{5}); grid on;

figure(1);title('482\circC Slug n/n_0');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');
figure(2);title('538\circC Slug n/n_0');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');
figure(3);title('593\circC Slug n/n_0');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');
figure(4);title('482\circC Slug Average Core Temperature');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');
figure(5);title('538\circC Slug Average Core Temperature');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');
figure(6);title('593\circC Slug Average Core Temperature');xlabel('Time (s)');ylabel('Fractional Power'); legend ('0.28 m^3', '0.57 m^3', '0.85 m^3');

% figure(4);title('Core Outlet Temperature');xlabel('Time (s)');ylabel('Temperature (C)');
% figure(5);title('Core Average Temperature');xlabel('Time (s)');ylabel('Temperature (C)');
% figure(6);title('Graphite Temperature');xlabel('Time (s)');ylabel('Temperature (C)');

p = figure('color','w');
p(1) = subplot(3,2,1); semilogy(new_tout{1}-1000,new_power{1},new_tout{2}-1000,new_power{2},new_tout{3}-1000,new_power{3}); grid on;
p(2) = subplot(3,2,3); semilogy(new_tout{4}-1000,new_power{4},new_tout{5}-1000,new_power{5},new_tout{6}-1000,new_power{6}); grid on;
p(3) = subplot(3,2,5); semilogy(new_tout{7}-1000,new_power{7},new_tout{8}-1000,new_power{8},new_tout{9}-1000,new_power{9}); grid on;
p(4) = subplot(3,2,2); plot(new_tout{1}-1000,(new_fuel_1{1}+new_fuel_2{1})/2,new_tout{2}-1000,(new_fuel_1{2}+new_fuel_2{2})/2,new_tout{3}-1000,(new_fuel_1{3}+new_fuel_2{3})/2); grid on;
p(5) = subplot(3,2,4); plot(new_tout{4}-1000,(new_fuel_1{4}+new_fuel_2{4})/2,new_tout{5}-1000,(new_fuel_1{5}+new_fuel_2{5})/2,new_tout{6}-1000,(new_fuel_1{6}+new_fuel_2{6})/2); grid on;
p(6) = subplot(3,2,6); plot(new_tout{7}-1000,(new_fuel_1{7}+new_fuel_2{7})/2,new_tout{8}-1000,(new_fuel_1{8}+new_fuel_2{8})/2,new_tout{9}-1000,(new_fuel_1{9}+new_fuel_2{9})/2); grid on;

