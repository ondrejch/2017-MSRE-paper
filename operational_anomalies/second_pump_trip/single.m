clear; close all; clc;
a='Parameters_MSRE_sec_pump_U235.m';
b='Parameters_MSRE_sec_pump_U233.m';
i=1;
for parm_file = [1]

    if parm_file==1
        run(a)
    elseif parm_file==2
        run(b)
    else
        disp('error')
    end

sim('MSRE_final.slx')

indices       = find(tout>900);
new_tout{i}   = tout(indices);
new_inlet{i}  = msre_core_mux(indices,1);
new_power{i}  = msre_core_mux(indices,4);
new_pc{i}     = msre_core_mux(indices,[5 7 9 11 13 15]);
new_grph{i}   = msre_core_mux(indices,17);
new_fuel_1{i} = msre_core_mux(indices,18);
new_fuel_2{i} = msre_core_mux(indices,19);
new_hx_out{i} = msre_he_mux(indices,12);
new_hx_in{i}  = msre_he_mux(indices,2);
i=i+1;
end

figure;plot(new_tout{1}-1000,new_power{1}); grid on;
figure;plot(new_tout{1}-1000,new_fuel_2{1}); grid on;
figure;plot(new_tout{1}-1000,(new_fuel_1{1}+new_fuel_2{1})/2); grid on;
figure;plot(new_tout{1}-1000,new_grph{1}); grid on;
figure;plot(new_tout{1}-1000,new_hx_out{1}); grid on;
figure;plot(new_tout{1}-1000,new_hx_in{1}); grid on;
figure;plot(new_tout{1}-1000,new_inlet{1}); grid on;

figure(1);title('n/n_0');xlabel('Time (s)');ylabel('Fractional Power');
figure(2);title('Core Outlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');
figure(3);title('Core Average Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');
figure(4);title('Graphite Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');
figure(5);title('Heat Exchanger Outlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');
figure(6);title('Heat Exchanger Inlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');
figure(7);title('Core Inlet Temperature');xlabel('Time (s)');ylabel('Temperature (\circC)');