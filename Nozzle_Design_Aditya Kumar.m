clc;
clear;
%% Inputs
F = input('Enter the Thrust you want generated(N): ');
P_st = input('Enter the stagnation pressure(psi): ');
P_e = input('Enter the exit pressure(psi): ');
M_th = input('Enter throat Mach: ');
T_st = input('Enter stagnation temperature(K): ');

%% Assumptions
Y = 1.4; % Constant pressure and volume ratio 
R = 287; %Universal Gas Constant
P_a = 14.7; %Ambient Pressure (psi)
T_in = T_st; %Stagnant temperature
P_in = P_st; %Stagnant temperature

%% Conversion
P_conv = 6.89476; %Converting psi to kpa
A_conv = 1550; %Converting m^2 to in^2
P_a_kpa = P_a * P_conv; %Converting Ambient Pressure from psi to KPa
P_e_kpa = P_e * P_conv; %Converting Exit Pressure from psi to KPa
P_in_kpa = P_in * P_conv; %Converting Inlet Pressure from psi to KPa
%% Calculations
M_e = sqrt((2/(Y-1)) * (((P_in_kpa/P_e_kpa)^((Y-1)/(Y))) - 1)); %Calculating exit Mach speed
A_ratio = (M_th/M_e)*sqrt(((1+(((Y-1)/2)*(M_e^2)))/(1+(((Y-1)/2)*(M_th^2))))^(((Y+1)/(Y-1)))); %Area ratio
T_e = T_st/(1+((Y/1)/2)*M_e^2); %Calculation for exit temperature
a = sqrt(Y*R*T_e); %Ca
U_e = M_e*a;

if P_e == P_a
    m_flow = F/U_e;
    A_th = (m_flow*sqrt(Y*R*T_in)/(1000*P_in_kpa*Y*(2/(Y+1))^((Y+1)/(2*(Y-1))))) * A_conv;
    A_e = A_th*A_ratio;
    r_e = sqrt(A_e/pi);
    r_th = sqrt(A_th/pi);

elseif P_e < P_a
    A_th = (F / ((U_e*Y*((2/(Y+1))^((Y+1)/(Y-1)))) + (A_ratio * (P_e_kpa - P_a_kpa)))) * A_conv;
    A_e = A_th*A_ratio;
    r_e = sqrt(A_e/pi);
    r_th = sqrt(A_th/pi);
    
elseif P_e > P_a
    A_th = (F / ((U_e*Y*((2/(Y+1))^((Y+1)/(Y-1)))) + (A_ratio * (P_e_kpa - P_a_kpa)))) * A_conv;
    A_e = A_th*A_ratio;
    r_e = sqrt(A_e/pi);
    r_th = sqrt(A_th/pi);
end

fprintf('\nThe Throat Diameter is (in):');
disp(2*r_th)
fprintf('The Exit Diameter is (in):');
disp(2*r_e)

L_cone = (r_e - r_th)/(tand(15));
Ax_slope = (r_e - r_th)/(L_cone);
x = -1.5*r_th:0.0001:L_cone;
Ax_At = r_th + (Ax_slope)*(x/L_cone).^2;
Ax_ve = -Ax_At;

plot(x,Ax_At)
title('2D Sketch of required Nozzle')
grid on
hold on
plot(x,Ax_ve)
xlabel('Length (in)')
ylabel('Diameter (in)')