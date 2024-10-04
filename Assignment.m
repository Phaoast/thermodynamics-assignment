clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder);
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set.
%% values should not be changed. These are used by all Nasa Functions.
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions.
%  For the final assignment take the ones from the specific case you are supposed to do.
v1=200;Tamb=250;P3overP2=8;Pamb=55000*kPa;mfurate=0.68*kg/s;AF=102.78;
cFuel='CH4';
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Here starts your part (compressor,combustor,turbine and nozzle)...

%% Compressor solution using interpolation
cMethod = 'Interpolation Method';
sPart = 'Compressor';
P3 = P2 * P3overP2;
S3 = S2;
% We find the temperature for which the thermal part of entropy is s2thermal + Rg*log(P3/P2)
% which was found using S_3 - S_2 = 0
T3 = interp1(sair_a, TR, s2thermal + Rg*log(P3/P2));
h3 = interp1(TR, hair_a, T3);
v3 = 0;

% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2, T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa, P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2, v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ, h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ, S3/kJ);

%% Combustor solution using interpolation
cMethod = 'Interpolation Method';
sPart = 'Combustor';
P4 = P3;

%{
Assumptions (terms neglected): \dot{Q_{cv}}, \dot{W_{cv}}, kinetic energy, potential energy

Conservation of energy:
\dot{m_air} h_3 + \dot{m_{fuel}} h_fuel = \dot{m_products} h_4      , where m_3 is the mass flow rate of air, m_{fuel} is the mass flow rate of fuel, m_4 is the mass flow rate of the products of combustion
%}
v4 = 0;
mAir = mfurate * AF;
mFuel = mfurate;
mProducts = mAir + mFuel;
hFuel = interp1(TR, Yfuel * hia', Tref);
h4 = (mAir * h3 + mFuel * hFuel) / mProducts;
% We use the data from the database for the fuel
x = Sp(iSp(1)).Elcomp(3);
y = Sp(iSp(1)).Elcomp(2);
MFuel = Sp(iSp(1)).Mass;
% We find the composition of the products of combustion using the method in Turns
a = x + y / 4;
AFStoic = (4.76 * a) * (MAir / MFuel);
equivalenceRatio = AFStoic / AF;
b = x;
c = y / 2;
d = (2 * (a / equivalenceRatio) - 2 * b - c) / 2;
e = (a * 3.76) / equivalenceRatio;
NTot = b + c + d + e;
XCO2 = b / NTot;
XH2O = c / NTot;
XO2 = d / NTot;
XN2 = e / NTot;
XProducts = [0 XO2 XCO2 XH2O XN2];
% We find the enthalpy of the products of combustion
MProducts = XProducts * Mi';
YProducts = XProducts .* Mi / MProducts;
hProducts_a = YProducts * hia';
% We use the enthalpy we found to find the temperature
T4 = interp1(hProducts_a, TR, h4);

sProducts_a = YProducts * sia';
s4thermal = interp1(TR, sProducts_a, T4);
RgProducts = Runiv / MProducts;
S4 = s4thermal - RgProducts * log(P4/Pref);

% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,3,4);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T3, T4);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P3/kPa, P4/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v3, v4);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h3/kJ, h4/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S3/kJ, S4/kJ);

%% Turbine solution using interpolation
h5 = h4 + h2 - h3;
v5 = 0;
S5 = S4;
T5 = interp1(hProducts_a, TR, h5);

% For P5, we use S = s_thermal - R_g*ln(P/P_ref)
s5thermal = interp1(TR, sProducts_a, T5);
lnP5 = ((s5thermal - S5)/(RgProducts)) + log(Pref);
P5 = exp(lnP5);

% Print to screen
cMethod = 'Interpolation Method';
sPart = 'Turbine';
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,4,5);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T4, T5);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P4/kPa, P5/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v4, v5);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h4/kJ, h5/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S4/kJ, S5/kJ);

%% Nozzle solution using interpolation
cMethod = 'Interpolation Method';
sPart = 'Nozzle';
P6 = Pamb;
S6 = S5;
% We use S = s_thermal - R_g*ln(P/P_ref), find s_thermal, and use it to find T6
s6thermal = S6 + RgProducts * log(P6/Pref);
T6 = interp1(sProducts_a, TR, s6thermal);
h6 = interp1(TR, hProducts_a, T6);
% Then, we use Turns' v6 = sqrt(2*(h5-h6))
v6 = sqrt(2 * (h5 - h6));

% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,5,6);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T5, T6);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P5/kPa, P6/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v5, v6);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h5/kJ, h6/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S5/kJ, S6/kJ);
