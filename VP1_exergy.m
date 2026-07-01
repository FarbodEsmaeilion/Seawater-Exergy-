*** Source: https://doi.org/10.1016/j.enconman.2016.04.106 
*** Dr. Farbod Esmaeilion

function [ex]=VP1_exergy(x,y)
T=(x+y)/2;                                                                 %K
cp=2.125e3-11.017*T+0.049862*T^2-7.7663e-5*T^3+4.394e-8*T^4;               %J/kg K
h_x=-18.34+1.498*(x-273.15)+.001377*(x-273.15)^2;                          %kJ/kg
h_y=-18.34+1.498*(y-273.15)+.001377*(y-273.15)^2;                          %kJ/kg
ex=(h_x-h_y)-298.15*cp*log(x/y)/1000;                                      %kJ/kg
end