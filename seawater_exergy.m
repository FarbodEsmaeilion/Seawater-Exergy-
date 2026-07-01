*** Source: Second law analysis of reverse osmosis desalination plants: An alternative design using pressure retarded osmosis
*** Dr. Farbod Esmaeilion

function [ex]=seawater_exergy(T,ws,P)
T0=25;                                                                     
ws0=0.035;                                                                 
a1=8.02e2;
a2=-2.001;
a3=1.677e-2;
a4=-3.06e-5;
a5=-1.613e-5;
dw=9.999e2+2.034e-2*T-6.162e-3*T^2+2.261e-5*T^3-4.657e-8*T^4;              %Kg/m3
dsw=dw+ws*(a1+a2*T+a3*T^2+a4*T^3+a5*ws*T^2);                               %Kg/m3   
b1=-2.348e4;
b2=3.152e5;
b3=2.803e6;
b4=-1.446e7;
b5=7.826e3;
b6=-4.417e1;
b7=2.139e-1;
b8=-1.991e4;
b9=2.778e4;
b10=9.728e1;
hw_star=141.355+4202.07*T0-0.535*T0^2+0.004*T0^3;                          %j/kg
hsw_star=hw_star-ws*(b1+b2*ws+b3*ws^2+b4*ws^3+b5*T0+b6*T0^2+b7*T0^3+b8*ws*T0+b9*ws^2*T0+b10*ws*T0^2); %j/kg
hw=141.355+4202.07*T-0.535*T^2+0.004*T^3;                                   %j/kg
hsw=hw-ws*(b1+b2*ws+b3*ws^2+b4*ws^3+b5*T+b6*T^2+b7*T^3+b8*ws*T+b9*ws^2*T+b10*ws*T^2); %j/kg @p0
hsw = hsw+(P-101325)/dsw;                                                  
c1=-4.231e2;
c2=1.463e4;
c3=-9.88e4;
c4=3.095e5;
c5=2.562e1;
c6=-1.443e-1;
c7=5.879e-4;
c8=-6.111e1;
c9=8.041e1;
c10=3.035e-1;
sw_star=0.1543+15.383*T0-2.996e-2*T0^2+8.193e-5*T0^3-1.37e-7*T0^4;         %j/kgK
ssw_star=sw_star-ws*(c1+c2*ws+c3*ws^2+c4*ws^3+c5*T0+c6*T0^2+c7*T0^3+c8*ws*T0+c9*ws^2*T0+c10*ws*T0^2); %j/kgK
sw=0.1543+15.383*T-2.996e-2*T^2+8.193e-5*T^3-1.37e-7*T^4;                  %j/kgK  
ssw=sw-ws*(c1+c2*ws+c3*ws^2+c4*ws^3+c5*T+c6*T^2+c7*T^3+c8*ws*T+c9*ws^2*T+c10*ws*T^2); %j/kgK @p0 
a=-1*(b1+2*b2*ws+3*b3*ws^2+4*b4*ws^3+b5*T0+b6*T0^2+b7*T0^3+2*b8*ws*T0+3*b9*ws^2*T0+2*b10*ws*T0^2);
b=-1*(c1+2*c2*ws+3*c3*ws^2+4*c4*ws^3+c5*T0+c6*T0^2+c7*T0^3+2*c8*ws*T0+3*c9*ws^2*T0+2*c10*ws*T0^2);
f=a-(T0+273.15)*b;
gsw=hsw_star-(T0+273.15)*ssw_star;
muw=gsw-ws*f;                                                              %j/kg
mus=gsw+(1-ws)*f;                                                          %j/kg
ex=(hsw-hsw_star)-298*(ssw-ssw_star)+ws*(mus-7.0814e+04)+(1-ws)*(muw+7.2785e+03);   %j/kg


