// Adiabatic
// Creates a model of an adiabatic plug-flow reactor
// Reaction : A + E --> R
// given the final conversion rate, and the thermodynamic and kinetic data

// Author : Mathieu Besançon
// UTC student
// Process Engineering

// Execute the enthalpy.sce first, which contains the enthalpy function

Xa = [0:0.001:0.8435914];
//Xa conversion rate, from 0 to the maximal conversion rate (given by the
//system at the end / equilibrium
T = ones(1,844); // T vector represents the temperature corresponding to each conversion rate
CpE=75.33;
CpA=189.7;
CpR=123.1;
// Cpx calorific capacity for the x component
Ca0=0.4;
Ce0=50;
// Initial concentrations in A and E component
Ea=46000;
k0=2160;
R=8.314;
//Kinetic constants of the reaction
T(1) = 50; //Temperature of the inflow
Qo=10^-3; // Initial volumic flow (m^3/s)
ratio = 200; // Ratio = L/D characterizes the geometry of the reactor, supposed cylindric

//-------------------------------------
// Generation of the T vector from the reaction enthalpy
for i=2:844
    T(i)= (-1)*enthalpy(T(i-1))*(Xa(i)-Xa(i-1))/(CpE*(Ce0-Ca0*Xa(i))+Ca0*CpA*(1-Xa(i))+2*CpR*Ca0*Xa(i)) + T(i-1);
end

//-------------------------------------
// Computation of the kinetic constant k = k0*exp(-Ea/RT)*Ce
// Ca and Ce depend on the conversion rate Xa
// Reaction is pseudo first-order in Ca if Ce remains constant
k=ones(1,844);
for i=1:844
    k(i)= k0*(Ce0-Xa(i)*Ca0)*exp(-Ea/(R*(T(i)+273.15)));
end

//-------------------------------------
// Computation of the time vector from the differential equation :
// dCa/dt = k*Ca
t=ones(1,844);
t(1)=0;
for i = 2:844
    t(i)=t(i-1)+(Xa(i)-Xa(i-1))/((1-Xa(i))*k(i));
end

//-------------------------------------
//Determination of the volume and geometry of the reactor
V=0;
rho = zeros(1,844);
for i=1:844
    rho(i)= Ca0*(1-Xa(i))*102.9 + Ca0*Xa(i)*0.5*60 + (Ce0-Xa(i)*Ca0)*18;
end
    for i=2:844
    V=V+Qo*rho(1)*(t(i)-t(i-1))*2/(rho(i)+rho(i-1))
end
L = (4*V*(ratio^2)/%pi)^(1/3);
D = L/ratio;
S = %pi*0.25*D^2;


//-------------------------------------
// Computation of the x vector, position for each temperature and conversion
// From the differential equation dx/dt = v = Qo/S = Qo/(pi*D^2/4))
x = zeros(1,844);
for i=2:844
    x(i)=x(i-1)+L/843;
end

//-------------------------------------
// Ploting the two functions T(x) and Xa(x)
scf();
plot(x,Xa,'Color','blue');
xgrid(1);
title('Visualization of (x,Xa)');
xlabel('x(m)');
ylabel('Xa');
scf();
plot(x,T,'Color','red');
xgrid(1);
title('(x,T)');
xlabel('x(m)');
ylabel('T(°C)');


//-------------------------------------
// We have to check the plug-flow reactor hypothesis
// Plug-flow : high Reynolds Number

// Computation of the volumic mass = sum of massic concentration (x)
// With massic concentration = M*molar concentration

//We observe a low variation of rho. However, for the accuracy of the results, 
//we compute Re as a vector too

//mu is the dynamic viscosity of water, we use a correlation 
//in order to compute mu(T)
Re = zeros(1,844);
mu = Re;

for i =1:844
mu(i)= 2.414*10^-5 *10^(247.8/(T(i)-133.15))
Re(i) = Qo*rho(i)*D/(mu(i)*S);
end

scf();
plot(x,Re,'Color','b');
title('Reynold''s number in the reactor');
xlabel('x (m)');
ylabel('Re');

//The Reynolds number is high, therefore the plug-flow hypothesis is true
//Rho variation is negligible


    
