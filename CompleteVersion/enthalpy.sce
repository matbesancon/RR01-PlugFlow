function dH = enthalpy(T)
%finds the reaction enthalpy based on the temperature 
% -18.83 corresponds to the sum of mu(j)*cp(j) for the j components
% mu is the coefficient for j in the reaction 
%(negative if reactive, positive if product)
dH = -57500 - 18.83*(T-25)

end

