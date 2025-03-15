function k = k_surfeffect(TR, theta, v, R)
% theta in degree
va = v/sind(theta);
z = (2*pi*R)/(va*TR);
k = ((va * TR)/(pi*R)) * besselj(1, z);
% fprintf("coef: %f\n", k);

end