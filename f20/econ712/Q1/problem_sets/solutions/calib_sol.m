function [ D ] = calib_sol( c_in, alpha, beta, delta, z, k_bar, k_bar2, c_bar2)
% this function gives the distance to the new steady state after 30
% periods given the initial capital k_bar and consumption c_star
T = 50;
Trj = zeros(2,T);
Trj(1,1) = k_bar;
Trj(2,1) = c_in;
for i = 2:T
    k = Trj(1,i-1);
    c = Trj(2,i-1);
    Trj(1,i) = z*k.^alpha + (1-delta)*k - c;
    Trj(2,i) = (1-delta+z*alpha*Trj(1,i).^(alpha-1))*beta*c;
end
D = norm(Trj(:,T)-[k_bar2; c_bar2]);

