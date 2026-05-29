function err = impact_error(r,Xe,theta,z_10_ret,G_10,W,idx,f_dyn,cmax,n,options)

u = [cos(theta);sin(theta)];

X_nom = Xe';

X_cible = X_nom(idx) +  r*u;

x_10_ret = z_10_ret(1:n);

c_10_r = inv(G_10(idx,:)'*G_10(idx,:)+W)*G_10(idx,:)'*r*u; % delta_cible = G*c
c_10_sat = c_10_r;

if norm(c_10_sat) > cmax
    c_10_sat = cmax * c_10_sat/norm(c_10_sat);
end 

delta_x10_r = zeros(n,1);

delta_x10_r(4:6) = c_10_sat;

state0_10_r = x_10_ret + delta_x10_r;

[~,~,~,Xe_10_r] = ode45(@(t,x) f_dyn(x), [10 200], state0_10_r,options); 

Xe_10 = Xe_10_r';

err = norm(X_cible - Xe_10(idx))/norm(X_cible-X_nom(idx));

end