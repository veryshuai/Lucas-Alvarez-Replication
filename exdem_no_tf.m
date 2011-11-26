function [out,p_f] = exdem_no_tf(w_orig,L,lam,p)
%This function calculates excess demand for Lucas Alvarez model

alp = p.alp;
bet = p.bet;
the = p.the;
eta = p.eta;
dist = p.dist;

%make omeg easier to work with
omeg = ones(size(w_orig,1));

A = gamma(1+the*(1-eta));
B = bet^-bet * (1-bet)^(-1+bet);

den = dist.*omeg;

err = 1;
w = w_orig;
v = 1;

while err>1e-6

options = optimset('Display','off','Jacobian','on');

pm = fsolve(@(p) int_p(p,A,B,log(w),lam,bet,the,den),ones(size(w)),options);
pm = exp(pm);

D = tsh(A,B,the,bet,w,pm,den,lam);

F = sum(D.*omeg,2);

s = (alp*(1-(1-bet)*F))./((1-alp)*bet*F+alp*(1-(1-bet)*F));

w_new = w.*(1+v*real_ex_dem(w,L,s,F,D,omeg)./L);

err = norm(w_new - w);

w = w_new;
end

out = w;

p_f = alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp.*pm;

end