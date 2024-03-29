function [empty,out,p_f] = exdem(x,p)
%This function calculates excess demand for Lucas Alvarez model

empty = [];

w = x(1:59);
L = x(59+1:59+58);
L = [1-sum(L);L];
lam = x(59+58+1:end);

alp = p.alp;
bet = p.bet;
the = p.the;
eta = p.eta;
dist = p.dist;
omeg = p.omeg;
relp = p.relp;
siz = p.siz;

%make omeg easier to work with
omeg = repmat(omeg,1,size(omeg,1));
omeg(1:size(omeg,1)+1:end) = 1;

%put in trade agreements to match AL
%EU: Finland,Sweden,Austria,Spain,Portugal,Greece,Denmark,UK,Ireland,Netherlands,Belgium,Italy,Germany,France
%NAFTA: Mexico,USA,Canada
%CEFTA: Poland,Hungary,Czech Rep, Slovakia, Romania
%MERCOSUR: Argentina,Brazil,Venezuela

EU = [30 19 13 10 32 31 23 5 39 14 18 6 3 4];
NAFTA = [11 1 9];
CEFTA = [27 48 46 57 51];
MERCOSUR = [16 8 36];

for k = EU
    for m = EU
        omeg(k,m) = 1;
    end
end

for k = NAFTA
    for m = NAFTA
        omeg(k,m) = 1;
    end
end

for k = CEFTA
    for m = CEFTA
        omeg(k,m) = 1;
    end
end

for k = MERCOSUR
    for m = MERCOSUR
        omeg(k,m) = 1;
    end
end

A = gamma(1+the*(1-eta));
B = bet^-bet * (1-bet)^(-1+bet);

den = dist.*omeg;

options = optimset('Display','off','Jacobian','on');

pm = fsolve(@(p) int_p(p,A,B,log(w),lam,bet,the,den),ones(size(w)),options);
pm = exp(pm);

D = tsh(A,B,the,bet,w,pm,den,lam);

F = sum(D.*omeg,2);

s = (alp*(1-(1-bet)*F))./((1-alp)*bet*F+alp*(1-(1-bet)*F));

Z = real_ex_dem(w,L,s,F,D,omeg);
Z = Z(1:58);

P = relp - alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp;

GDP = siz(1:59)-L(1:59).*w(1:59).*(1+(1-s(1:59)).*(1-F(1:59))./(bet*F(1:59)));

out = [Z;P;GDP];

p_f = alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp.*pm;

end