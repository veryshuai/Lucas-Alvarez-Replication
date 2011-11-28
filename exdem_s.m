function [w,L,lam,p_f,pm,eps] = exdem_alt_alt_alt(p)
%This function calculates excess demand for Lucas Alvarez model

load results.mat 
clearvars -except p lam

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
%lam = siz.*relp.^(bet/(alp*bet));
%lam = max(lam,.1);
den = dist.*omeg;
    
options = optimset('Display','iter','MaxIter',2e3,'MaxFunEvals',2e3);
lam = fminsearch(@(x) lam_solve_s(x,A,B,bet,the,den,alp,relp,omeg,siz),log(lam),options);
%lam = fsolve(@(x) lam_solve_s(x,A,B,bet,the,den,alp,relp,omeg,siz),log(lam),options);
lam = exp(lam);
[~,w,L,pm,eps] = lam_solve_s(log(lam),A,B,bet,the,den,alp,relp,omeg,siz);

p_f = alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp.*pm;

end