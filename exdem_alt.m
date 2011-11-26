function [w,L,lam,p_f] = exdem_alt(p)
%This function calculates excess demand for Lucas Alvarez model

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

out_err = 1;

w = ones(59,1);
L = p.siz/sum(p.siz);
lam = ones(59,1);

while out_err > 1e-6
    
    options = optimset('Display','iter','Jacobian','off','Algorithm','interior-point');
    %options = optimset('Display','iter');
    %lam = fmincon(@(x) fakeobj(x),lam,[],[],[],[],ones(59,1)*1e-8,ones(59,1)*100,@(x) lam_solve(x,A,B,bet,the,den,alp,relp,omeg,L),options);
    %lam = fsolve(@(x) lam_solve(x,A,B,bet,the,den,alp,relp,omeg,L),lam,options);
    lam = ktrlink(@(x) fakeobj(x),lam,[],[],[],[],ones(59,1)*1e-8,[],@(x) lam_solve(x,A,B,bet,the,den,alp,relp,omeg,L),options);
    [~,w,F,S,pm] = lam_solve(lam,A,B,bet,the,den,alp,relp,omeg,L);
    
    eps = w.*(1+(1-S).*(1-F)./(bet*F));
    new_L = (siz./eps)/sum(siz./eps);
    out_err = norm(new_L-L)
    L = max(new_L,1e-9);
    L = L/sum(L);
end

p_f = alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp.*pm;

end