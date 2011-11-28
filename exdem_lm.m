function [w,L,p_f,eps] = exdem_lm(w_orig,L_orig,lam,p)
%This function calculates excess demand for Lucas Alvarez model

alp = p.alp;
bet = p.bet;
the = p.the;
eta = p.eta;
dist = p.dist;
omeg = p.omeg;

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

err_out = 1;
maxits = 500;
w = w_orig;
L = L_orig;
v = 1;
u = 1;
pm = ones(size(w))*2;
while err_out>1e-9
err = 1;
it = 0;
    while err > 1e-3
            
            options = optimset('Display','off','Jacobian','on','MaxIter',1000,'TolX',1e-10,'TolFun',1e-10);

            %pm = ktrlink(@(x) fakeobj(x),pm,[],[],[],[],zeros(59,1),ones(59,1)*inf,@(x) int_p(x,A,B,log(w),lam,bet,the,den),options);
            pm = fsolve(@(x) int_p(x,A,B,log(w),lam,bet,the,den),pm,options);
            pm = exp(pm);

            D = tsh(A,B,the,bet,w,pm,den,lam);
            F = sum(D.*omeg,2);
            s = (alp-alp*(1-bet)*F)./(alp-(alp-bet)*F);
            
            new_w = w.*(1+real_ex_dem(w,L,s,F,D,omeg)./L);
            
            err = norm(real_ex_dem(w,L,s,F,D,omeg));
            %err = norm(new_w - w);
            w = max(abs(new_w),1e-2);
            it = it + 1;
            if it > maxits
                %display('did not converge')
                %w
                %L
                %lam
                %den
                break;
            end
end
eps = w.*(1+((1-s).*(1-F))./(bet*F));
eps = abs(eps); %just in case
    
p_f = alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp.*pm;

L_new = L.*(1+u*(eps./p_f-mean(eps./p_f)));
L_new = L_new/sum(L_new);
err_out = norm(L_new-L);
L = L_new;
end

end