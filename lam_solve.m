function  [lam,w,F,s,pm] =lam_solve(A,B,bet,the,den,alp,relp,omeg,L);

%function [diff,w,F,s,pm] = lam_solve(lam,A,B,bet,the,den,alp,relp,omeg,L)

w = ones(59,1);
lam = ones(59,1);
err = 1;
err_mid = 1;
maxits = 500;
it = 0;
while err_mid > 1e-12
while err > 1e-12
            options = optimset('Display','off','Jacobian','on','MaxIter',50);

            pm = fsolve(@(p) int_p(p,A,B,log(w),lam,bet,the,den),ones(size(w)),options);
            pm = exp(pm);

            D = tsh(A,B,the,bet,w,pm,den,lam);
            F = sum(D.*omeg,2);
            s = (alp*(1-(1-bet)*F))./((1-alp)*bet*F+alp*(1-(1-bet)*F));
            
            new_w = w.*(1+real_ex_dem(w,L,s,F,D,omeg)./L);
            
            err = norm(new_w - w);
            w = max(new_w,1e-12);
            it = it + 1;
            if it > maxits
                display('did not converge')
                w
                L
                lam
                %den
                break;
            end
end 

diff = relp - alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp; 

C = alp^(-alp)*(1-alp)^(-1+alp)*w.^(alp*(1-bet)).*pm.^(alp*(bet-1))./(A*B*diag(den));

new_lam = max(lam.^(the*alp).*(1+.1*(diff./C)),1e-12);

err_mid = norm(new_lam - lam);
lam = new_lam;
end