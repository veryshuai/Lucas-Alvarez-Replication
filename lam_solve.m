function [diff,w,F,s,pm] = lam_solve(lam,A,B,bet,the,den,alp,relp,omeg,L)

w = ones(59,1);
maxits = 500;
it = 0;
pm = ones(size(w));
err = 1 ;
while err > 1e-12
            options = optimset('Display','off','Jacobian','on','MaxIter',50);

            %pm = ktrlink(@(x) fakeobj(x),pm,[],[],[],[],zeros(59,1),ones(59,1)*inf,@(x) int_p(x,A,B,log(w),lam,bet,the,den),options);
            pm = fsolve(@(x) int_p(x,A,B,log(w),lam,bet,the,den),pm,options);
            pm = exp(pm);

            D = tsh(A,B,the,bet,w,pm,den,lam);
            F = sum(D.*omeg,2);
            s = (alp*(1-(1-bet)*F))./((1-alp)*bet*F+alp*(1-(1-bet)*F));
            
            new_w = w.*(1+real_ex_dem(w,L,s,F,D,omeg)./L);
            
            err = norm(real_ex_dem(w,L,s,F,D,omeg));
            %err = norm(new_w - w);
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

end