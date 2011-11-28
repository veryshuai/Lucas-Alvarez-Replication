function [diff,w,L,pm,eps] = lam_solve_s(lam,A,B,bet,the,den,alp,relp,omeg,siz)

lam = exp(lam);
maxits = 50;
its = 0;
pm = ones(size(relp));
L = siz.*relp.^(-1/alp);
w = (lam./L).^(the/(bet+the));
err_mid = 1;
while err_mid > 1e-3
err = 1 ;
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
    eps = real(eps);
    new_L = (siz./eps)/sum(siz./eps);
    err_mid = norm(L-new_L);
    L = max(new_L,1e-12);
    L = L/sum(L);
    its = its + 1;
    if its > maxits
                display('did not converge')
                %w
                %L
                %lam
                %den
                break;
    end
end
%display('finished one evaluation');
%diff =norm(relp - alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp);
diff =norm(relp - alp^(-alp)*(1-alp)^(-1+alp)*(w./pm).^alp);
end