function [out,jac] = int_p(p,A,B,w,lam,bet,the,den)
%function [empty1,out,jac] = int_p(p,A,B,w,lam,bet,the,den)
%empty1 = [];
%out = A*B*(sum(bsxfun(@times,bsxfun(@rdivide,w'.^bet.*p'.^(1-bet),den).^(-1/the),lam'),2).^(-the)))-p;
out = log(A*B) - the*log(sum(bsxfun(@times,den.^(1/the),exp(-(1/the)*((1-bet)*p'+bet*w')).*lam'),2))-p;
jac_den = sum(bsxfun(@times,den.^(1/the),(exp(p)'.^(1-bet).*exp(w)'.^bet).^(-1/the).*lam'),2);
ksi = bsxfun(@rdivide,bsxfun(@times,den.^(1/the),(exp(p)'.^(1-bet).*exp(w)'.^bet).^(-1/the).*lam'),jac_den);
jac = (1-bet)*ksi-eye(size(p,1));
end