function D = tsh(A,B,the,bet,w,pm,den,lam)
D = (A*B)^(-1/the)*bsxfun(@times,bsxfun(@rdivide,w'.^bet.*pm'.^(1-bet),bsxfun(@times,pm,den)).^(-1/the),lam');
end