function Z = real_ex_dem(w,L,s,F,D,omeg)

Z = 1./w.*(sum(bsxfun(@times,L'.*w'.*(1-s')./F',D'.*omeg'),2)-L.*w.*(1-s));

end