%Lucas Alvarez 2006 replication program

clear;
diary on;
tic

%input parameters (notation from paper)
p = struct;

p.alp = .75;
p.bet = .5;
p.the = .15; %three values in AL .1,.15,.25
p.eta = 2;

%import bilateral distance data
dist = csvread('dist_dat.txt');
dist(:,1) = dist(:,1)/mean(dist(:,1));
dist(:,1) = exp(-.05*(dist(:,1)-1)); %there is a mysterious kappa hanging around in the paper's notation.
dist_mat = reshape(dist(:,1),59,59);
dist_mat = [zeros(1,60);zeros(59,1),dist_mat];
dist_mat(2:end,1) = dist(1:59,3);
dist_mat(1,2:end) = dist(1:59,3);
dist_mat = sortrows(dist_mat,1);
dist_mat = sortrows(dist_mat',1);
dist_mat = dist_mat(2:end,2:end);
p.dist = dist_mat;

%data
siz = [27.99;15.69;7.35;5.25;4.86;4.3;3.85;2.86;2.29;2.06;1.93;1.4;1.34;...
    1.32;1.31;1.09;.94;.91;.85;.80;.72;.61;.59;.58;.52;.50;.49;.47;.47;.47;...
    .41;.40;.36;.32;.32;.31;.29;.29;.29;.26;.25;.25;.23;.20;.19;.18;.18;...
    .16;.15;.15;.14;.12;.11;.11;.08;.08;.07;.06;.06;.05]/100;
relp = [1.37;1.65;1.48;.7;1.55;1.28;1.27;.7;1;1.3;1.35;.72;.7;1.13;1.51;...
    .48;1.02;1.7;1.61;1.61;1.57;.62;.6;1.5;1.7;1.67;.64;.65;.7;.7;1.39;...
    1.08;.97;1.65;.48;.7;.71;.7;2.02;1.29;.22;.67;.82;.76;1.29;.72;.52;...
    .7;.52;.24;.68;.38;.5;.61;.28;.3;.37;.43;.4;.57];
omeg = (100-[5.4;5.48;5.86;11.49;5.86;5.86;5.86;18.58;13.73;5.1;5.86;...
    14.26;33.44;5.3;5.86;12.47;12.4;.68;5.86;5.86;5.86;12.2;9.88;5.86;...
    0;4.04;18;15.55;12.48;8.3;5.86;5.86;5.86;7.55;4.9;11.7;12.28;9.18;...
    .16;5.86;27.6;11.22;10.25;39.9;4.84;13.3;7.04;24.6;13.42;10.2;21.3;...
    16.03;30.77;24.06;15.15;12.63;1.15;6.88;31.17;7.67])/100;

%doesn't make sense to talk about distance with rest of the world
siz_NOROW = [siz(1:3);siz(5:end)]/sum([siz(1:3);siz(5:end)]);
relp_NOROW = [relp(1:3);relp(5:end)];
omeg_NOROW = [omeg(1:3);omeg(5:end)];

p.siz = siz_NOROW;
p.relp = relp_NOROW;
p.omeg = omeg_NOROW;

%[w,L,lam,pf] = exdem_alt_alt(p);
%[w,L,lam,pf] = exdem_alt_alt_alt(p);
%[w,L,lam,pf] = exdem_ga(p);
[w,L,lam,pf,pm,eps] = exdem_s(p);

% 
% pat = ones(58+2*59);
% pat(59+1:58+59,59+1:58+59)=0;
% 
% options = optimset('Display','iter','JacobPattern',pat,'Algorithm','active-set','MaxIter',60);
% %options = optimset('Display','iter');
% 
% init = [repmat(p.siz(2:59)/sum(p.siz),3,1);1;1];
% 
% sol = ktrlink(@(x) fakeobj(x),init,[],[],[],[],zeros(58+2*59,1),[ones(59,1)*inf;ones(58,1);ones(59,1)*inf],@(x) exdem(x,p),options);
% %sol = fsolve(@(x) exdem(x,p),ones(59*3,1),options);
% 
% w = [1;sol(1:59)]
% L = sol(59+1:59+58)/sum(sol(59+1:59+58));
% L = [1-sum(L);L]
% lam = sol(59+58+1:end)
% [~,~,pf] = exdem(sol,p);

%now reduce tariffs to zero (replication of AL experiment)
%options = optimset('Display','iter');
%sol_nt = ktrlink(@(x) fakeobj(x),sol(1:59),[],[],[],[],zeros(59,1),ones(59,1)*inf,@(x) exdem_no_tf(x,p,L,lam),options);

[w,L,lam]

[w_nt,pf_nt] = exdem_no_tf(w,L,lam,p);

wel_nt = ((w_nt./pf_nt)-(eps./pf))./(eps./pf);

tot_wel_nt = (sum((w_nt./pf_nt).*L)-sum((eps./pf).*L))./sum(eps./pf.*L);

%now keep tariffs, but let labor move freely

[w_lm,L_lm,pf_lm,eps_lm] = exdem_lm(w,L,lam,p);

wel_lm = ((eps_lm./pf_lm)-(eps./pf))./(eps./pf);

tot_wel_lm = (sum((eps_lm./pf_lm).*L_lm)-sum((eps./pf).*L))./sum(eps./pf.*L);

out_table = [eps./pf,lam,L*100,w_nt./pf_nt,eps_lm./pf_lm,L_lm*100,wel_nt*100,wel_lm*100];

col_labs = {'$\xi$/ pf','$\lam$','L','$w_{nt}$/$pf_{nt}$','$\xi_{lm}$/$pf_{lm}$','$L_{lm}$','$\Delta$wel\_nt','$\Delta$wel\_lm'};

row_labs = {'United States','Japan','Germany','France','United Kingdom','Italy',...
    'China','Brazil','Canada','Spain','Mexico','India','Australia','Netherlands',...
    'Russian Federation','Argentina','Switzerland','Belgium','Sweden','Austria',...
    'Turkey','Indonesia','Denmark','Hong Kong PRC','Norway','Thailand','Poland','Saudia Arabia',...
    'South Africa','Finland','Greece','Portugal','Israel','Iran','Colombia','Venezuela',...
    'Malaysia','Singapore','Ireland','Egypt','Philippines','Chile','Pakistan','New Zealand',...
    'Peru','Czech Republic','Algeria','Hungary','Ukraine','Bangladesh','Romania',...
    'Morocco','Nigeria','Vietnam','Belarus','Kazakhstan','Slovak Republic','Tunisia','Sri Lanka'};  

matrix2latex(out_table,'lat_table.tex','rowLabels',row_labs,'columnLabels',col_labs,'format','%.2f','alignment','c','size','small');

save results

toc;
diary off

