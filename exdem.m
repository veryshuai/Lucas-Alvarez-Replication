function out = sys(w,L,lam,p)
%This function calculates excess demand for Lucas Alvarez model

alp = p.alp;
bet = p.bet;
the = p.the;
eta = p.eta;

A = gamma(1+the*(1-eta));
B = bet^-bet * (1-bet)^(-1+bet);




end