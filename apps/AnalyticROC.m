%File: AnalyticROC.m
%Author: Chris Tralie
%Purpose: Given a dSquared (Es/Sigma_n^2) value, compute the 
%ROC curve ((Pf, Pd) pairs) by integrating the conditional density
%functions of the likelihood ratio
function [Pf, Pd] = AnalyticROC(dSquared)
    variance = dSquared;
    mu = dSquared/2;
    NPoints = 400;
    k = -2*dSquared:4*dSquared/NPoints:2*dSquared;
    Pf = 0.5*erfc((k+mu)./sqrt(2*variance));
    Pd = 0.5 + 0.5*erf((mu-k)./sqrt(2*variance));
end