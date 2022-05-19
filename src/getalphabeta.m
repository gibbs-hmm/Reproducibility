%take mean and variance of desired distribution, and find alpha, beta

function [alpha,beta] = getalphabeta(mu,vari)
alpha = ((mu^2*(1-mu))/vari) - mu;
beta = ((1-mu)*alpha)/mu;
end