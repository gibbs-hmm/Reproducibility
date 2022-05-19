%Given a total number of experiments and number of desired replicates, give
%a vector of experiments per replicate (if doesn't divide evenly, make them
%as close as possible - numbers only differ by 1 or 0)
function Nperrep = split_experiments(Ntotal,Nrep)

val1 = floor(Ntotal/Nrep);
numdiff = Ntotal-val1*Nrep;
Nperrep = [repmat(val1,1,Nrep-numdiff),repmat(val1+1,1,numdiff)];

end