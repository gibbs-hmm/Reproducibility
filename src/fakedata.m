%generate data for given Nrep (Number of replicates, a scalar), Nsites
% (Number of sites per replicate, a vector) for given alpha and beta. Generate M independent data sets

function out = fakedata(Nrep,Nsites,alpha,beta,M)
% reset(RandStream.getDefaultStream,sum(100*clock));
out =cell(1,M);

    for i=1:M
        outv = zeros(1,Nrep);
        for j=1:Nrep
            p = betarnd(alpha,beta);
            k = binornd(Nsites(1,j),p);
            outv(1,j) = k;
        end

    out{i} = outv;
    end
end




