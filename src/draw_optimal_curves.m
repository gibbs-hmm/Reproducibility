function draw_optimal_curves(N,inputmean,stdev,ps,Nper)


%rng('default');  %uncomment if running multiple times


Nreps = choose_nreps(N);
SumValidN = 0;
[alpha,beta] = getalphabeta(inputmean,stdev.^2);
stdvalsforplot = zeros(1,length(Nreps));
stderrsforplot = zeros(1,length(Nreps));
p10valsforplot = zeros(1,length(Nreps));
p10errsforplot = zeros(1,length(Nreps));
fileid = fopen(['N',num2str(N),'_Mean',num2str(inputmean),'_StandardDeviation',num2str(stdev),'_PS',num2str(ps),'.txt'],'w');

parpool('threads')

for i=1:length(Nreps)
    disp(strcat(num2str(i),' out of  ',num2str(length(Nreps))))
    fprintf(fileid,[num2str(Nreps(1,i)),' Replicates\n\n']);
    pavgvector = zeros(1,Nper);
    pstdevvector = zeros(1,Nper);
    ppercentile10vector = zeros(1,Nper);
    nrep = Nreps(1,i);
    nsites1 = split_experiments(N,nrep);
    ValidN = Nper;
 
    parfor j = 1:Nper
    %for j = 1:Nper
        nsites = zeros(1,length(nsites1));
        numzeros = 0;
        for k = 1:length(nsites1)
            nsites(1,k) = binornd(nsites1(1,k),ps);
            if nsites(1,k) == 0
                numzeros = numzeros + 1;
            end
        end
        newnrep = nrep - numzeros;
        newnsites = nsites(nsites > 0);
        res = fakedata(newnrep,newnsites,alpha,beta,1);
        dat = res{1};
        
        if length(find(newnsites==dat)) + length(find(dat==0)) ~= length(dat)
            [pavg,pstdev,ppercentile10] = hierarchical_pred_inf_FAST(newnsites,dat,0,1200);
            %pavgvector(1,j) = pavg;
            pstdevvector(1,j) = pstdev;
            ppercentile10vector(1,j) = ppercentile10;
        else %if EVERYTHING is either 100% or 0%
            ValidN = ValidN - 1;
            pavgvector(1,j) = -Inf;
            pstdevvector(1,j) = -Inf;
            ppercentile10vector(1,j) = -Inf;
        end
    end
    
    SumValidN = SumValidN + ValidN;
    pstdevvector = pstdevvector(find(pstdevvector>=0));
    ppercentile10vector = ppercentile10vector(find(ppercentile10vector>=0));
    pstdmean = mean(pstdevvector);
    pstdstd = std(pstdevvector);
    pstderr = pstdstd./sqrt(ValidN);
    p10mean = mean(ppercentile10vector);
    p10std = std(ppercentile10vector);
    p10err = p10std./sqrt(ValidN);
    stdvalsforplot(1,i) = pstdmean;
    stderrsforplot(1,i) = pstderr;
    p10valsforplot(1,i) = p10mean;
    p10errsforplot(1,i) = p10err;
    fprintf(fileid,['Standard deviation of p has Mean = ',num2str(pstdmean),', Standard Deviation = ',num2str(pstdstd),', SEM = ',num2str(pstderr),'\n10th percentile of p has Mean = ',num2str(p10mean),', Standard Deviation = ',num2str(p10std),', SEM = ',num2str(p10err),', Valid N = ',num2str(ValidN),'\n\n\n']);
    
end
AVGPERCENTVALID = (SumValidN./length(Nreps))/Nper;
if AVGPERCENTVALID < 0.95
    fprintf(fileid,'WARNING: RESULTS MAY BE UNRELIABLE DUE TO HIGH MEAN/WIDE VARIANCE (SEE PAPER FOR DETAILS)\n');
end
fprintf(fileid,[num2str(Nper),' simulations per data point, valid simulations per data point has average = ',num2str(SumValidN./length(Nreps))]);
fclose(fileid); 
h = figure('visible', 'on');
errorbar(Nreps,stdvalsforplot,stderrsforplot,'r')
hold on
errorbar(Nreps,p10valsforplot,p10errsforplot,'b')
hold off
legend('Standard Deviation','10th Percentile','Location','East')
xlabel('Number of Replicates')
ylabel('Standard Deviation / 10th Percentile of Predictive Distribution')
titleline0 = 'Predictive Distribution Statistics:';
titleline1 = ['Total Number of Experiments = ',num2str(N),', Expected Mean = ',num2str(inputmean),','];
titleline2 = ['Expected Standard Deviation = ',num2str(stdev),', Expected % Successful Experiments = ',num2str(ps)];
title({titleline0,titleline1,titleline2});

highnrep = Nreps(end); 
lownrep = Nreps(1);
d = (highnrep-lownrep)./30;
xt = Nreps(1);
for i=2:length(Nreps)
    if Nreps(i) > xt(end) + d
        xt = [xt,Nreps(i)];
    end
end
set(gca,'XTick',xt)
print(h,'-dpdf',strcat('Optimal_Curve_N',num2str(N),'_Mean',num2str(inputmean),'_StandardDeviation',num2str(stdev),'_PS',num2str(ps),'.pdf'));
print(h,'-dpng','-r75', strcat('Optimal_Curve_N',num2str(N),'_Mean',num2str(inputmean),'_StandardDeviation',num2str(stdev),'_PS',num2str(ps),'.png'));

poolobj = gcp('nocreate');
delete(poolobj);
end

function values = choose_nreps(N)
count = 0;
val = floor(N/3);
for i=4:N
    newval = floor(N/i);
    if newval ~= val
        count = count + 1;
    end
    val = newval;
end
values = zeros(1,count);
ind = 1;
val = floor(N/3);
for i=4:N
   newval = floor(N/i);
   if newval ~= val
       values(1,ind) = i-1;
       ind = ind+1;
   end
   val = newval;
end
end


