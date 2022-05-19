

%This script does the following:
%1) It computes a contour plot on a user-defined grid that represents the
%posterior distribution of the beta distribution parameters on transformed
%axes.
%2) Given this contour plot, it draws NUMSAMPLES samples of p|alpha,beta
%3) Given the samples of p, it draws NUMSAMPLES samples of k|p,N=100

function [pavg,pstdev,ppercentile10] = hierarchical_pred_inf_FAST(trialnum,trialsuc,verbose,NUMSAMPLES)

%Trialnum is a vector of sample size per replicate
%Trialsuc is a vector of successful validations per replicate
%Example: If you have 3 replicates with 20 sites each, then
%Trialnum=[20,20,20] and Trialsuc could be something like [15,16,12]

%NUMSAMPLES is the number of samples used for predictive inference -
%default for single run should be 10,000

%gridnum determines how fine the grid is - currently it is 200x200. If
%things are taking too long to run, it can be set to 100 basically without
%any differences in output


gridnum=200;
xstepsize = 0.1;
ystepsize = 0.1;
[xcenter,ycenter] = estimatecenters(trialnum,trialsuc);


%X and Y are the axes for the contour plot
X = xcenter-xstepsize*floor(gridnum/2):xstepsize:xcenter+xstepsize*ceil(gridnum/2-1);
Y = ycenter-ystepsize*floor(gridnum/2):ystepsize:ycenter+ystepsize*ceil(gridnum/2-1);
 
Z = makeZ(trialnum,trialsuc,xcenter,ycenter,xstepsize,ystepsize,gridnum,verbose);


%FIX CONTOUR PLOT SO IT'S PERFECT:

%does the contour plot go off the edges?
overlap = checkoverlap(Z,gridnum);


if overlap == 0
    %expand plot to fill window
    [xc,yc,xstepsize,ystepsize] = expand(Z,gridnum,X,Y);
	Z = makeZ(trialnum,trialsuc,xc,yc,xstepsize,ystepsize,gridnum,verbose); 
    overlap = checkoverlap(Z,gridnum);
end
if overlap == 1
    %increase stepsize until it no longer goes off the edges
    while overlap == 1
        xstepsize = xstepsize + .1*xstepsize;
        ystepsize = ystepsize + .1*ystepsize;
        Z = makeZ(trialnum,trialsuc,xcenter,ycenter,xstepsize, ystepsize,gridnum,verbose);
        overlap = checkoverlap(Z,gridnum);
    end
    %take new Z and increased stepsizes and expand to fit the window
    X = xcenter-xstepsize*floor(gridnum/2):xstepsize:xcenter+xstepsize*ceil(gridnum/2-1);
    Y = ycenter-ystepsize*floor(gridnum/2):ystepsize:ycenter+ystepsize*ceil(gridnum/2-1);
    [xc,yc,xstepsize,ystepsize] = expand(Z,gridnum,X,Y);
    Z = makeZ(trialnum,trialsuc,xc,yc,xstepsize,ystepsize,gridnum,verbose);

end  
%for future calculations, this has the final centers and stepsizes
X = xc-xstepsize*floor(gridnum/2):xstepsize:xc+xstepsize*ceil(gridnum/2-1);
Y = yc-ystepsize*floor(gridnum/2):ystepsize:yc+ystepsize*ceil(gridnum/2-1);
    

    
%APPROXIMATE THE POSTERIOR DISTRIBUTION OF ALPHA AND BETA BY SAMPLING
%FROM THE GRID IN PROPORTION TO Z, THEN JITTERING UNIFORMLY (ALSO EXPLAINED IN GELMAN)

%get marginal distribution of log(alpha/beta) (horizontal axis)
M = zeros(1,gridnum);
for i = 1:gridnum
    sum1 = 0;
    for j = 1:gridnum
        sum1 = sum1 + Z(i,j);
    end
    M(1,i) = sum1;
end
M = M./sum(M);

%DRAW NUMSAMPLES SAMPLES:
sampledx = zeros(1,NUMSAMPLES);
sampledy = zeros(1,NUMSAMPLES);
for k = 1:NUMSAMPLES
    %draw log(alpha/beta)
    r = rand;
    cumul = M(1,1);
    i = 1;
    while r > cumul
        i = i+1;
        cumul = cumul + M(1,i);
    end
    alphaoverbeta_index = i;
    %draw log(alpha+beta)
    condit = Z(i,:)/sum(Z(i,:));
    r = rand;
    cumul = condit(1,1);
    i = 1;
    while r > cumul
        i = i+1;
        cumul = cumul + condit(1,i);
    end
    alphaplusbeta_index = i;
    alphaoverbeta = X(1,alphaoverbeta_index);
    alphaplusbeta = Y(1,alphaplusbeta_index);


    %jitter:
    %xgrid = xstepsize;
    %ygrid = ystepsize;
    alphaoverbeta = alphaoverbeta - xstepsize./2 + xstepsize*rand;
    alphaplusbeta = alphaplusbeta - ystepsize./2 + ystepsize*rand;
    sampledx(1,k) = alphaoverbeta;
    sampledy(1,k) = alphaplusbeta;


end
%figure(11);
%plot(sampledx,sampledy,'.')


%put parameters back to un-transformed axes
sampledalpha = zeros(1,length(sampledx));
sampledbeta = zeros(1,length(sampledx));

for i=1:length(sampledx)
    [alpha,beta] = makealphabeta(sampledx(1,i),sampledy(1,i));
    sampledalpha(1,i) = alpha;
    sampledbeta(1,i) = beta;
end
meanofbetadistributions = zeros(1,length(sampledx));
logstdofbetadistributions = zeros(1,length(sampledx));
for i=1:length(sampledx)
    meanofbetadistributions(1,i) = sampledalpha(1,i)/(sampledalpha(1,i)+sampledbeta(1,i));
    logstdofbetadistributions(1,i) = log10(sqrt(sampledalpha(1,i)*sampledbeta(1,i)/((sampledalpha(1,i)+sampledbeta(1,i))^2*(sampledalpha(1,i)+sampledbeta(1,i)+1))));
end

if verbose == 1
    scatter2contour(meanofbetadistributions,logstdofbetadistributions);
end


%for each (alpha,beta) draw a proportion p
props = zeros(1,length(sampledalpha));
for i = 1:length(props)
    props(1,i) = betarnd(sampledalpha(1,i),sampledbeta(1,i));
end
%sample mean, standard deviation, print standard deviation
pvariance = 0;
pavg = 0;
for i=1:length(props)
    pavg = pavg + props(1,i);
end
pavg = pavg./length(props);
for i=1:length(props)
    pvariance = pvariance + (props(1,i)-pavg)^2;
end
pvariance = pvariance/(length(props)-1);
pstdev = sqrt(pvariance);
if verbose == 1
    fileid = fopen('Predictive_Distribution_Statistics.txt','w');
    fprintf(fileid, ['Number of Experiments: ', num2str(trialnum), '\n']); % BT
    fprintf(fileid, ['Number of Validations: ', num2str(trialsuc), '\n']);
    fprintf(fileid,['Mean of p: ',num2str(pavg),'\n']);
    fprintf(fileid,['Standard Deviation of p: ',num2str(pstdev),'\n']);
end
%Get 10th percentile:
sortedprops = sort(props);
ppercentile10 = sortedprops(round(.1*length(props)));
ppercentile90 = sortedprops(round(.9*length(props)));
ppercentile05 = sortedprops(round(.05*length(props)));
ppercentile95 = sortedprops(round(.95*length(props)));
ppercentile025 = sortedprops(round(.025*length(props)));
ppercentile975 = sortedprops(round(.975*length(props)));

%Histogram of sampled p:
if verbose == 1
    g = figure('visible', 'on'); % BT
    hist(props,25);
    h = findobj(gca,'Type','patch', 'visible', 'on');
    set(h,'EdgeColor','b','FaceColor','b');
    xlabel('Sampled Proportions')
    ylabel('Relative Frequency')
    title('Sampled Proportions from Posterior (alpha,beta) Distribution')
    %[histbottom,histtop] = ylim;
    hold on
    h95 = plot([ppercentile025,ppercentile025],ylim,'-r');
    hold on
    plot([ppercentile975,ppercentile975],ylim,'-r');

    hold on
    h90 = plot([ppercentile05,ppercentile05],ylim,'-g');
    hold on
    plot([ppercentile95,ppercentile95],ylim,'-g');

    hold on
    h10 = plot([ppercentile10,ppercentile10],ylim,'-k');
    hold on
    plot([ppercentile90,ppercentile90],ylim,'-k');

    hold off
    legend([h95,h90,h10],'95% Confidence Limits','90% Confidence Limits','80% Confidence Limits','Location','NorthWest');
    Yticks = get(gca,['y','Tick']);
    newticks = Yticks./NUMSAMPLES;
    newticks = strread(num2str(newticks),'%s');
    set(gca,['y','TickLabel'],newticks);
    print(g,'-dpdf',strcat('Sampled Proportions from Posterior Beta Distribution.pdf'));
    print(g,'-dpng', '-r75', 'Sampled_Proportions_from_Posterior_Beta_Distribution.png');
    %set(gcf,'ResizeCmd','fixLabels');
end

if verbose == 1   
    fprintf(fileid,['10th percentile of p: ',num2str(ppercentile10),'\n']);
    fprintf(fileid,['0.95 Confidence Interval: (',num2str(ppercentile025),', ',num2str(ppercentile975),')\n']);
    fprintf(fileid,['0.90 Confidence Interval: (',num2str(ppercentile05),', ',num2str(ppercentile95),')\n']);
    fprintf(fileid,['0.80 Confidence Interval: (',num2str(ppercentile10),', ',num2str(ppercentile90),')\n']);
    fclose(fileid);
end 
end

function [ALPHA_MATRIX,BETA_MATRIX] = makealphabetamatrix(X,Y,gridnum)
    %FIND A WAY TO VECTOR-IZE
    z = ones(1,gridnum);
    XPLUSY = X'*z+z'*Y;
    XMATRIX = repmat(X',1,gridnum);
    YMATRIX = repmat(Y,gridnum,1);
    ALPHA_MATRIX = exp(XPLUSY)./(exp(XMATRIX)+1);
    BETA_MATRIX = exp(YMATRIX)./(exp(XMATRIX)+1);
end
function [alpha,beta] = makealphabeta(a,b)
    alpha = exp(a+b)./(exp(a)+1);
    beta = exp(b)./(exp(a)+1);
end

function Z = makeZ(trialnum,trialsuc,xcenter,ycenter,xstepsize,ystepsize,gridnum,verbose)
%X and Y are the axes for the contour plot
X = xcenter-xstepsize*floor(gridnum/2):xstepsize:xcenter+xstepsize*ceil(gridnum/2-1);
Y = ycenter-ystepsize*floor(gridnum/2):ystepsize:ycenter+ystepsize*ceil(gridnum/2-1);

if any(~isfinite(X)) || any(~isfinite(Y))
    disp('not finite')
end

[ALPHAMATRIX,BETAMATRIX] = makealphabetamatrix(X,Y,gridnum);
same = trialnum == trialsuc;
sum1s = zeros(gridnum,gridnum);
if length(same) == length(find(same==1))
    for k=1:length(trialnum)
        sum1s = sum1s + gammaln(ALPHAMATRIX+BETAMATRIX) + gammaln(ALPHAMATRIX + trialsuc(k)) - gammaln(ALPHAMATRIX) - gammaln(ALPHAMATRIX+BETAMATRIX+trialnum(k));
    end
else
  for k=1:length(trialnum)
    sum1s = sum1s + gammaln(ALPHAMATRIX+BETAMATRIX) + gammaln(ALPHAMATRIX + trialsuc(k)) + gammaln(BETAMATRIX + trialnum(k) - trialsuc(k)) - gammaln(ALPHAMATRIX) - gammaln(BETAMATRIX) - gammaln(ALPHAMATRIX+BETAMATRIX+trialnum(k));
  end  
end
    
Z = sum1s + log(ALPHAMATRIX)+log(BETAMATRIX)-(5./2)*log(ALPHAMATRIX+BETAMATRIX);
Zcontour = Z';

%Zcontour is for the contour plot, Z is for the rest of the analysis
%(Zcontour is Z reflected over i=j - this is due to some stupid matlab
%idiosyncracy)


%exponentiate:
maximum = max(max(Zcontour));
Zcontour = Zcontour - maximum;
Zcontour = exp(Zcontour);
Z = Z - maximum;
Z = exp(Z);

if verbose == 1
    %create contour plot:
    v = 0.05:0.1:0.95;
    figure(1);
    contour(X,Y,Zcontour,v);
    xlabel('log(alpha/beta)');
    ylabel('log(alpha+beta)');
    title('Posterior Distribution of Hyperparameters, Transformed Axes') 
    colorbar
end
        
end

function [xc,yc,xstep,ystep] = expand(Z,gridnum,X,Y)
L = gridnum+1;
B = gridnum+1;
R = 0;
T = 0;
%find L,R,T,B (leftmost non-zero point, etc...):
%SPED UP - NO NESTED FOR LOOPS
for i=1:gridnum
    Lind = find(Z(:,i)>=0.05,1,'first');
    if Lind < L
        L = Lind;
    end
    Rind = find(Z(:,i)>=0.05,1,'last');
    if Rind > R
        R = Rind;
    end
end
L = X(L);
R = X(R);
for j=1:gridnum
    Tind = find(Z(j,:)>=0.05,1,'last');
    Bind = find(Z(j,:)>=0.05,1,'first');
    if Bind < B
        B = Bind;
    end
    if Tind > T
        T = Tind;
    end
end
B = Y(B);
T = Y(T);


%Calculate the boundaries of the new plot so that the curve takes up
%8/10ths of the left-right and top-bottom space 
top = (9/8)*T - (1/8)*B;
bottom = B + T - top;
right = (9/8)*R - (1/8)*L;
left = L + R - right;
xc = (left+right)/2;
yc = (top+bottom)/2;
%new stepsizes
xstep = (right-left)/gridnum;
ystep = (top-bottom)/gridnum;
end

function overlap = checkoverlap(Z,gridnum)
%Go around the boundary and check for non-zero elements (here non-zero
%means <0.05, since that is the lowest contour
overlap = 0;
i = 1;
while overlap == 0 && i <= gridnum
    if Z(i,1) >= 0.05 || Z(i,gridnum) >= 0.05 || Z(1,i) >= 0.05 || Z(gridnum,i) >= 0.05
        overlap = 1;
    end
    i = i+1;
end 

  
end
