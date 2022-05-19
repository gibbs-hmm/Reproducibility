%Create a contour plot from a scatter plot:

function scatter2contour(xvals,yvals)
gridnum = 200;
%create 100x100 grid:
Xmax = max(xvals);
Xmin = min(xvals);
Ymax = max(yvals);
Ymin = min(yvals);
Xgrid = linspace(Xmin-0.1,Xmax+0.1,gridnum);
Ygrid = linspace(Ymin-0.1,Ymax+0.1,gridnum);
%Xgridforcontour falls in midpoints of Xgrid
Xforcontour = linspace((Xgrid(1)+Xgrid(2))/2,(Xgrid(end)+Xgrid(end-1))/2,gridnum-1);
Yforcontour = linspace((Ygrid(1)+Ygrid(2))/2,(Ygrid(end)+Ygrid(end-1))/2,gridnum-1);
%Make grid of squares w/ numbers of scatter points in each square...
Scatternum = zeros(gridnum-1,gridnum-1);
for i = 1:gridnum

    xind = findxindex(Xgrid,xvals(1,i));
    yind = findxindex(Ygrid,yvals(1,i));
    Scatternum(xind,yind) = Scatternum(xind,yind) + 1;
end

windowWidth = 70;
halfWidth=windowWidth/2;
gaussFilter = gausswin(windowWidth);
gaussFilter = gaussFilter/sum(gaussFilter);

Z = conv2(gaussFilter,gaussFilter,Scatternum,'same');


maxm = max(max(Z));
v = 0.05*maxm:0.1*maxm:0.95*maxm;
h = figure('visible', 'off');
contour(Xforcontour,Yforcontour,Z',v);
   
xlabel('Mean of Beta Distribution');
ylabel('Log_{10} Standard Deviation of Beta Distribution')
title('Posterior Beta Distribution')
colorbar('YTick',[.05*maxm:0.1*maxm:0.95*maxm],'YTickLabel',{'5%','15%','25%','35%','45%','55%','65%','75%','85%','95% of max'})
print(h,'-dpdf',strcat('Contour Representation of Posterior Beta Distribution.pdf'));
print(h,'-dpng','-r80', 'Contour_Representation_of_Posterior_Beta_Distribution.png');
  
end
%----------------------------------------------------------------
%Given Xgrid, find where an x-value lies (above what index)
function ind = findxindex(Xgrid,xval)
ind = 2;
while xval > Xgrid(1,ind) && ind < length(Xgrid)
    ind = ind + 1;
end
ind = ind - 1;

end
%----------------------------------------------------------------    
    