%N and trialsuc are both vectors, N is the vector of sample sizes, trialsuc
%is the vector of successful validations

function [xc,yc] = estimatecenters(N,trialsuc)
min_var = 1.0e-06;  % BT

propv = trialsuc./N;
valid = 0;
firstval = propv(1,1);
i = 2;
while valid == 0 && i <= length(N)
    if propv(1,i) ~= firstval
        valid = 1;
    end
    i = i + 1;
end
if valid == 0
    if trialsuc(1,1) == N(1,1)
        trialsuc(1,1) = trialsuc(1,1) -1;
    else
        trialsuc(1,1) = trialsuc(1,1) + 1;
    end
end
cont = 0;
while cont == 0
    exp_p = mean(trialsuc./N);
    var_p = max(var(trialsuc./N), min_var);  %BT 
    alphaplusbeta = exp_p*(1-exp_p)/var_p - 1;
    thing1 = round(10000*(exp_p*(1-exp_p)));
    thing2 = round(10000*var_p);
    if thing1 ~= thing2
        cont = 1;
    else
        if trialsuc(1,length(trialsuc)) == N(1,length(trialsuc))
            trialsuc(1,length(trialsuc)) = trialsuc(1,length(trialsuc)) - 1;
        else
            trialsuc(1,length(trialsuc)) = trialsuc(1,length(trialsuc)) + 1;
        end
    end
end
xc1 = alphaplusbeta*exp_p;
yc1 = alphaplusbeta*(1-exp_p);
xc = log(xc1/yc1);
yc = log(xc1+yc1);
if ~isfinite(xc) || ~isfinite(yc)
    disp('not finite')
end
end