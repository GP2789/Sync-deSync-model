function [CI_L, CI_U, M] = bootstrap( data, rep, bins)

    h2 = waitbar(0, 'Time Points', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);
    x = size(data,2); y = size(data,1); m = zeros(rep, x); x2 = 0:x-1;
    if(isempty(bins)==1); bins = ones(1, x)*y; end
    
    for b=1:rep % choose random sets of data
        m(b,:) = mean(data(ceil(rand(y,x).*bins) + x2*y),1);
        waitbar(b/rep,h2);
    end
    m(isnan(m)==1)=0;
    for xi = 1:x % sort data through time points
        m(:,xi) = sort(m(:,xi));
    end
    m = m(ceil(rep*0.05)+1:floor(rep-rep*0.05),:); % cut of top and bottom 5%
    CI_L = m(1,:); CI_U = m(end,:); M = mean(m,1);
    close(h2);
    
end