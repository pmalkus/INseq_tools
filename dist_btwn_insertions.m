function dist_btwn_insertions(data)

%find distribution of distances betewenn insertion sites

%add begining and end?
use = 0;
use(2:length(data)+1) = data;
%end?

sep = []; %distances between sites
for i=1:length(use)-1
    sep(i) = use(i+1)-use(i);
end

f1 = hist(sep, 40);
f2 = cdfplot(sep);