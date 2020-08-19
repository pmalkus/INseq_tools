%get insertions per locus

%from uniqueinidices

%take difference between entries in col2 and col1

mui=myUniqueIndices;
mus=mui(:,2)-mui(:,1);
myTAhits=mus+1;

ui=uniqueindices;
us=ui(:,2)-ui(:,1);
TAhits=us+1;