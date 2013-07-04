function [W,ISR] = fcombi(X,varargin)
% fcombi - FCOMBI algorithm. See [1] for details.
%
% Usage:
%   [W, ISR] = fcombi(X [,'parameter name',parameter value])
%
% Input/output parameters:
%   W   - (m x m) Estimated mixing matrix
%   ISR - (m x m) Estimated ISR matrix
%   X   - (m x T) Observed mixtures (samples are columnwise)
%
% Optional input parameters:
%   'ar_order' - AR order for EWASOBI (default is: 10)
%
% Reference:
% [1] Gomez-Herrero, G., Koldovsky, Z., Tichavsky, P., Egiazarian, K.,
% "A fast algorithm for Blind Separation of Non-Gaussian and Time-Correlated
% sources", In Proc. EUSIPCO 2007, Poznan, Poland, 2007.
%
% See also
%   IWASOBI, EFICA

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

% CONTANTS
MINQOF = 10; % reasonable values in the range 10-20.
% default arguments
p = 10;  % order of the AR models

% process varargin
i = 1;
while i <= length(varargin),
    argok = 1;
    if ischar(varargin{i}),
        switch lower(varargin{i}),
            % argument IDs
            case 'ar_order', i = i+1; p = varargin{i};
        end
    else
        argok = 0;
    end
    if ~argok,
        disp(['(fcombi) Ignoring invalid argument #' num2str(i+1)]);
    end
    i = i+1;
end

% preprocessing (make data have zero mean and unit variance)
X = X-mean(X,2)*ones(1,length(X));
CC = cov(X')^(-1/2);
X = CC*X;
d = size(X,1);

% run WASOBI
[Wwa, tmp, ISRwa] = iwasobi(X,p,.99);

% search for clusters within WASOBI ISR matrices
[owa,cluswa] = hclus(ISRwa,MINQOF);

% run EFICA on each still unresolved cluster
W = Wwa;
ISR = ISRwa;
ord = owa;
clus = cluswa;
if size(clus,1)<d,
for i = 1:size(clus,1),
    if clus(i,2)>1,
        cindex = ord(clus(i,1):clus(i,1)+clus(i,2)-1);
        [Wef,ISRef] = efica(W(cindex,:)*X,eye(length(cindex)));
        [oef,clusef] = hclus(ISRef,MINQOF);
        % update the de-mixing matrix only if EFICA did better than WASOBI
        % old line
        %if size(clusef,1)>1 || (mean(mean(ISR(cindex,cindex)))>mean(mean(ISRef))),            
        % new line:
        if size(clusef,1)>1 || (min(sum(ISR(cindex,cindex),2))>min(sum(ISRef,2))),
            % now run again WASOBI in each unresolved cluster
            for j = 1:size(clusef,1),
                if clusef(j,2)>1,
                    cindex2 = oef(clusef(j,1):clusef(j,1)+clusef(j,2)-1);
                    [Wwa,tmp,ISRwa2] = iwasobi(Wef(cindex2,:)*W(cindex,:)*X,p,.99);
                    [owa,cluswa] = hclus(ISRwa2,MINQOF);
                    % old line
                    %if size(cluswa,1) > 1 || (mean(mean(ISRef(cindex2,cindex2)))>mean(mean(ISRwa2)))
                    % new line:
                    if size(cluswa,1) > 1 || (min(sum(ISRef(cindex2,cindex2),2))>min(sum(ISRwa2,2)))
                        Wef(cindex2,:) = Wwa*Wef(cindex2,:);
                        ISR(cindex2,cindex2) = ISRwa2;
                    end
                end
            end
            W(cindex,:) = Wef*W(cindex,:);
            ISR(cindex,cindex) = ISRef;
        end
    end
end
end
% dewhiten
W=W*CC;

% -------------------------------------------------------------------------
% HCLUS - Clustering sub-function
% -------------------------------------------------------------------------
function [order,clus_out,qof] = hclus(ISR,dd)
% Agglomerative hierarchical clustering with single linkage strategy

d = size(ISR,1);

% distance matrix
D = -.5*(ISR+ISR');
D = D-diag(diag(D));
%D = -ISR;
% initial clusters are the individual components
cluster = 1:d;

% merge clusters until there is no more clusters to merge
for i = 2:d,
    % initialize the clusters
    cluster(i,:) = cluster(i-1,:);
    % find pair of clusters that are closest to each other
    [mval,mrow] = min(D);
    [mval,mcol] = min(mval);
    mrow = mrow(mcol);
    % merge the clusters
    C1 = find(cluster(i,:)==cluster(i,mrow));
    C2 = find(cluster(i,:)==cluster(i,mcol));
    cindex = min(cluster(i,[C1 C2]));
    cluster(i,[C1 C2]) = cindex;
    % Make sure that we will not try to merge them again
    D(C1,C2) = Inf;
    D(C2,C1) = Inf;
end

% select the best partition level (THIS IS THE CRUCIAL STEP!)
% Note: This is a bit based on heuristics. Furthermore, two thresholds need
% to be settled in order to detect the situation where there is only 1
% cluster and the situation where all clusters are 1-dimensional.
D = -ISR;
Dmin = min(min(D));
qof = zeros(1,d);
for i = 2:d-1,
    cindex = unique(cluster(i,:));
    numt = 0;    
    nintra = 0;    
    for j = 1:length(cindex)
        Cin = find(cluster(i,:)==cindex(j));
        nin = length(Cin);        
        if length(Cin) == 1,
            num = 0;                
        else            
            num = sum(sum(D(Cin,Cin)));                        
            nintra = nintra+nin.^2-nin;
        end
        numt = numt + num;        
    end
    % the QOF is a ratio between the average distance within clusters and
    % the average clusters between clusters. A high QOF means a good fit.
    qof(i) = (numt/nintra)/(((sum(D(:))-numt)/(d^2-nintra-d)));    
end
% find local maxima of the QOF function
qofdiff = diff(qof(2:end-1));
v1 = qofdiff>0;
v1 = [0 0 v1(2:end)-v1(1:end-1)];
maxloc = find(v1<0);
[maxval,index] = sort(qof(maxloc),'descend');
maxloc = maxloc(index);
% simple way of discarding spurious local maxima
maxdiff = zeros(1,length(maxloc));
for i = 1:length(maxloc),
    maxdiff(i) = .5*(-qof(maxloc(i)-1)+2*qof(maxloc(i))-qof(maxloc(i)+1));
end
if ~isempty(maxloc),
    maxloc(maxdiff<.75*max(maxdiff))=[];
    % select the maximum from the local maxima
    plevel = maxloc(1);
    qof = qof(plevel);
else
    [qof,plevel] = max(qof);
end
% check wether we have 1-dimensional clusters (perfect separation)
if qof < -.1/(Dmin),
    plevel = 1;
    qof = -.1/(Dmin);
end
% check wether we have just 1 cluster (all sources are still mixed)
% A reasonable value of dd would be in the range: 10-20
if qof < dd,
    plevel = d;
    qof = dd;
end

% ordering of components and setting up the output format in the same way
% as in function hcsort
[cluster,order] = sort(cluster(plevel,:));
cindex = unique(cluster);
D = D(order,order);
clus_out = zeros(length(cindex),4);
for i = 1:length(cindex)
    tmp = find(cluster==cindex(i));
    clus_out(i,1) = tmp(1);
    clus_out(i,2) = tmp(end)-tmp(1)+1;
    Cin = find(cluster==cindex(i));
    Cout = setdiff(1:d,Cin);
    if length(Cin) == 1,
        num = -1;
        den = min(min(D(Cin,Cout)));
        dlen = 1;
        clus_out(i,3) = sum(sum(ISR(order(Cin),order(Cout))))*(d-1)/(dlen*(d-dlen));
    elseif isempty(Cout),
        den = -1;
        num = max(max(D(Cin,Cin)));
        clus_out(i,3) = NaN;
    else
        num = max(max(D(Cin,Cin)));
        den = min(min(D(Cin,Cout)));
        dlen = length(Cin);
        clus_out(i,3) = sum(sum(ISR(order(Cin),order(Cout))))*(d-1)/(dlen*(d-dlen));
    end
    clus_out(i,4) = num/den;
end