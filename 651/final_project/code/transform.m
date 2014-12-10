function scale=transform(vectors,by,newscale,oldscale)
if ~exist('newscale','var')                              % Default to standardize normal
    newscale='normal';
end
if ~exist('by','var')                                    % Default to transform by columns
    by='col';
end
if ~strcmp(newscale,'normal') && length(newscale) ~= 2   % If invalid 'type' set to 'normal'
    newscale='normal';
end
if ~(strcmp(by,'row') || strcmp(by,'col'))               % If invalid 'by' set to 'col'
    by='col';
end
if strcmp(by,'row')                                      % If 'row' then transpose matrix
    vectors = vectors';
end
n = size(vectors, 2);                                    % 'n' is number of columns
scale=zeros(size(vectors));
if ~exist('oldscale','var')
    oldscale=zeros(n,2);
    for ii=1:n
        oldscale(ii,1)=min(vectors(:,ii));
        oldscale(ii,2)=max(vectors(:,ii));
    end
end
for ii=1:n
    if strcmp(newscale,'normal')
        scale(:,ii)=(vectors(:,ii)-mean(vectors(:,ii)))/...
            std(vectors(:,ii));
    else
        if range(vectors(:,ii)) == 0
            scale(:,ii)=vectors(:,ii);
        else
            scale(:,ii)=(vectors(:,ii)-oldscale(ii,1))*...
                range(newscale)/range(oldscale(ii,:))+newscale(1);
        end
    end
end
if strcmp(by,'row')
    scale = scale';
end
    
end
