function score = gridscore(ac, ind)
    try
        score = GridnessRadius(ac, FindROuter(ac));
    catch ME
        fprintf('%d: error calculting grid score %s\n', ind, ME.message);
        score = -2;
    end
end

function gridness2 = GridnessRadius(org_Acor ,R_outer)
    if length(org_Acor) == 1 || isnan(max(org_Acor(:)))
        gridness2 = -2;
        return %EXIT
    end
    
    % compute distance from center
    [val,center_y]=max(max(org_Acor));
    [val,center_x]=max(max(org_Acor'));
    [Y,X] = ndgrid(1:1:size(org_Acor,1), 1:1:size(org_Acor,2));
    dist_from_center=sqrt((Y-center_y).^2+(X-center_x).^2);
    % making sure that outer radius of the ring (R_outer) is not bigger than the distance matrix
    %if R_outer ~= -1%NOAM
    R_outer=min([min(dist_from_center(1,:)),min(dist_from_center(:,1)),...
        min(dist_from_center(size(dist_from_center,1),:))...
        min(dist_from_center(:,size(dist_from_center,2))),R_outer]);
    % compute inner radius of the anulus (ring)
    R_inner=ceil(min(dist_from_center(org_Acor<0.1))); %this is the problem <<
    
    % else %NOAM
    %    acg = imgaussfilt(org_Acor, 2,'FilterDomain','spatial'); %imshow
    %     [k,l] = find(imregionalmax(acg));
    %     dist = pdist2([l k],[center_x, center_y] ); [dist,ind]=sort(dist); %l=l(ind(2:7)); k=k(ind(2:7)); %CENTER FLIPPED?
    %     R_inner = dist(2)*0.5
    %     R_outer = dist(2)*1.5
    % end
    %extract the original anulus (ring) from Acor
    org_Ring=org_Acor(dist_from_center<=R_outer & dist_from_center>=R_inner);
    % make sure that after rotation and interpolation the center will remain the maximum point.
    
    org_Acor(center_x,center_y)=10;
    for jj = 2:6
        % rotate the auto-correlation
        rot_Acor=imrotate(org_Acor,(jj-1)*30,'bicubic');
        % compute distance from new center
        [val,tmp_center_x]=max(max(rot_Acor));
        [val,tmp_center_y]=max(max(rot_Acor'));
        [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
        tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
        % extract the anulus(ring)
        rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
        if length(rot_Ring)~=length(org_Ring)
            gridness2=nan;
            return
        end
        % compute pearson correlation between rotate Acor and original Acor
        corrValues(jj) = PointCorr(org_Ring,rot_Ring);
        clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X
    end
    % min of higher correlation at 60 and 120 degree rotation
    min_rot_60_120 = min([corrValues(3),corrValues(5)]);
    % max of lower correlations 30,90,150 rotation
    max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
    % calculate gridness min(60,120)-max(30,90,150)
    gridness2 = min_rot_60_120 - max_rot_30_90_150;
    % different way to calculate gridness
    gridness1=mean(([corrValues(3),corrValues(5)]))-...
        mean([corrValues(2),corrValues(4),corrValues(6)]);
    
    %{
    figure
    hold off;
    imagesc(org_Acor)
    set(gca,'ydir','normal')
    hold on
    %c =[length(org_Acor)/2 length(org_Acor)/2];
    c=[center_x center_y];
    viscircles(c,R_outer);
    viscircles(c,R_inner,'Color','b');
        %}
        
end


function out_mat=PointCorr(Rxx,RxxR)
    nan_mat=Rxx .* RxxR;
    notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
    n=length(notnan_inds);
    if n < 2
        out_mat = NaN;
    end
    sigma_x_y =sum(nan_mat(notnan_inds));
    sigma_x =      sum(Rxx(notnan_inds));
    sigma_y =      sum(RxxR(notnan_inds));
    sigma_x2 = sum(Rxx(notnan_inds).^2);
    sigma_y2 = sum(RxxR(notnan_inds).^2);
    
    out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
        sqrt(n*sigma_x2-sigma_x.^2) ./ ...
        sqrt(n*sigma_y2-sigma_y.^2);
end


function R_outer = FindROuter(acorr)
    if length(acorr) == 1 || isnan(max(acorr(:)));
        R_outer = -1;
        return;   %EXIT
    end
    % calculate all the extrema points in the spatial autocorrelation
    %[zmax,imax,zmin,imin]= Extrema2(acorr);%NOAM
    %[i,j]=ind2sub(size(acorr),imax); %NOAM
    [i,j] = find(imregionalmax(acorr));
    %put all extrema points in dist
    dist(:,1)=j;
    dist(:,2)=i;
    n=length(i);
    %calculate the distance of all extrema to the central peak and put them in
    %column 3
    dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
    % sort the hexonal peaks by distance to the centeral peak
    [score,ind]=sort(dist(:,3));
    dist=dist(ind,:);
    %zmax=zmax(ind);
    R=dist(2,3);
    count=1;
    i=2;
    hex_peaks(1,:,:)=dist(1,:,:);
    % finds the first 6 closest peaks to the central peak
    while count<7 && i<=length(dist)
        % calculate the min distance of point i from all other already chosen
        % points
        min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
            (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
        % point i needs to be on the right side (we choos only half the point cause its semetrical)
        % and the distance of point i from all other already chosen points
        % needs to be higher than R/2
        if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
            hex_peaks(count,:,:)=dist(i,:,:);
            count=count+1;
            hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
            hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);
            
            count=count+1;
        end
        i=i+1;
    end
    R_outer=max(hex_peaks(:,3))*1.15;
    %{
figure
imagesc(acorr)
set(gca,'ydir','normal')
hold on
viscircles([50 50],R_outer)
%plot([50 50+R_outer],[50 50])
    %}
end