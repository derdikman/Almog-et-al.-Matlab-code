function [gridness2,cang] = gridscore2(ac, sigma) %R_inner,R_outer,[cenx ceny]
if length(ac) == 1 || sum(~isnan(ac(:)) )<25
    'bad ratemap nonnan<25' 
    gridness2 = 0; cang = 0;
else
    %try
    %finding peaks
    acg = imgaussfilt(ac, sigma,'FilterDomain','spatial');ac = acg; %INSERTED HERE!!
    
    %only look for max in middle area of matrix
    ms = min(size(ac))/2; mss = ms/4; mss = round((-mss:1:mss) + ms); %box around center of matrix?
    t = ac; t(mss,mss) = t(mss,mss)+ max(ac(:));
    [cmaxy, cmaxx]=ind2sub(size(ac),find(t==max(t(:)) ) );%center on max to make ring? assumes one max
    
    %assert(len(cmaxy)==1);
    if len(cmaxy)~=1
        disp('error in gridscore: len(cmaxy)~=1, ring set to center.');
        t = round(flip(size(ac)/2));
        cmaxx=t(1);cmaxy=t(1);
        %gridness2 = 0; cang = 0;
        %return;
    end
    
    ac(cmaxy,cmaxx)= max(ac(:))+0.01; %set to max incase there are others out of middle area;
    %acg = imgaussfilt(ac, sigma,'FilterDomain','spatial');
    acg(isnan(acg))=0; % only used to find peaks.
    [py,px] = find(imregionalmax(acg)); %find peacks of gaussian smoothed autocorrelation
    dist = pdist2([px py],[cmaxx cmaxy] ); [dist,ind]=sort(dist); % distances of peaks to center
    px=px(ind);py=py(ind);
    cang = 0; t = round(flip(size(ac)/2)); cax=t(1); cay=t(2); cx1=px(1);cy1=py(1);
    
    %[px(1) py(1) cmaxx cmaxy cax cay]
    if length(dist) >= 2
        %cang = rad2deg(atan2(py(2)-ceny, px(2)-cenx));
        cang = rad2deg(atan2(py(1)-cay, px(1)-cax));
         %{
        % PLOT FIRING FIELDS)
        figure(7221); clf;
        hold off;
        imagesc(acg);colormap jet;
        set(gca,'ydir','normal') %FLIP XY EHEN GRAPHING
        hold on
        %plot(cenx, ceny,'mo')
        %a = [1:len(px)]'; b = num2str(a); nstr = cellstr(num2str(1:len(px)));
        text(px(2:7), py(2:7), cellstr(num2str([1:len(px(2:7))]')),'color','w');
        plot(px(1:2), py(1:2),'mo');
        plot(cmaxx, cmaxy,'go')
        plot(px(2:7), py(2:7),'ko');
        plot(cax,cay,'ro'); axis square;
        %}
    end
    
    
    if length(dist) >= 7
        R_inner = dist(2)*0.5; %inner ring dist(2) is the closest peak to the center
        R_outer = min(min(size(ac))/2-1,dist(7)*1); %outer ring
        %making ring
        [Y,X] = ndgrid(1:1:size(ac,1), 1:1:size(ac,2));
        dist_from_center=sqrt((Y-cmaxy).^2+(X-cmaxx).^2); %grid of distances of indexes to the center;
        %extract the original anulus (ring) from Acor
        org_Ring=ac(dist_from_center<=R_outer & dist_from_center>=R_inner);
        %{
        figure(92322); clf;
        hold off; colormap jet;
        imagesc(acg)
        set(gca,'ydir','normal'); axis square;
        hold on
        plot(cmaxx, cmaxy,'mo')
        viscircles([cmaxx, cmaxy],R_outer);
        viscircles([cmaxx, cmaxy],R_inner,'Color','b');
        %}
        % make sure that after rotation and interpolation the center will remain the maximum point.
        %ac(cmaxx,cmaxy)= max(ac(:)) + max(ac(:))*0.2; %NOT SURE THIS WORKS
        [y x] = meshgrid(1:len(ac), 1:len(ac)); del=(x-cmaxx).^2+(y-cmaxy).^2<2^2; %rad=b;
        ac(del)=max(ac(:))*1.1; ac(cmaxx,cmaxy)=ac(cmaxx,cmaxy)*1.1;
        for i = 2:6
            % rotate the auto-correlation
            rot_Acor=imrotate(ac,(i-1)*30,'bicubic'); %max(rot_Acor(:))
            % compute distance from new center
            [~,tmp_center_x]=max(max(rot_Acor)); 
            [~,tmp_center_y]=max(max(rot_Acor'));
            [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
            tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
            % extract the anulus(ring)
            rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
            if length(rot_Ring)~=length(org_Ring)
                gridness2=0;
                return
            end
            % compute pearson correlation between rotate Acor and original Acor
            corrValues(i) = PointCorr(org_Ring,rot_Ring);
            clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X
        end
        % min of higher correlation at 60 and 120 degree rotation
        min_rot_60_120 = min([corrValues(3),corrValues(5)]);
        % max of lower correlations 30,90,150 rotation
        max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
        % calculate gridness min(60,120)-max(30,90,150)
        gridness2 = min_rot_60_120 - max_rot_30_90_150;
        %different way to calculate gridness
        %gridness1=mean(([corrValues(3),corrValues(5)]))-...
        %mean([corrValues(2),corrValues(4),corrValues(6)]);
    else
        gridness2 = 0;
        fprintf('%bad cell sigma%d gridscore -2',sigma);
    end
    %catch ME
    %fprintf('%d: error calculting grid score %s\n', ind, ME.message);
    %gridness2 = -2;
    %end
    
end
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

