function [rate_mat, max_r] = createRateMap(px, py, pt, sx, sy, st, newmethod, nbins)
    if nargin==2
        nbins= py;
        pt=    px.pt;
        py=    px.py;
        sx=    px.sx;
        sy=    px.sy;
        st=    px.st;
        px=    px.px;
        newmethod = true;
    end
    
    if newmethod
        if ~isempty(st)
            px = toCol(px); py = toCol(py); %pt = toCol(pt); 
            sx = toCol(sx); sy = toCol(sy); %st = toCol(st); 
            %si = discretize(st, [-Inf; mean([pt(2:end) pt(1:end-1)],2); +Inf]);
            mx = max(px);
            my = max(py);
            pxi = discretize(px, 0:mx/nbins:mx);
            pyi = discretize(py, 0:my/nbins:my);
            rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
            rmt(rmt<1) = 1e6; %NO NANS??
            sxi = discretize(sx, 0:mx/nbins:mx); %px(si)
            syi = discretize(sy, 0:my/nbins:my); %py(si)
            rms = accumarray([syi,sxi], 1, [nbins nbins]);
            %!!!!
            rm = rms./(rmt*.02);
            max_r = max(rm(:));
            t = rms-rmt; t(t<0)=0; rmt = rmt + t; %add in extra timestep ONLY when spikes occured faster than timestep
            if sum(t)~=0
                disp('$$ STILL ADDING TIME TO BINS');
            end
            %max_r = (max(t(:))+1)/.02;
            rate_mat = (rms ./ rmt*0.02); 
            rate_mat(isnan(rate_mat)) = 0; %DO THIS?
        else
            rate_mat = zeros(nbins);
        end
    else
        [rate_mat, max_r] = Create_Rate_Map(px, py, pt, sx, sy, st, new, nbins)
    end
end

function [rate_mat, max_r] = Create_Rate_Map(px, py, pt, sx, sy, st, new, nbins)
    posx = px; posy = py; post = pt;
    spkx = sx; spky = sy; spkt = st;
    parms.sigma = 3; % gaussiam smoothing factor
    parms.time_per_bin=0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
    % Minimum radius used in the auto-correlogram when finding the best
    parms.bin_size = 3; % size of spacial bin (for create the rate map)

    min_x = min(floor(posx)); max_x = ceil(max(posx));  
    min_y = min(floor(posy)); max_y = ceil(max(posy));

    % divide the environment into spatial bins
    axis_x = min_x : parms.bin_size : max_x;
    axis_y = min_y : parms.bin_size : max_y;

    time_mat =  zeros(length(axis_y),length(axis_x));
    spike_mat = zeros(length(axis_y),length(axis_x));
    rate_mat =  zeros(length(axis_y),length(axis_x));

    %create time mat (compute how much time the rat spends in each bin)
    % find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
    for i = 1:length(post)
        if ~isnan(posx(i)) && ~isnan(posy(i))
            [min_val,x_ind] =  min(abs(posx(i) - axis_x));
            [min_val,y_ind] =  min(abs(posy(i) - axis_y));
            time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind) + parms.time_per_bin;
            %if ~isnan(min_val)
        end
    end
    %create count mat( count the num of spikes in each bin)
    for i = 1:length(spkt)
        if ~isnan(spkx(i)) && ~isnan(spky(i))
            [min_val,x_ind] =  min(abs(spkx(i) - axis_x));
            [min_val,y_ind] =  min(abs(spky(i) - axis_y));
            spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
            %if ~isnan(min_val)
        end
    end
    % create rate mat
    rate_mat=spike_mat./time_mat;
    rate_mat(rate_mat==inf)=NaN;
    rate_mat_unsmoothed = rate_mat;
    %create window
    h=fspecial('gaussian',3*[parms.sigma,parms.sigma],parms.sigma);
    rate_mat = nanconv2(rate_mat,h);
    max_r = max(rate_mat(:));
    rate_mat(isnan(rate_mat))=0; %bad??????nza
end

function out_mat = nanconv2(mat,h)
out_mat = mat;
nan_mat = isnan(mat);
% dilate nan_mat
SE = strel('disk', 2);
nan_mat =  ~imdilate(~nan_mat,SE);
i_size = size(h,1); j_size = size(h,2);
[work_mat,npad_i,npad_j] = pad_edges(mat,h,2);
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        % for each i and j, choose the correct sub-mat (size of h) to multiply with h
        sub_mat = work_mat(npad_i+i-floor(i_size/2):npad_i+i+floor(i_size/2), ...
            npad_j+j-floor(j_size/2):npad_j+j+floor(j_size/2)  ); % assumes h is odd in number
        
        notnan_inds = find(~isnan(sub_mat));
        if ~isempty(notnan_inds)
            sum_h = sum(h(notnan_inds));   % normalize to the places without a NaN
            out_mat(i,j) = nansum(nansum(sub_mat .* h));
            out_mat(i,j) = out_mat(i,j)/sum_h;
        end
    end % for j
end % for i
out_mat(nan_mat) = NaN;
disp('')
end


