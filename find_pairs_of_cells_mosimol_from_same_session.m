%%%%%%%%%

% N Z A % 

% grid cells: 228
% pairs: 749

%%% DATASTATS() !!

%%%%%%%%%


function find_pairs_of_cells_mosimol


params.dir_load =...
    'C:\Noam\Data\muscimol\DB_musc_MEC\';
    %'C:\Users\noam\Desktop\proj\muscimol\DB_musc_MEC\';

params.dir_save =....
    'C:\Noam\Output\muscimol\noam_output\';
    %'C:\Users\noam\Desktop\proj\muscimol\noam_output';
%
tic
 files = dir(strcat(params.dir_load,'DB*.mat'));
 for i=1 %1:length(files) %HOW MUCH DATA TO LOAD
    data{i} = load(strcat(params.dir_load,files(i).name));
    fprintf('%.1f \n',(i*100/length(files)));
 end
 toc
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data = load(strcat(params.dir_load,'ALL_CELLS.mat')); 
%data = data.data;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

grid_cells = 0;
for i=1:length(data) %a = load(strcat(params.dir_load,files(i).name));
    a = data(i); a = a{1}.db;
    if (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ')) % '_GRID' ?????      
        for j=i+1:length(data)
            b = data(j); b = b{1}.db;
            if ( strcmp(a.rat, b.rat) && strcmp(a.date, b.date)...
                    && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.area, b.area)  ...  %area????
                    && ( a.tetrode ~= b.tetrode || a.cell ~= b.cell) ...
                )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                    %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i; 
                pairs_by_file_index(grid_cells,2) = j; 
            end
        end
    else
        disp (a.cell_type_subj);
    end
end
disp(grid_cells);
disp('fin!');

%spike plot
%
a = data{1}.db.B(2) % B session is combination of all recordings. has 3 structs: before, during, after muscomol injection
%figure();
%plot( mean([a.pos_data.x1, a.pos_data.x2] ,2) , mean([a.pos_data.y1, a.pos_data.y2] ,2) );
%hold on;
%plot( mean([a.spike_data.x, a.spike_data.x2] ,2) , mean([a.spike_data.y, a.spike_data.y2] ,2), 'r.' );

temp = data{1}.db.B(2);
c1.pos_x =  mean([temp.pos_data.x1, temp.pos_data.x2], 2)';
c1.pos_x = c1.pos_x - min(c1.pos_x) + 0.00001;
c1.pos_y =  mean([temp.pos_data.y1, temp.pos_data.y2], 2)';
c1.pos_y = c1.pos_y - min(c1.pos_y) + 0.00001;

c1.pos_time = temp.pos_data.t;
c1.spike_x =  mean([temp.spike_data.x, temp.spike_data.x2], 2)';
c1.spike_x = c1.spike_x - min(c1.spike_x) + 0.00001;
c1.spike_y =  mean([temp.spike_data.y, temp.spike_data.y2], 2)';
c1.spike_y = c1.spike_y - min(c1.spike_y) + 0.00001;
c1.spike_time = temp.spike_data.ts;

figure();
plot(c1.pos_x, c1.pos_y);
hold on
plot(c1.spike_x, c1.spike_y,'.');

% fix spike position NaN's
for i = 1:length(c1.spike_time);
    if isnan(c1.spike_x(i)) %assumes NaN's happen in x and y at same time
        %find closest time value to spike time, use that as x and y
        [~, ind] = min(abs(c1.pos_time(:)-c1.spike_time(i)));
        c1.spike_x(i) = c1.pos_x(ind);
        c1.spike_y(i) = c1.pos_y(ind);  
    end
end

plot( c1.pos_x ,c1.pos_y );
hold on;
plot( c1.spike_x, c1.spike_y, 'r.' );

%%%%%
num_bins = 100;
rate_map_time =   zeros(num_bins, num_bins);
rate_map_spikes = zeros(num_bins, num_bins);
for i = 1:length(c1.pos_time);
    col =            ceil(c1.pos_x(i)/max(c1.pos_x) * num_bins);
    row = num_bins - ceil(c1.pos_y(i)/max(c1.pos_y) * num_bins) +1; %flips y
    rate_map_time(row, col) = rate_map_time(row, col) + 1;
end

for i = 1:length(c1.spike_time);
    col =            ceil(c1.spike_x(i)/max(c1.spike_x) * num_bins);
    row = num_bins - ceil(c1.spike_y(i)/max(c1.spike_y) * num_bins) + 1;
    rate_map_spikes(row, col) = rate_map_spikes(row, col) + 1;
    % because some spikes are in other positions by interpolation
    if(rate_map_time(row, col) < rate_map_spikes(row, col))
         rate_map_time(row, col) = rate_map_time(row, col) + 1;
    end
    %rate_map_time(row, col) = rate_map_time(row, col) + 1; %% IS THIS VALID? double counting???
end
rate_map = (rate_map_spikes ./ rate_map_time);
rate_map(isnan(rate_map)) = 0;
max(rate_map(:));
rate_map = rate_map/max(rate_map(:)); %normalize

% GAUSS SMOOTHE
 N = num_bins;
 sigma = 2;
 [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 f=f./sum(f(:));
 figure();
 imagesc(conv2(rate_map,f,'same'));
 colormap jet;
%
%figure('Position', [100, 100, 400, 400]);
figure();
imshow(imgaussfilt(rate_map, 1));
colormap jet;

imshow([0,0.2,0.3;0.4,0.5,0.6;0.7,0.8,0.9;-1,0.1,-10])
colormap jet;


% CORRELATION
imagesc(corrcoef(rate_map));
imagesc(xcorr2(rate_map,rate_map));

a = xcorr2(rate_map);
a(a > 0.1) = 0.1; imagesc(a);
figure();
imshow(imgaussfilt( a, 1));
colormap jet;

%
for i=1:length(pairs_by_file_index)
    
end


            %{
            ift strcmp(rat{i},rat{j}) && strcmp(date{i},date{j}) && strcmp(session{i},session{j})...
                    && (tetrode(i)~=tetrode(j) || cell(i)~=cell(j)) && ...
                    (p_gridness(i)<0.05  && p_gridness(j)<0.05)
            %} 

end

function f=gaussian2d(N,sigma)
  % N is grid size, sigma speaks for itself
 [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 f=f./sum(f(:));
end
 

            