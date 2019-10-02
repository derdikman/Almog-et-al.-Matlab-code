function module = Find_Module(acorr)
    
module.major_ax = [];
module.minor_ax = [];
module.angle = -1;
module.x0 = -1;
module.y0 = -1;
module.x = 0;
module.y = 0;
module.exists = 0;
    
    
%REBEKKAHS
%module.hex_peaks2 = find_six_points(acorr);
%noam

[k,l] = find(imregionalmax(acorr));
dist = pdist2([l k], flip(size(acorr))/2); [~,ind]=sort(dist);
if length(l) >= 7
    l=l(ind(2:7)); k=k(ind(2:7)); 
end
hex_peaks = [l k];
module.hex_peaks = [l k];    

if length(hex_peaks) == 6
    
    % fitting the elipse to the 6 peaks that we have found
    [major_ax, minor_ax, x0, y0, angle] = ellipse_fit(hex_peaks(1:6,1), hex_peaks(1:6,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% fixing a bug in ellipse fit %%%%%%%%%%%%%%%%
    [x1,y1]=ellipse(major_ax,minor_ax,angle,x0,y0);%'k',100);
    [x2,y2]=ellipse(major_ax,minor_ax,-angle,x0,y0);%'k',100);
    
    %%%%%%%%%%%%%%%%%%%%% calculate which ellipse fits better (distance to points)
    i=1:length(x1);
    d1=0;d2=0;
    for j=1:6
        d1=d1+min((hex_peaks(j,1)-x1(i)).^2+(hex_peaks(j,2)-y1(i)).^2);
        d2=d2+min((hex_peaks(j,1)-x2(i)).^2+(hex_peaks(j,2)-y2(i)).^2);
    end
    module.x = x1; module.y = y1;
    %%%%%%%%%%%%%%%%% if d2 smaller put -angel and the bug is fixed
    if d2<d1
        angle=-angle;
        module.x = x2; module.y = y2;
    end
    
    if ~isreal(major_ax) ||  ~isreal(minor_ax)
        major_ax=nan;
        minor_ax=nan;
    end
    
    hex_peaks(7,1)=x0;
    hex_peaks(7,2)=y0;
    
    module.major_ax = major_ax;
    module.minor_ax = minor_ax;
    module.angle = angle;
    module.hex_peaks = hex_peaks;
    module.x0 = x0;
    module.y0 = y0;
    module.exists = 1;
end
disp('')
end

function [semimajor_axis, semiminor_axis, x0, y0, phi] = ellipse_fit(x, y)
%
% Programmed by: Tal Hendel <thendel@tx.technion.ac.il>
% Faculty of Biomedical Engineering, Technion- Israel Institute of Technology
% 12-Dec-2008

x = x(:);
y = y(:);

%Construct M
M = [2*x.*y y.^2 2*x 2*y ones(size(x))];

% Multiply (-X.^2) by pseudoinverse(M)
e = M\(-x.^2);

%Extract parameters from vector e
a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);

%Use Formulas from Mathworld to find semimajor_axis, semiminor_axis, x0, y0
%, and phi

delta = b^2-a*c;

x0 = (c*d - b*f)/delta;
y0 = (a*f - b*d)/delta;

phi = 0.5 * acot((c-a)/(2*b));

%phi = 0.5 * atan2((2*b),(c-a));

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);

if (a_prime < b_prime)
    phi = pi/2 - phi;
end

disp('')
end
function [x,y]=ellipse(major,minor,phi,x_center,y_center)
    sinphi = sin(phi);
    cosphi = cos(phi);
    alpha =0:pi/100:2*pi;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
    x = x_center + (major * cosalpha * cosphi - minor * sinalpha * sinphi);
    y = y_center + (major * cosalpha * sinphi + minor * sinalpha * cosphi);
end





%UNUSED
function [x,y,h]=ellipse2(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same.
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value.
% If no color is specified, it makes automatic use of the colors specified by
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%

% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original
% written by Peter Blattner, Institute of Microtechnology, University of
% Neuchatel, Switzerland, blattner@imt.unine.ch


% Check the number of input arguments

if nargin<1,
    ra=[];
end;
if nargin<2,
    rb=[];
end;
if nargin<3,
    ang=[];
end;

%if nargin==1,
%  error('Not enough arguments');
%end;

if nargin<5,
    x0=[];
    y0=[];
end;

if nargin<6,
    C=[];
end

if nargin<7,
    Nb=[];
end

% set up the default values

if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
%if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);

%if isstr(C),C=C(:);end;

if length(ra)~=length(rb),
    error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
    error('length(x0)~=length(y0)');
end;

% how many inscribed elllipses are plotted

if length(ra)~=length(x0)
    maxk=length(ra)*length(x0);
else
    maxk=length(ra);
end;

% drawing loop

for k=1:maxk
    
    if length(x0)==1
        xpos=x0;
        ypos=y0;
        radm=ra(k);
        radn=rb(k);
        if length(ang)==1
            an=ang;
        else
            an=ang(k);
        end;
    elseif length(ra)==1
        xpos=x0(k);
        ypos=y0(k);
        radm=ra;
        radn=rb;
        an=ang;
    elseif length(x0)==length(ra)
        xpos=x0(k);
        ypos=y0(k);
        radm=ra(k);
        radn=rb(k);
        an=ang(k)
    else
        rada=ra(fix((k-1)/size(x0,1))+1);
        radb=rb(fix((k-1)/size(x0,1))+1);
        an=ang(fix((k-1)/size(x0,1))+1);
        xpos=x0(rem(k-1,size(x0,1))+1);
        ypos=y0(rem(k-1,size(y0,1))+1);
    end;
    
    co=cos(an);
    si=sin(an);
    the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
    x=radm*cos(the)*co-si*radn*sin(the)+xpos;
    y=radm*cos(the)*si+co*radn*sin(the)+ypos;
    h = [];
    %TO DRAW
    %warning('off','MATLAB:plot:IgnoreImaginaryXYPart')
    %h(k)=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
    %warning('on','MATLAB:plot:IgnoreImaginaryXYPart')
    %set(h(k),'color',C(rem(k-1,size(C,1))+1,:));
    %close; %noam
    
end
end

function return_value=Distance(x,y,x2,y2)
return_value = sqrt((x2-x).^2+(y2-y).^2);
end
function [six_orientation_pts] = find_six_points(autocorr)

%find autocorrelation mat maximum points
[auto_max_inds] = FindMaxIndsRateMap(autocorr);

if length(auto_max_inds) >= 7
    
    [size_x, size_y]=size(autocorr);
    
    %finds distances from center of auto_corr_map to all max peaks
    cen = 1:length(auto_max_inds);
    auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),(size_x/2)+0.5,(size_y/2)+0.5);
    
    
    
    %find point closest to center
    min_distance_index = find(auto_distances==min(auto_distances));
    middle_pt(1) = auto_max_inds(min_distance_index(1),1); %NZA Added(1) because rarely there are two vals
    middle_pt(2) = auto_max_inds(min_distance_index(1),2); %same
    
    %repeats find distances with more accurate center point
    cen = 1:length(auto_max_inds);
    
    auto_distances = Distance(auto_max_inds(cen, 1),auto_max_inds(cen, 2),middle_pt(1),middle_pt(2));
    
    [new_distances, six_inds] = sort(auto_distances (:));
    
    % finds the points for the module ellipse
    % finds the 6 closest pts to the center
    auto_dist_inds1 = find(auto_distances == new_distances(2));
    auto_dist_inds2 = find(auto_distances == new_distances(3));
    auto_dist_inds3 = find(auto_distances == new_distances(4));
    auto_dist_inds4 = find(auto_distances == new_distances(5));
    auto_dist_inds5 = find(auto_distances == new_distances(6));
    auto_dist_inds6 = find(auto_distances == new_distances(7));
    
    union1 = union(auto_dist_inds1, auto_dist_inds2);
    union2 = union(auto_dist_inds3, auto_dist_inds4);
    union3 = union(auto_dist_inds5, auto_dist_inds6);
    union1 = union(union1, union2);
    auto_dist_inds = union(union1, union3);
    
    %
    six_orientation_pts=nan(length(auto_dist_inds),2);
    for k= 1:length(auto_dist_inds);
        six_orientation_pts (k,1) = auto_max_inds(auto_dist_inds(k),1);
        six_orientation_pts (k,2) = auto_max_inds(auto_dist_inds(k),2);
    end
    %}
    
    %adds center point
    
    six_orientation_pts (7,1) = middle_pt(1);
    six_orientation_pts (7,2) = middle_pt(2);
    
    % to check for accuracy:
    % figure; imagesc(autocorr); hold on;
    % plot(six_orientation_pts(:,2), six_orientation_pts(:,1), 'x')
    
    disp('');
else
    six_orientation_pts = -1;
end

end
function [max_inds] = FindMaxIndsRateMap(rate_mat)

max_inds = 0;

% turn nans in rate_mat to zeros
rate_mat(isnan(rate_mat))=0;

% pad rate mat with zero edges to catch border maxs
rm_size= size(rate_mat);
rm_size= rm_size+2;
rate_mat_new= zeros(rm_size);
rate_mat_new(2:end-1, 2:end-1)= rate_mat(1:end, 1:end);
rate_mat=rate_mat_new;

max_inds_len = 0;

[size_x, size_y]= size(rate_mat);

for fig_i = 2:size_x-1
    for j = 2:size_y-1
        if rate_mat(fig_i,j) > rate_mat(fig_i+1,j) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i-1,j) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i,j+1) && ...
                rate_mat(fig_i,j) > rate_mat(fig_i,j-1)
            %hold on %WHY????????????
            max_inds_len = max_inds_len+1;
            max_inds(max_inds_len,1) = fig_i-1;     %list of maximum pt indices
            max_inds(max_inds_len,2) = j-1;
        end
    end
end
end