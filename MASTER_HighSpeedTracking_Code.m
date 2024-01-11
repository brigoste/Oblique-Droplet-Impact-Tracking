%
close all
clc
clear all
%

% ---------------------HEADER------------------------------------------
% Author: Brigham C Ostergaard
% Date: 6/13/2023
% Title: MASTER_HighSpeedTracking_Code.m
% Info: This code is intended to track the droplet spread on hydrophobic
% surfaces. It has the capability of tracking the droplets extreme-most
% points (top, bottom, left, right) and interpolating the area based on an
% ellipsoid fit of these points the major/minor diameters made with them
%
%Notes: This is untested on SH surfaces and also on angled impact.
%
% Along with saving figures of spread over time and area, I also
% incorrporated a "Master Save File" for the data across all experiments.
% It will add a line to the file every time we run this code and I am ready
% to save (boolean controlled). It will give max diameters, timestamps for
% each of the maxes, and the testing conditions (temp, surf, angle, etc.)
% It doesn't check to see if that test was saved previously.
%
% 6/23/23
% I have incorrporated a measure of non-dimentional diameters, a graph for
% them, and added them to the save file along with the total surface time
% and slide distance.

%--------------Testing Parameters----Edit between dotted lines-------------
test_path = "D:\Grad Research - HDD\RESEARCH COMBINED TESTING\0 Degrees\Smooth\WE 40\80 C\T2\";    %Change name here
images_dir = dir(fullfile(test_path,'*.tif'));

binary_threshold = 50;
start_fn = 112;
end_fn = 213;
start_n = 1;               %number of first frame in the names, change here too

save_both = false;           %Save graphs and save data to file
Test_num = 1;
Surface_type = 1;       %1 = smooth, 55 = 55 cf, 86 = 86 cf
Temp = 80;               %0 = room temp, 60 = 60, 80 = 80
angle = 0;              %0,15,30,45,60
WE = 40;                %40, 100, 140

old_cam = true;

%Background Subtraction parameters
subtract_back = true;
first_im = false;       %True = First image filter, False = Average image filter
toggle_first_last = true;      %Use last image vs first image as filter

if(start_fn > end_fn)
    temp = start_fn;
    start_fn = end_fn;
    end_fn = temp;
end

recalibrate = false;
% pix2mm = 0.064206137674026;             %0 degrees
% pix2mm = 0.082278670645751;             %15 degrees
pix2mm = 0.077218273339431;             %30 degrees


impact_param = [binary_threshold, start_fn, end_fn];
photo_dim = [1024,1024,1];


%Saving parameters
if save_both == true
    save_tracking = true;
else 
    save_tracking = false;
end
% save_tracking = false;
% save_fig_name = "15_degrees_We_40_80_Smooth_T2";         %Chnage name here
save_path = append(test_path,'Analysis\');
if(save_tracking)
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
end

if save_both == true
    save_data = true;      %Whether to write analysis results to master file
else
    save_data = false;
end
datafile = "D:\Grad Research - HDD\KINGSTON\Tests\Data.txt";     %Save in top file of directory
datafile2 = "D:\Grad research - HDD\KINGSTON\Tests\Data2.txt";   %New file to hold additional info


test_Date = 21062023;    %M/DD/YYYY --> MM/DD/YYYY if Oct->Dec
fps = 3000;

dt = 1/fps;            %1/fps of camera;

d0 = 5;            %Initial droplet diameter for non-dimensionalizing

%--------------------------------------------------------------------------
%------------------------Calibration Points--------------------------------
if(recalibrate)
    calibration_file_path_0_degrees = "D:\Grad Research - HDD\KINGSTON\Tests\Brighams_Tests(6_6_2023 - 6_7_2023)\0 degrees\0_degrees_calibration.tif";
    calibration_file_path_15_degrees = "D:\Grad Research - HDD\KINGSTON\Tests\Brighams_Tests(6_6_2023 - 6_7_2023)\15 degrees\15_degrees_calibration.tif";
    calibration_file_path_30_degrees = "D:\Grad Research - HDD\KINGSTON\Tests\Brighams_Tests(6_6_2023 - 6_7_2023)\30 degrees\30_degrees_calibration.tif";

    figure('name','Calibration Pane')
    if angle == 0
        imshow(calibration_file_path_0_degrees);                                                            %change for other angles
    elseif angle == 15
        imshow(calibration_file_path_15_degrees);
    elseif angle == 30
        imshow(calibration_file_path_30_degrees);
    end
    title('Select along two measurment lines next to each other at same disntace from rulers edge')
    
    [x,y] = ginput(2);

    dist = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);

    pix2in = 0.1/dist; % in/pix
    in2mm = 25.4;

    pix2mm = pix2in * in2mm;            %This is the calibration distance.

end

%--------------------------------------------------------------------------
threshold = impact_param(1);
start_fn = impact_param(2) - start_n;
end_fn = impact_param(3) - start_n;

n_images = abs(end_fn-start_fn)+1;
Area_vec = zeros(n_images,1);

impact_frames = zeros(photo_dim(1),photo_dim(2),photo_dim(3),n_images);
impact_frames = uint8(impact_frames);

impact_bw = uint8(zeros(photo_dim(1),photo_dim(2),n_images));

if(~old_cam)
    for i = 1:n_images
        impact_frames(:,:,:,n) = uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name)));
        impact_bw(:,:,i) = rgb2gray(impact_frames(:,:,:,i));
    end
else
    for i = 1:n_images
        impact_frames(:,:,i) = uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name)));
        impact_bw(:,:,i) = impact_frames(:,:,i);
    end
end

if(subtract_back)
    if(first_im == true)
        if toggle_first_last == true
            if(old_cam)
                background = imread(append(images_dir(end).folder,'\',images_dir(end).name));
            else
                background = rgb2gray(uint8(imread(append(images_dir(end).folder,'\',images_dir(end).name))));
            end
        else
            if(old_cam)
                background = imread(append(images_dir(1).folder,'\',images_dir(1).name));
            else
                background = rgb2gray(uint8(imread(append(images_dir(1).folder,'\',images_dir(1).name))));
            end

        end
        n = 0;
        for i = start_fn:end_fn
            n = n+1;
            if(~old_cam)
                impact_frames(:,:,:,n) = uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name)));
                impact_bw(:,:,n) = rgb2gray(uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name))));
            else
                impact_frames(:,:,n) = imread(append(images_dir(i).folder,'\',images_dir(i).name));
                impact_bw(:,:,n) = impact_frames(:,:,1);
            end
        end
    elseif(first_im == false)
        background = zeros(photo_dim(1),photo_dim(2));
        background = uint8(background);
        
        show_background_prog = false;
        
        n = 0;
        if(show_background_prog)
            figure('name','background creation')
        end
        for i = start_fn:end_fn
            n = n+1;
            if(~old_cam)
                impact_frames(:,:,:,n) = uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name)));
                impact_bw(:,:,n) = rgb2gray(uint8(imread(append(images_dir(i).folder,'\',images_dir(i).name))));
            else
                impact_frames(:,:,n) = imread(append(images_dir(i).folder,'\',images_dir(i).name));
                impact_bw(:,:,n) = impact_frames(:,:,1);
            end
            background = background + impact_bw(:,:,n)/n_images;
            if(show_background_prog)
                imshow(background)
                pause(0.2)
                clf
            end
        end
    end
    n = 0;
    for i = start_fn:end_fn
        n = n+1;
        impact_bw(:,:,n) = impact_bw(:,:,n) - background;
    end
end

%show the images


show_progress = false;

if(show_progress)
    figure("name","Droplet impact")
    for i = 1:n_images
        clf
        subplot(2,1,1)
        imshow(impact_bw(:,:,i));
        subplot(2,1,2)
        imshow(append(images_dir(i+start_fn).folder,'\',images_dir(i+start_fn).name));
%         pause(0.25)
    end
end

binary_images = zeros(photo_dim(1),photo_dim(2),n_images);

show_binary_comp = false;

%convert to binary?
for i = 1:n_images
    binary_images(:,:,i) = impact_bw(:,:,i)>threshold;
end
if(show_binary_comp)
    figure("name","Binary Comparison")
    for i = 1:n_images
        clf
        subplot(2,1,1)
        imshow(impact_bw(:,:,i));
        subplot(2,1,2)
        imshow(binary_images(:,:,i));
        pause(0.25)
    end
end


edge_track = figure('name','edge_track');

lateral_spread = zeros(n,1);
longitudinal_spread = zeros(n,1);
center =  zeros(n,2);

for i = 1:n_images
    edge_pos = Find_Edge(binary_images(:,:,i));

%     imshow(binary_images(:,:,i));
%     imshow(impact_bw(:,:,i));
    imshow(impact_frames(:,:,:,i));
  

    lateral_spread(i) = abs(edge_pos(2,1)-edge_pos(4,1));
    longitudinal_spread(i) = abs(edge_pos(1,2)-edge_pos(3,2));

    center(i,:) = [edge_pos(3,2)-floor(longitudinal_spread(i)/2),edge_pos(4,1) + floor(lateral_spread(i)/2)];

    hold on
    plot(edge_pos(:,2),edge_pos(:,1),'ro','LineWidth',3)
    plot(center(i,1),center(i,2),'go','LineWidth',2)
    plotEllipse(lateral_spread(i),longitudinal_spread(i),center(i,:),'r',3); %majorAxes,minorAxes,centers,color,width

    Area_vec(i,1) = lateral_spread(i) * longitudinal_spread(i) * pi * (pix2mm)^2;

    title(append('Image ', num2str(i)))
    hold off

    if(save_tracking)
        save_name = append(num2str(i,'%.5d'),'.tif');
        save_path_temp = append(save_path,'\Tracking Images\');
        if ~exist(save_path_temp)
            mkdir(save_path_temp)
        end
        save_track_path = append(save_path_temp,save_name);
        saveas(edge_track,save_track_path)
    end

    pause(0.5)
    if i < n_images; clf; end
end

plot_frames = 1:1:n_images;

lateral_spread = lateral_spread*pix2mm;
longitudinal_spread = longitudinal_spread*pix2mm;

%non-dimensionalize by d0.
nd_lat_spread = lateral_spread./d0;
nd_long_spread = longitudinal_spread./d0;


max_y = ceil(max(max([longitudinal_spread,lateral_spread])));
min_y = floor(min(min([longitudinal_spread,lateral_spread])));

spread_graph = figure('name','Spread Dynamics');
hold on
plot(plot_frames,lateral_spread,'-ob')
plot(plot_frames,longitudinal_spread,'-or')
xlabel('frame')
ylabel('spread distance (mm)')
ylim([min_y, max_y]);
legend('lateral spread','longitudinal spread','location','south')

max_y = max(max([nd_lat_spread, nd_long_spread]));

nd_spread_graph = figure('name','nd Spread Dynamics');
hold on
plot(plot_frames,nd_lat_spread,'--ob')
plot(plot_frames,nd_long_spread,'--or')
xlabel('frame')
ylabel('spread distance (mm)')
ylim([0, 1.1*max_y]);
legend('lateral spread','longitudinal spread','location','south')

center2mm = center.*pix2mm;

centroid_tracking_graph = figure('name','Centroid Track');
subplot(2,1,1)
plot(plot_frames,center2mm(:,1))
ylabel('Position along axis (mm)')
legend('Centroid x','location','south')
subplot(2,1,2)
plot(plot_frames,center2mm(:,2),'r')
hold off
legend('Centroid y','location','south')
xlabel('Frame')
ylabel('Position along axis (mm)')

y_translation = abs(center2mm(1,2) - center2mm(end,2));
x_translation = abs(center2mm(1,1) - center2mm(end,1));


if(save_tracking)
    save_name = append(save_path,'Spread Diameters.jpg');
    saveas(spread_graph,save_name);
    save_name = append(save_path,'ND Spread Diameters.jpg');
    saveas(nd_spread_graph,save_name);
    save_name = append(save_path,'Centroid Tracking.jpg');
    saveas(centroid_tracking_graph,save_name);
end

max_y = ceil(max(Area_vec));
min_y = floor(min(Area_vec));

area_graph = figure('name','Area Approximation');
plot(plot_frames,Area_vec,'o-r')
xlabel('Frame')
ylabel('Area (mm^2)')
ylim([0.9*min_y,1.1*max_y]);
legend('Estimated Area of elliptical spread','location','south')

if(save_tracking)
    save_name = append(save_path,'Area Tracking.jpg');
    saveas(area_graph,save_name);
end

% Now, we have all the data from the tests, we want to save the following:
%   Max Diameter
%   Time to max Diameter
%   Testing Parameters(Test number, surface, WE, Temperature, Angle, date)

if(save_data)
    %Extract the time and distnace values for maximum spread and initial
    %diameter in both lateral/longitudinal directions
    [max_long_spread,index_long] = max(longitudinal_spread);
    [min_long_spread,index_min_long] = min(longitudinal_spread);
    [max_lat_spread,index_lat] = max(lateral_spread);
    [min_lat_spread,index_min_lat] = min(lateral_spread);
    max_long_spread_D0 = max(nd_long_spread);
    max_lat_spread_D0 = max(nd_lat_spread);

    %determine time to max spreading based on impact and fps of camera
    time_long = abs(index_long-index_min_long)*dt;
    time_lat = abs(index_lat - index_min_lat)*dt;

    %Determine toatl slide distance and total slide time
    slide_dist = sqrt(y_translation^2 + x_translation^2);        %This is the hypotenus of the triangle, accounts for any tilt while filming.
    slide_time = ((end_fn - start_fn)+1) * dt;

    %append data to the file. the file must have existing data in order for
    %this to work though.

    currData_new = [max_long_spread, max_lat_spread, time_long, time_lat, max_long_spread_D0, max_lat_spread_D0, slide_dist, slide_time, Test_num, Surface_type, Temp, angle, WE, test_Date];
%     currData = [max_long_spread, max_lat_spread, time_long, time_lat, Test_num, Surface_type, Temp, angle, WE, test_Date];

%     dlmwrite(datafile,currData,'-append','delimiter',',','precision','%10f')
    dlmwrite(datafile2,currData_new,'-append','delimiter',',','precision','%10f');
%     A = dlmread(datafile,',',1,0);
%     T = [A;currData];
%     dlmwrite(datafile,T,1,0);

end

disp("Tracking Complete")
if(save_tracking)
    disp("Images Saved")
end

if(save_data)
    disp("Data Saved to File")
end

function pos = Find_Edge(image)
    %The assumption is that image is already binary, but maybe that isn't
    %what we want.

    [width,length] = size(image);

    droplet_hit = false;

    largest_x = 0;
    smallest_x = width;
    largest_y = 0;
    smallest_y = length;

    top = zeros(1,2);
    bottom = zeros(1,2);
    left = zeros(1,2);
    right = zeros(1,2);
        
    for i = 1:width
        for j = 1:length
            val = image(i,j);
            if(val > 0)
                if(smallest_x > i)
                    smallest_x = i;
                    left = [smallest_x, j];
                end
                if(largest_x < i)
                    largest_x = i;
                    right = [largest_x,j]; 
                end
                if(smallest_y > j)
                    smallest_y = j;
                    top = [i,smallest_y];
                end
                if(largest_y < j)
                    largest_y  = j;
                    bottom = [i,largest_y];
                end
            end 
        end
    end

    pos = [top;right;bottom;left];
end

function plotEllipse(minorAxes,majorAxes,centers,color,width)
    theta = linspace(0,2*pi,100);
    for i=1:length(majorAxes)
        a = majorAxes(i)/2;
        b = minorAxes(i)/2;
        B = a/b;
        r = a./sqrt(cos(theta).^2 + B^2*sin(theta).^2);
        x = r.*cos(theta)+centers(i,1);
        y = r.*sin(theta)+centers(i,2);
        plot(x,y,'-','LineWidth',width,'Color',color)
    end
    set(gca,'XTick',[], 'YTick', [])
end
