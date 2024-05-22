function Process_Cyan_WhiteBkgrd_X_Other_setup_correctWhen2
% This function processes images (.jpg) in the folder it is currently in, assuming
% a white background (3M paper) and saves the data processed in an m-file

%% Get & store files
rect = [2067.51000000000,1443.51000000000,485.980000000000,395.980000000000]; %coordinates obtained from the picture, specific to the current imaging setup

current_folder = pwd; %Get current working directory
All_jpg_files = imageDatastore(current_folder,'FileExtensions',{'.jpg'});

%% Block to start ImageJ

javaaddpath 'C:\Program Files\MATLAB\R2020a\java\jar\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2020a\java\jar\ij.jar'

MIJ.start('C:\Program Files\MATLAB\R2020a\java\jar\'); %start ImageJ
% 

%% Begin Image Processing

Cyan_values = cell(length(All_jpg_files.Files),3); 
% Cell with C1 = File name, C2 = Top circle intensity, C3 = Bottom circle
% intensity

for i = 1:length(All_jpg_files.Files)
    current_imageData = readimage(All_jpg_files,i);
    [x,y,~] = size(current_imageData); %check size, only implement crop IF size is large
    cropped = false;
    if (x == 3264) && (y == 4896)
        current_imageData = imcrop(current_imageData,rect); %crops image to predefined rectangle
        cropped = true; %To send to bottom code
    end
    
%     imshow(current_imageData)
    
    current_imageDir = All_jpg_files.Files{i}; %gets current opening file directory for ImageJ
    
    name_index_start = strfind(current_imageDir,'\'); %Find all the backslashes
    name_index_start = name_index_start(end); %get the last slash, indicative of where name is
    
    current_imageName = current_imageDir(name_index_start+1:end);
    
    Cyan_values{i,1} = current_imageName;
    
    Temp_Cyan = process_Mean_Cyan_white(current_imageData,current_imageDir,rect,cropped);
    
    % Export data
    Cyan_values{i,2} = Temp_Cyan(1);
    Cyan_values{i,3} = Temp_Cyan(2); 
    
end %Repeat for all files in folder

MIJ.exit; %close ImageJ
save('Processed_Cyan_data.mat','Cyan_values');

end

function Cyan_val = process_Mean_Cyan_white(imageData,imageDir,cropping_rect,cropped)

%This function takes the image datafile (as processed by imread/readimage)
%and processes it by first detecting the circles, then using ImageJ to
%process the circle data for cyan intensity

% Inputs: imageData = image datafile as processed by imread/readimage
%         imageDir = image directory for ImageJ to open the file  

% Output: Cyan_val = vector corresponding to mean cyan intensities of the
%                    top and bottom wells (in image)

% MIJ.setRoi(coordinates of bounding rectangle, 1) - 1 is for oval
% [Top left, Top Right, Bottom Right Bottom Left]
% For now, take a radius of 30 for the white background version

Cyan_val = zeros(1,2); %create a vector to store values
R = 22; %Circles are typically 22 pixels in radius
%% For input image, determine where the circles are

% Tune to detect circles
Binarize_threshold = 0.76; %starting value for binarization (Higher = more foreground --> circle becomes larger/darker/more obvious)
Sensitivity_threshold = 0.89; %starting value for binarization (Higher = Detect more circles)

Binarized_image = imbinarize(rgb2gray(imageData),Binarize_threshold); %Get a binary image for circle detection 
[centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);


iter_count = 0; %stopping from going on an infinite loop

if length(radii) == 2
    dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
    x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
    while ((dist_between_centers < 30) || (x_dist>40)) && (length(radii) == 2)
        if iter_count>100
            centers = [50 51; 50 100]; % Visible after visualization to be noted as an error for reprocessing
            break;
        end
        Binarize_threshold = Binarize_threshold+0.005;
        Binarized_image = imbinarize(rgb2gray(imageData),Binarize_threshold); %Get a binary image for circle detection 
        [centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);
        if length(radii) ~= 2
            break;
        end
        dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
        x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
        iter_count = iter_count+1;
    end
end

iter_count = 0; %reset
    

while length(radii)~= 2
    if iter_count>100 % avoid infinite loop
        centers = [50 51; 50 100]; % Visible after visualization to be noted as an error for reprocessing
        break;
    end
    if (length(radii) < 2) 
        Binarize_threshold = Binarize_threshold+0.005;
        Binarized_image = imbinarize(rgb2gray(imageData),Binarize_threshold); %Get a binary image for circle detection 
        [centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);
     
        if (length(radii) == 2) 
            dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
            x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
            while ((dist_between_centers < 30) || (x_dist>40)) && (length(radii) == 2)
                Binarize_threshold = Binarize_threshold+0.005;
                Binarized_image = imbinarize(rgb2gray(imageData),Binarize_threshold); %Get a binary image for circle detection 
                [centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);
                if length(radii) ~= 2
                    break;
                end
                dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
                x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
            end
        end
    elseif length(radii) > 2
        Sensitivity_threshold = Sensitivity_threshold-0.005;
        [centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);
        if (length(radii) == 2) 
            dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
            x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
            while ((dist_between_centers < 30) || (x_dist>40)) && (length(radii) == 2)
                Binarize_threshold = Binarize_threshold+0.005;
                Binarized_image = imbinarize(rgb2gray(imageData),Binarize_threshold); %Get a binary image for circle detection 
                [centers, radii] = imfindcircles(Binarized_image,[20 40],'ObjectPolarity','dark','Sensitivity',Sensitivity_threshold);
                if length(radii) ~= 2
                    break;
                end
                dist_between_centers = sqrt(sum((centers(1,:) - centers(2,:)).^2)); %calculate distance between centers found - guards against finding the same circle
                x_dist = abs(centers(1,1) - centers(2,1)); %find distance between x's
            end
        end
    end
    iter_count = iter_count+1;
end


y_pos = centers(:,2);

[~,pos_index] = sort(y_pos);

centers(:,1) = centers(pos_index,1);
centers(:,2) = centers(pos_index,2);

%% QCing images - visualize to confirm correct identification of circles

name_index_start = strfind(imageDir,'\'); %Find all the backslashes
name_index_start = name_index_start(end); %get the last slash, indicative of where name is
imageName = imageDir(name_index_start+1:end);

figure;
imshow(imageData);
h = viscircles(centers,[R R].');
title(imageName);


%% ImageJ Processing
rect_coords = [(cropping_rect(1)) (cropping_rect(1)+cropping_rect(3)) (cropping_rect(1)+cropping_rect(3)) (cropping_rect(1))
               (cropping_rect(2)) (cropping_rect(2)) (cropping_rect(2)+cropping_rect(4)) (cropping_rect(2)+cropping_rect(4))];
                 
for i = 1:2 %only 2 circles in image. Note that centers are ouputted wrt to row (1 row corresponds to 1 circle
    ij.IJ.open(imageDir);
    if cropped
        MIJ.setRoi(rect_coords,2);
        MIJ.run("Crop");
    end
    circle_coords = [(centers(i,1)-R) (centers(i,1)+R) (centers(i,1)+R) (centers(i,1)-R)
                     (centers(i,2)-R) (centers(i,2)-R) (centers(i,2)+R) (centers(i,2)+R)];
    MIJ.setRoi(circle_coords,1); 
    
    % Cyan processing
    MIJ.run("Crop");
    
    MIJ.run("RDT Cyan MIJ macro"); % Run inbuilt macro
    pause(1); %pause 2.5s for macro
    Cyan_val(i) = MIJ.getColumn("Mean");
    MIJ.run("Closing macro");
    
end
    

end