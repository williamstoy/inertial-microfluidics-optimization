close all; clc;


% how long - select
% how wide - x axis
% what flow rate - select

% output what pressure - dual y 1
% output how good focusing? - dual y 2

path = '../';

d = dir([path, 'step_data_*.csv']);
fprintf(1, 'Found %d Files\n', length(d));

shape_db = [];
data_db = {};
data_std_db = {};
all_stat_db = [];
for i = 1:length(d)
    fprintf(1, 'Loading %s\n', d(i).name);
    [shape_db, all_stat_db, data_db, data_std_db] = display_flow_progress(path, d(i).name, shape_db, all_stat_db, data_db, data_std_db, 0);
end

for i = 1:size(shape_db, 1)
   figure(shape_db(i, 1)); 
   subplot(1,2,1); hold on;
   plot([1, 150], [0.01, 0.01], 'r--', 'HandleVisibility', 'off');
   width = shape_db(i, 2);
   height = shape_db(i,3);
   notch_height = shape_db(i,4);
   xlim([0, 100]);
   xlabel('Pressure (PSI)');
   ylabel('3 STDs Y and Z (\mum)');
   set(gca(), 'FontSize', 16);
   
   subplot(1,2,2); hold on;
   plot([1, 150], [0.01, 0.01], 'r--', 'HandleVisibility', 'off');
   %ylabel('3 STDs Y and Z (\mum)');
   xlabel('Step #');
   xlim([1,100]);
   
   position = get(gcf(), 'Position');
   position(1) = position(1) - position(3);
   position(3) = position(3) * 2;
   set(gcf(), 'Position', position);
   set(gca(), 'FontSize', 16);
   lgd = legend();
   set(lgd, 'FontSize', 9);
   
   h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', [num2str(shape_db(i,1)) '.pdf']);
end

function [shape_db, all_stat_db, data_db, data_std_db] = display_flow_progress(folder, filename, shape_db, all_stat_db, data_db, data_std_db, fileoffset)
    min_re = 7;
    max_re = ceil(80*1.33);

    opts = detectImportOptions([folder, filename]);
    opts.DataLines = [fileoffset+1, Inf];
    data = table2array(readtable([folder, filename], opts));
    if isempty(data)
        return
    end
    
    if iscell(data)
        data = str2double(data);
    end
    data(:,1:4) = data(:,1:4) * 1e6;
    
    % drop the file extension
    extension = 'csv';
    filename_without_extension = filename(1:end-(length(extension) + 1));
    fnamefrag = split(filename_without_extension, '_');
    width = str2double(fnamefrag{3}(3:end));
    height = str2double(fnamefrag{4}(3:end));
    notch_height = str2double(fnamefrag{6}(3:end));
    notch_length = str2double(fnamefrag{7}(3:end));
    q = str2double(fnamefrag{8}(2:end));
    re = str2double(fnamefrag{9}(3:end));
    
    exc = 1;
    if (length(fnamefrag) >= 10)
        exc = str2double(fnamefrag{10}(4:end));
    end
    
    shape_id = str2double([num2str(width) num2str(height) num2str(notch_height), num2str(notch_length)]);
    
    if sum(shape_db == shape_id) == 0
        shape_db(end+1, :) = [shape_id, width, height, notch_height];
    end
    
    if mod(size(data,1),100) == 0
        for k = 1:(size(data,1)/100)
            d2(:,:,k) = data((1+(k-1)*100):(1+(k)*100-1),:);
        end
        
        d2 = d2(1:end,:,:);
        
        data = mean(d2,3);
        data_std = std(d2,[],3);
        
        figure(shape_id);
        subplot(1,2,1); hold on;
        %plot(2 * 3 * sqrt(data(:,4).^2 + data(:,2).^2), 'Color', [0 1-(re-50)/(166-50) 1-(re-50)/(166-50)]);
        %title(sprintf('BY PRESSURE: Width: %d\\mum, Height: %d\\mum, H Under Notch: %d\\mum', width, height, height-notch_height));
        set(gca, 'YScale', 'log');
        ylim([1e-3, 1e2]);

        % some files could be named with implicit zexclusion zone of 1,
        % otherwise explicitly include the zexclusion zone value in the
        % filename
        d = dir(sprintf([folder 'sim_params_cw%d_ch%d_*_nh%d_*_re%d*'], width, height, notch_height, floor(re)));
        if(length(fnamefrag) >= 10)
            d = dir(sprintf([folder 'sim_params_cw%d_ch%d_*_nh%d_*_re%d_exc%d*'], width, height, notch_height, floor(re), exc));
        end

        fid = fopen([d(1).folder '\' d(1).name], 'rt');
        C = textscan(fid, '%s %f %s', 'HeaderLines', 3, 'CollectOutput', true);
        pressure_drop = C{2}(1);
        focus_quality = 2* 3*sqrt(data(:,4).^2 + data(:,2).^2); %3 standard deviations (times 2 to get both positive and negative values)
        focus_quality_std = 2*3*sqrt(data_std(:,4).^2 + data_std(:,2).^2);
        c = 1-(re-min_re)/(max_re-min_re);
        color = [0 c c];
        xdata = (1:length(focus_quality))*pressure_drop;
        fill([xdata'; flipud(xdata')], [focus_quality - focus_quality_std; flipud(focus_quality + focus_quality_std)], color, 'FaceAlpha', 0.4, 'LineStyle', 'none', 'HandleVisibility', 'off');
        plot3(xdata, focus_quality, 1:length(focus_quality), 'Color', color, 'LineWidth', 3);
        fclose(fid);

        %lgd = legend();
        %lgd.String{end} = sprintf('NL=%.0f, Re=%.0f, V_{mean}=%.2f[m/s]', notch_length, re, data(end,5));

        % PLOT FOCUS VS STEP NUMBER
        subplot(1,2,2); hold on;
        set(gca, 'YScale', 'log');
        ylim([1e-3, 1e2]);

        xdata = (1:length(focus_quality));
        plot3(xdata, focus_quality, (1:length(focus_quality))*pressure_drop, 'Color', color, 'LineWidth', 3);
        fill([xdata'; flipud(xdata')], [focus_quality - focus_quality_std; flipud(focus_quality + focus_quality_std)], color, 'FaceAlpha', 0.4, 'LineStyle', 'none', 'HandleVisibility', 'off');
        lgd = legend('Units', 'normalized', 'Position', [0.75, 0.7, 0.2, 0.2]);
        lgd.String{end} = sprintf('Q=%.0f[ul/min], Re=%.0f, V_{mean}=%.2f[m/s]', q, re, data(end,5));
        
        sgtitle(sprintf('Width: %d\\mum, Height: %d\\mum, H Under Notch: %d\\mum', width, height, height-notch_height));

        %set(gcf(), 'FontSize', 20);
        
        data_db{end+1} = data;
        data_std_db{end+1} = data_std;
        all_stat_db(end+1, :) = [width, height, notch_height, notch_length, q, re, pressure_drop];
    end
    
    
    

%     figure; subplot(1,7,1:3); hold on;
%     for i = 1:size(data,1)
%         [X,Y] = calculateEllipse(data(i,1), data(i,3), data(i,2), data(i,4), 0, 1000);
%         plot(Y,X, 'Color', [0 1-i/size(data,1) 1-i/size(data,1)]);
%     end
%     plot([-width/2 width/2], [height/2-notch_height height/2-notch_height], 'r--');
%     axis equal;
%     xlim([-width/2 width/2]);
%     ylim([-height/2 height/2]);
%     
%     subplot(1,7,4);
%     axis off;
%     title(sprintf('Width: %d\\mum, Height: %d\\mum, H Under Notch: %d\\mum, Re: %.1f', width, height, height-notch_height, re));
% 
%     subplot(1,7,5:7); hold on;
%     plot(2 * 3 * data(:,2));
%     plot(2 * 3 * data(:,4));
%     xlabel('Step #');
%     ylabel('3 STD z (\mum)');
%     ylim([0 10]);
%     axis square
end

function [X,Y] = calculateEllipse(x, y, a, b, angle, steps)
    %# This functions returns points to draw an ellipse
    %#
    %#  @param x     X coordinate
    %#  @param y     Y coordinate
    %#  @param a     Semimajor axis
    %#  @param b     Semiminor axis
    %#  @param angle Angle of the ellipse (in degrees)
    %#

    narginchk(5, 6);
    if nargin<6, steps = 36; end

    beta = -angle * (pi / 180);
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha = linspace(0, 360, steps)' .* (pi / 180);
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
    Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);

    if nargout==1, X = [X Y]; end
end