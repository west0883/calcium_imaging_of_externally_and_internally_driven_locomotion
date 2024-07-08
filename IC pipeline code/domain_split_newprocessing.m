% domain_split_newprocessing.m
% 8/24/20
% Sarah West
% splits the color mask ICs with multiple regions into
% individual regions. Use on the SVd compressed then ICA data.

%% load data
%close all;
clear all

num_sources=50;
z_score=3.5; 
a_thr=150;
downsample=0;
mice={'13'}; %'14';'15'; '20'; '21'; '22'; '23'; '24'}; 
dir_in=['Y:\Sarah\Analysis\motorized treadmill experiments\3 days only ICA test\consecutive days\thresholded ICs\'];
dir_out=['Y:\Sarah\Analysis\motorized treadmill experiments\3 days only ICA test\consecutive days\IC split domains\'];                 
mkdir(dir_out); 
for mousei=1:size(mice,1) 
    mouse=mice{mousei};
        
        ICs=[];       % put all spatial ICs into one 3d matrix (pixels by pixels by number of ICs) 
        if downsample==1
            load([dir_in 'm' mouse '_ICs_downsampled_' num2str(num_sources) 'sources_' num2str(z_score) 'z score thresh_' num2str(a_thr) '_pixels.mat'])
        elseif downsample==0
            load([dir_in 'm' mouse '_ICs_' num2str(num_sources) 'sources_' num2str(z_score) 'z score thresh_' num2str(a_thr) '_pixels.mat'])
        end   
        pixels=size(color_mask_pos,1);
            color_mask_pos_extend=[];
            for j=1:size(dom_cat,3)                %size(color_mask_pos,3)

              mx=max(max(dom_cat(:,:,j)));
             
                 for jj=1:mx
                     zeros_hold=zeros(pixels, pixels);
                     [r,c]=find(dom_cat(:,:,j)==jj);
                     zeros_hold(r,c)=color_mask_pos(r,c,j); 
                     color_mask_pos_extend=cat(3,color_mask_pos_extend, zeros_hold);

                 end

            end
            figure; for i=1:size(color_mask_pos_extend,3); subplot(7,8,i); imagesc(color_mask_pos_extend(:,:,i)); end
            if downsample==1
                save([dir_out 'm' mouse '_split_domains_downsampled_' num2str(num_sources) '_sources_'  num2str(z_score) 'zscore_' num2str(a_thr) 'pixels.mat'], 'color_mask_pos_extend'); 
            elseif downsample==0
                save([dir_out 'm' mouse '_split_domains_' num2str(num_sources) '_sources_'  num2str(z_score) 'zscore_' num2str(a_thr) 'pixels.mat'], 'color_mask_pos_extend'); 
                savefig([dir_out 'm' mouse '_split_domains_' num2str(num_sources) '_sources_'  num2str(z_score) 'zscore_' num2str(a_thr) 'pixels.fig']);
                 overlay=zeros(size(color_mask_pos,1), size(color_mask_pos,2));
                for j=1:size(color_mask_pos_extend,3)
                   ind2=find(color_mask_pos_extend(:,:,j)>0); 
                   overlay(ind2)=j;
                end 
                 figure; imagesc(overlay);
                 savefig([dir_out 'm' mouse '_overlay_split_domains_' num2str(num_sources) '_sources_'  num2str(z_score) 'zscore_' num2str(a_thr) 'pixels.fig']);
            end
    end
            
            
