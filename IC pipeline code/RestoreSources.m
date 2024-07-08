% RestoreSources.m
% Sarah West
% 3/30/22

% Plots all the sources you removed with the RemoveArtifacts.m code and lets
% you restore some that you might have removed by accident. 

function [parameters] = RestoreSources(parameters)
% Inputs:
% sources -- a 3D matrix (pixels, pixels, source #) of all the sources you drew.
% bRep -- a representative image (pixels, pixels) that you used to draw 
    
    sources_to_remove = parameters.sources.sources_removed; 

     % Get a colormap for qualitative data.
     mymap=[1 1 1 ; cbrewer('qual', 'Paired', size(sources_to_remove,3), 'linear')];

     % Make blank images for holding sources 
     holder=zeros(size(sources_to_remove,1), size(sources_to_remove,2));
     all_sources=zeros(size(sources_to_remove,1), size(sources_to_remove,2));
     
     % Add each source in its own color to  holder, and also put all sources into single image.
     for i=1:size(sources_to_remove,3)
        holder(find(sources_to_remove(:,:,i)))=i; 
        all_sources(find(sources_to_remove(:,:,i)))=1;
     end
     

     % Put all sources onto bRep, keep original bRep for later.
     % Make a representative image with just the brain mask, if present.
     if isfield(parameters, 'indices_of_mask')
         bRep = 
     bRep_masked = bRep;
     bRep_masked(find(all_sources))=NaN; 
     
     % Create figure.
     figure; 
     
     % Plot masked bRep in first subplot
     subplot(1,2,1); imagesc(bRep_masked); colormap(gca,[0 0 0; parula(1000)]); 
     title('brain with sources');
     
     % Plot the sources overlays with their own colorscheme in second subplot.
     subplot(1,2,2); imagesc(holder); colormap(gca, mymap); colorbar; 
     title('sources numbers for restoring');
     
     % Ask user if they want to restore a sources
     user_answer1= inputdlg('Do you want to restore one of these sourcess? 1=Y, 0=N');
     
     %Convert the user's answer into a value
     answer1=str2num(user_answer1{1});

     % While user keeps saying yes,
     while answer1 == 1 
         
         % Ask user which sources to restore.
         user_answer2= inputdlg('Which source would you like to restore? Enter number.');
         
         %Convert the user's answer into a value
         answer2=str2num(user_answer2{1});

         % Remove that sources.
         sources_to_remove(:,:,answer2) = [];
         
         % **Plot remaining sourcess as above.**

         % Make blank images for holding sourcess 
         holder=zeros(size(sourcess_to_remove,1), size(sources_to_remove,2));
         all_sources=zeros(size(sources_to_remove,1), size(sources_to_remove,2));

         % Add each source in its own color to  holder, and also put all sources into single image.
         for i=1:size(sources_to_remove,3)
            holder(find(sources_to_remove(:,:,i)))=i; 
            all_sources(find(sources_to_remove(:,:,i)))=1;
         end
         
         % Put all sources onto bRep, keep original bRep for later.
         bRep_masked = bRep;
         bRep_masked(find(all_sources))=NaN; 
          
         % Plot bRep in first subplot
         subplot(1,2,1); imagesc(bRep_masked); colormap(gca,[0 0 0; parula(1000)]); 
         title('brain with removed sources');

         % Plot the source overlays with their own colorscheme in second subplot.
         subplot(1,2,2); imagesc(holder); colormap(gca, mymap); colorbar; 
         title('source numbers for restoring');
         
         % Ask user if they want to restore another source
         user_answer1= inputdlg('Do you want to restore another source? 1=Y, 0=N');
         
         %Convert the user's answer into a value
         answer1=str2num(user_answer1{1});

     end 
end