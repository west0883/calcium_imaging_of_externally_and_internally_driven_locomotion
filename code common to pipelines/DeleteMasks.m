% DeleteMasks.m
% Sarah West
% 10/7/21

% Plots all the masks you drew with ManualMasking.m and lets you pick some
% for deletion.

function [masks, indices_of_mask] = DeleteMasks(masks, indices_of_mask, bRep)
% Inputs:
% masks -- a 3D matrix (pixels, pixels, mask #) of all the masks you drew.
% bRep -- a representative image (pixels, pixels) that you used to draw 

     % Get a colormap for qualitative data.
     mymap=[1 1 1 ; cbrewer('qual', 'Paired', size(masks,3), 'linear')];

     % Make blank images for holding masks 
     holder=zeros(size(masks,1), size(masks,2));
     all_masks=zeros(size(masks,1), size(masks,2));
     
     % Add each mask in its own color to  holder, and also put all masks into single image.
     for i=1:size(masks,3)
        holder(find(masks(:,:,i)))=i; 
        all_masks(find(masks(:,:,i)))=1;
     end
     
     % Put all masks onto bRep, keep original bRep for later.
     bRep_masked = bRep;
     bRep_masked(find(all_masks))=NaN; 
     
     % Create figure.
     figure; 
     
     % Plot masked bRep in first subplot
     subplot(1,2,1); imagesc(bRep_masked); colormap(gca,[0 0 0; parula(1000)]); 
     title('brain with masks');
     
     % Plot the mask overlays with their own colorscheme in second subplot.
     subplot(1,2,2); imagesc(holder); colormap(gca, mymap); colorbar; 
     title('mask numbers for removal');
     
     % Ask user if they want to delete a mask
     user_answer1= inputdlg('Do you want to delete one of these masks? 1=Y, 0=N');
     
     %Convert the user's answer into a value
     answer1=str2num(user_answer1{1});

     % While user keeps saying yes,
     while answer1 == 1 
         
         % Ask user which mask to delete.
         user_answer2= inputdlg('Which mask would you like to delete? Enter number.');
         
         %Convert the user's answer into a value
         answer2=str2num(user_answer2{1});

         % Remove that mask.
         masks(:,:,answer2) = [];
         
         % **Update indices of mask.**
         indices_of_mask =[];
         
         % For each remaining mask
         for i = 1:size(masks,3)
             
             % Find the indices of the new mask
             indices_new=find(masks(:,:,i));
             indices_of_mask = [indices_of_mask; indices_new];
             
         end
         
         % **Plot remaining masks as above.**

         % Make blank images for holding masks 
         holder=zeros(size(masks,1), size(masks,2));
         all_masks=zeros(size(masks,1), size(masks,2));

         % Add each mask in its own color to  holder, and also put all masks into single image.
         for i=1:size(masks,3)
            holder(find(masks(:,:,i)))=i; 
            all_masks(find(masks(:,:,i)))=1;
         end
         
         % Put all masks onto bRep, keep original bRep for later.
         bRep_masked = bRep;
         bRep_masked(find(all_masks))=NaN; 
          
         % Plot bRep in first subplot
         subplot(1,2,1); imagesc(bRep_masked); colormap(gca,[0 0 0; parula(1000)]); 
         title('brain with masks');

         % Plot the mask overlays with their own colorscheme in second subplot.
         subplot(1,2,2); imagesc(holder); colormap(gca, mymap); colorbar; 
         title('mask numbers for removal');
         
         % Ask user if they want to delete another mask
         user_answer1= inputdlg('Do you want to delete another mask? 1=Y, 0=N');
         
         %Convert the user's answer into a value
         answer1=str2num(user_answer1{1});

     end 
end