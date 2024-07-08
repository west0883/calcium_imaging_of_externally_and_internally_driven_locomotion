
% Inputs:
% day -- string; the "stacks" directory where all the data is saved
% Ex: day = 'Y:\Sarah\Data\Random Motorized Treadmill\012122\1088\stacks'
% roi_1_mask = 

function [ ] = Data_verification_v3_short_Automated(day, roi_1_mask, roi_2_mask)
    
    % Change Current folder display to your stacks folder
    cd(day)
    
    % Create directory of stacks (will work for both rigs)
    day_dir=[dir(fullfile(day,'/0*/*.ome.tif'));dir(fullfile(day,'/0*.tif'))];
     
    % Flag to check if this is your first stack
    first=0;
    
    % Necessary to creat ppt file
    import mlreportgen.ppt.*
    
    % Iterate through each stack in your day_dir
    for s=1:length(day_dir)
        % Necessary to creat ppt file
        if s==1
            ppt=Presentation('Data Verification 1.pptx');
        else
            ppt=Presentation(sprintf('Data Verification %d.pptx',s),sprintf('Data Verification %d',s-1));
        end
        open(ppt);
        % Increment check (only cares later if it is 1)
        first=first+1;
        
        % Generate filename of stack to open
        fileName=fullfile(day_dir(s).folder,'/',day_dir(s).name);
        
        if first == 1
            % Load in first frame to calculate low/high
            im_list=tiffreadAltered(fileName, (1200), 'ReadUnknownTags',1);
     
            image=double(im_list.data);
            low=mean(image,'all')-(std(image,0,'all')*2);
            high=mean(image,'all')+(std(image,0,'all')*2);

            % Get ROI borders
            roi_1_border=cell2mat(bwboundaries(roi_1_mask));
            roi_2_border=cell2mat(bwboundaries(roi_2_mask));
        end

        try 
        % Load in stack and time it
            tic
            im_list=tiffreadAltered(fileName, (1:2400), 'ReadUnknownTags',1);
            disp('Time to load stack')
            toc
    
            % Extract x, y, and total pixel size from data
            [yDim, xDim]=size(im_list(1).data);
            im_pixels=yDim*xDim;
    
            % Calculate the number of images from stack
            im_list_length=length(im_list);
    
            % Sets up removal of 1st 30 seconds of data (1200 frames from
            % interleaved (13200 frame) stack)
            skip=1200;    
    
            % Separate out full stack into blue and violet stacks
            tic
            im1=im_list(skip).data;
            im2=im_list(skip+1).data;
            lev1=mean(im1(roi_1_mask==1));
            lev2=mean(im2(roi_1_mask==1));
    
            if lev1~=lev2
                if lev1>lev2
                    sel470=skip:2:im_list_length;
                    sel405=skip+1:2:im_list_length;
                else
                    sel470=skip+1:2:im_list_length;
                    sel405=skip:2:im_list_length;
                end
    
                len=min(length(sel470),length(sel405));
                sel470=sel470(1:len);
                sel405=sel405(1:len);
            end
    
            bData=zeros(len,im_pixels);
            vData=zeros(len,im_pixels);
    
            for t=1:len
    
                im=double(im_list(sel470(t)).data);
                im=reshape(im,[1 im_pixels]);
                bData(t,:) = im;
    
                im=double(im_list(sel405(t)).data);
                im=reshape(im,[1 im_pixels]);
                vData(t,:) = im;
            end
            disp('Time to separate data')
            toc
    
            % Calculate average blue and violet images
            tic
            avb=mean(bData);
            avv=mean(vData);  
    
            % Hemo correction
            hemo_Data=detrend(bData)./avb - detrend(vData)./avv; 
    
            % Reshape 2D stack into 3D stack
            bData_3D=reshape(bData',xDim,yDim,len);
            vData_3D=reshape(vData',xDim,yDim,len);
            hemo_Data_3D=reshape(hemo_Data',xDim,yDim,len);
            disp('Time to hemocorrect')
            toc    
    
            % **Start extraction**
            tic

            % Reshape roi masks.
            roi_1_mask_vector = reshape(roi_1_mask, 1,  []);
            roi_2_mask_vector = reshape(roi_2_mask, 1, []);

           % Extract data with matrix multiplication.

            bROI_1 = (roi_1_mask_vector * bData)./sum(roi_1_mask_vector);
            bROI_2 = (roi_2_mask_vector * bData)./sum(roi_2_mask_vector);

            vROI_1 = (roi_1_mask_vector * vData)./sum(roi_1_mask_vector);
            vROI_2 = (roi_2_mask_vector * vData)./sum(roi_2_mask_vector);

            hROI_1 = (roi_1_mask_vector * hemo_Data)./sum(roi_1_mask_vector);
            hROI_2 = (roi_2_mask_vector * hemo_Data)./sum(roi_2_mask_vector);
            
         
            % Collect background image to check/plot later
            bCheck=bData_3D(:,:,1);
            vCheck=vData_3D(:,:,1);
    
            bCheck_all(:,:,s)=bCheck;
    
            disp('Time to extract data')
            toc
    
            % Adds new slide to ppt
            slide=add(ppt,'Two Content',s);
    
            % Inserts title of slide to fileName
            replace(slide,'Title',fileName);
    
            fig=figure;
            % Displays blue image along with ROIs
            subplot(4,2,1)
            imshow(bCheck,[low high])
            colorbar;
            title('Image 1: should be 470 nm Image')
            hold on
            plot(roi_1_border(:,2),roi_1_border(:,1),'b','LineWidth',1);
            plot(roi_2_border(:,2),roi_2_border(:,1),'b','LineWidth',1);
            hold off
    
            % Displays violet image along with ROIs
            subplot(4,2,2)
            imshow(vCheck,[low high])
            colorbar;
            title('Image 2: should be 405 nm Image')
            hold on
            plot(roi_1_border(:,2),roi_1_border(:,1),'m','LineWidth',1);
            plot(roi_2_border(:,2),roi_2_border(:,1),'m','LineWidth',1);
            hold off
    
            % Displays raw blue and violet fluorescence data
            subplot(4,2,[3,4])
            time=0:0.05:(size(vData_3D,3)-1)*0.05;
            plot(time,(bROI_1*100),'b','LineWidth',1);
            hold on
            plot(time,(vROI_1*100),'m','LineWidth',1);
            ylabel('Raw Fluorescence (AU)')
            xlabel('Time (s)')
            legend('470nm','405nm','Location','northeastoutside');
            title('Signal from Brain ROI')
    
            % Displays max blue image to check for tail flicks/flashes
            subplot(4,2,5)
            max_blue_image=max(bData_3D,[],3);
            imshow(max_blue_image,[low high])
            colorbar;
            title('Max 470 Image')
    
            % Displays max violet image to check for tail flicks/flashes
            subplot(4,2,6)
            max_violet_image=max(vData_3D,[],3);
            imshow(max_violet_image,[low high])
            colorbar;
            title('Max 405 Image')
    
            % Displays hemo-corrected fluorescence data
            subplot(4,2,[7,8])
            time=0:0.05:(size(hemo_Data_3D,3)-1)*0.05;
            plot(time,(hROI_1*100),'b','LineWidth',1);
            ylabel('GCaMP %\DeltaF/F_{0}')
            xlabel('Time (s)')
            title('Hemo Corrected Signal from ROI')
    
            % Pause for 1 sec to allow figure to fully populate
            pause(1);
    
            % Save figure
            image=sprintf('Stack %d.png',s);
            saveas(fig,image);
    
            % Pause for 1 sec to allow figure to fully save
            pause(1);
    
            % Insert picture into left column content area
            replace(slide,'Left Content',Picture(image));
            close(gcf)
    
            % Calculate mean fluorescence levels inside both ROIs for both blue and
            % violet images

            % [Why only the first 600 frames? You're not plotting this, and you need to know the value across the whole stack, right?]
            b_brain_mean=round(mean(bROI_1(1:600)));
    
            v_brain_mean=round(mean(vROI_1(1:600)));
    
            b_out_mean=round(mean(bROI_2));
    
            v_out_mean=round(mean(vROI_2));
    
            % Perform check to determine what condition the data is in
            if b_brain_mean>v_brain_mean && b_out_mean<v_out_mean 
                check='Correct Frame Identification';
            elseif b_brain_mean<v_brain_mean && b_out_mean>v_out_mean 
                check='Incorrect Frame Identification';
            elseif b_brain_mean>v_brain_mean && b_out_mean>v_out_mean 
                check='470 nm LED left on OR 405 brighter than 470';
            elseif b_brain_mean<v_brain_mean && b_out_mean<v_out_mean
                check='405 brighter than 470';
            end
    
            % Replaces right content with outputs of the checking step
            rights = find(slide, 'Right Content');
            replace(rights, {...
                sprintf('Stack Length = %d',len),...
                sprintf('470 Brain Mean = %d',b_brain_mean),...
                sprintf('470 Outside Mean = %d',b_out_mean),...
                sprintf('405 Brain Mean = %d',v_brain_mean),...
                sprintf('405 Outside Mean = %d',v_out_mean),...
                check});
    
            % Pause to allow image to fully insert into slide
            pause(1)
        catch e
            % Adds new slide to ppt
            slide=add(ppt,'Title and Content',s);
    
            % Inserts title of slide to fileName
            replace(slide,'Title',fileName);
            
            % Ineserts description of error
            content = find(slide, 'Content');
            replace(content, {...
                sprintf('The identifier was: %s',e.identifier),...
                sprintf('There was an error! The message was: %s',e.message),...
                });
            
        end   
    % Close the ppt-->This also fully saves it
    close(ppt);
    end
    
    ppt=Presentation('Final Data Verification.pptx',sprintf('Data Verification %d',s));
       
    open(ppt);
    
    % Pause to allow image to fully insert into slide
    pause(1)
    
    % Add new slide at the end
    slide=add(ppt,'Title and Content',s+1);
    
    % Insert title of slide
    replace(slide,'Title','Alignment check between stacks');
    
    % Create new figure compairing the alingment of each stack to the 1st stack
    fig2=figure;
    fig2.WindowState='maximized';
    for b=1:size(bCheck_all,3)
        if size(bCheck_all,3)<7 
            if sum(sum(bCheck_all(:,:,b)))~=0
                subplot(2,3,b)
                first_image=bCheck_all(:,:,1);
                second_image=bCheck_all(:,:,b);
                comparison=second_image-first_image;
                imshow(comparison,[]);
                title(sprintf('Stack %d to 1st Stack',b))
            else
                subplot(2,3,b)
                box off
                axis off
                title(sprintf('Stack %d not loaded properly',b))
            end
        elseif size(bCheck_all,3)>6 && size(bCheck_all,3)<11
            if sum(sum(bCheck_all(:,:,b)))~=0
                subplot(2,5,b)
                first_image=bCheck_all(:,:,1);
                second_image=bCheck_all(:,:,b);
                comparison=second_image-first_image;
                imshow(comparison,[]);
                title(sprintf('Stack %d to 1st Stack',b))
            else
                subplot(2,5,b)
                box off
                axis off
                title(sprintf('Stack %d not loaded properly',b))
            end
        elseif size(bCheck_all,3)>10 && size(bCheck_all,3)<19
            if sum(sum(bCheck_all(:,:,b)))~=0
                subplot(3,6,b)
                first_image=bCheck_all(:,:,1);
                second_image=bCheck_all(:,:,b);
                comparison=second_image-first_image;
                imshow(comparison,[]);
                title(sprintf('Stack %d to 1st Stack',b))
             else
                subplot(3,6,b)
                box off
                axis off
                title(sprintf('Stack %d not loaded properly',b))
            end
        elseif size(bCheck_all,3)>18 && size(bCheck_all,3)<26
            if sum(sum(bCheck_all(:,:,b)))~=0
                subplot(5,5,b)
                first_image=bCheck_all(:,:,1);
                second_image=bCheck_all(:,:,b);
                comparison=second_image-first_image;
                imshow(comparison,[]);
                title(sprintf('Stack %d to 1st Stack',b)) 
             else
                subplot(5,5,b)
                box off
                axis off
                title(sprintf('Stack %d not loaded properly',b))
            end
        end
     
    end
    
    % Save new figure
    image=sprintf('Alignment check.png');
    saveas(fig2,image);
    
    % Pause to allow for saving
    pause(1)
    
    % Insert new figure into content area of slide
    replace(slide,'Content',Picture(image)); 
    
    % Pause to allow image to fully insert into slide
    pause(1)
    close(gcf)
    
    % Close the ppt-->This also fully saves it
    close(ppt);
    
    % Clean up steps to delete all image files and extra ppts generated above. Can comment out
    % if you want to keep said files.
    image_dir=dir(fullfile(cd,'/Stack *.png'));
    
    for p=1:length(image_dir)
        delete(fullfile(image_dir(p).folder,'/',image_dir(p).name));
    end
    delete(fullfile(cd,'/Alignment check.png'));
    
    ppt_dir=dir(fullfile(cd,'/Data Verification *.pptx'));
    
    for p=1:length(ppt_dir)
        delete(fullfile(ppt_dir(p).folder,'/',ppt_dir(p).name));
    end

end