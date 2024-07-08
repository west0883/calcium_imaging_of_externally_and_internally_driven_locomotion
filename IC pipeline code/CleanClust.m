function [CReg,count]=CleanClust(map)

% function will eliminate pixels that form excursions(attached only by 1
% side) or filaments (attached by two parallel sdes only)

[imax,jmax]=size(map);

list=zeros(imax*jmax,2);
count=0;

for i=1:imax
    for j=1:jmax
        
        if map(i,j)>0
            crt=[0 0 0 0];
            i1=i-1;
            i2=i+1;
            j1=j-1;
            j2=j+1;
            
            if i1>=1
                if map(i1,j)>0
                    crt(1)=1;
                end
            end
            
            if i2<=imax
                if map(i2,j)>0
                    crt(2)=1;
                end
            end
            
            if j1>=1
                if map(i,j1)>0
                    crt(3)=1;
                end
            end
            
            if j2<=jmax
                if map(i,j2)>0
                    crt(4)=1;
                end
            end
            
            if sum(crt)==1
                count=count+1;
                list(count,:)=[i,j];
            else
                if sum(crt(1:2))==0 || sum(crt(3:4))==0
                    count=count+1;
                    list(count,:)=[i,j];
                end
            end
        end
    end
end

list(count+1:end,:)=[];
CReg=map;
for ix=1:count
    i=list(ix,1);
    j=list(ix,2);
    CReg(i,j)=0;
end
    


                
                
                
                
            
            
        
        