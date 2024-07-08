function [Reg,FinSize,DomId] = ClustRegFill(map,thr)

map(map~=0)=1;
map=map-1;
map=abs(map);

s=size(map);
s1=s(1);
s2=s(2);
Reg=zeros(s1,s2);


ClusSize=zeros(1,700000);

if map(1,1)>0
    ClusNum=1;
    Reg(1,1)=ClusNum;
    ClusSize(ClusNum)=ClusSize(ClusNum)+1;
else
    ClusNum=0;
end

for n=2:s2
    if map(1,n)>0
        if Reg(1,n-1)>0
            Reg(1,n)=Reg(1,n-1);
            ClusSize(Reg(1,n-1))=ClusSize(Reg(1,n-1))+1;
        else
            ClusNum=ClusNum+1;
            Reg(1,n)=ClusNum;
            ClusSize(ClusNum)=ClusSize(ClusNum)+1;
        end
    end
end

for m=2:s1
    if map(m,1)>0
        if Reg(m-1,1)>0
            Reg(m,1)=Reg(m-1,1);
            ClusSize(Reg(m-1,1))=ClusSize(Reg(m-1,1))+1;
        else
            ClusNum=ClusNum+1;
            Reg(m,1)=ClusNum;
            ClusSize(ClusNum)=ClusSize(ClusNum)+1;
        end
    end
end

for m=2:s1
    for n=2:s2
        if map(m,n)>0
            if Reg(m-1,n)==Reg(m,n-1)
                if Reg(m-1,n)+Reg(m,n-1)==0
                    ClusNum=ClusNum+1;
                    Reg(m,n)=ClusNum;
                    ClusSize(ClusNum)=ClusSize(ClusNum)+1;
                else
                    Reg(m,n)=Reg(m,n-1);
                    ClusSize(Reg(m,n-1))=ClusSize(Reg(m,n-1))+1; 
                end
            else
                if Reg(m-1,n)*Reg(m,n-1)==0
                    mmm=max(Reg(m-1,n),Reg(m,n-1));
                    Reg(m,n)=mmm;
                    ClusSize(mmm)=ClusSize(mmm)+1;
                else
                    mmm=max(Reg(m-1,n),Reg(m,n-1));
                    mm=min(Reg(m-1,n),Reg(m,n-1));
                    for i=1:s1
                        for j=1:s2
                            if Reg(i,j)==mmm
                                Reg(i,j)=mm;
                            else
                                if Reg(i,j)>mmm
                                    Reg(i,j)=Reg(i,j)-1;
                                end
                            end
                        end
                    end
                    
                    ClusSize(mm)=ClusSize(mm)+ClusSize(mmm);
                    Reg(m,n)=mm;
                    ClusSize(mm)=ClusSize(mm)+1;
                    
                    for s=mmm:ClusNum
                        ClusSize(s)=ClusSize(s+1);
                    end
                    
                    ClusNum=ClusNum-1;
                end
            end
        end
    end
end
ClusNum=ClusNum;
CluRegNum=max(max(Reg));

%for i=1:s1
%    for j=1:s2
%        if Reg(i,j)==0
%            Reg(i,j)=-1000;
%        end
%    end
%end

for t=ClusNum:-1:1
    if ClusSize(t)>=thr
        for i=1:s1
            for j=1:s2
                if Reg(i,j)==t
                    Reg(i,j)=0;
                end
            end
        end
    end
end
f=unique(Reg);
DomId=f(2:length(f));
FinSize=ClusSize(DomId);

% for d=DomId
%     Regf=zeros(s1,s2);
%     Regd=Reg;
%     Regd(Regd~=d)=0;
%     
%     ClusSize=zeros(1,70000);
% 
%     if Regd(1,1)==0
%         ClusNum=1;
%         Regf(1,1)=ClusNum;
%         ClusSize(ClusNum)=ClusSize(ClusNum)+1;
%     else
%         ClusNum=0;
%     end
% 
%     for n=2:s2
%         if Regd(1,n)==0
%             if Regf(1,n-1)~=0
%                 Regf(1,n)=Regf(1,n-1);
%                 ClusSize(Regf(1,n-1))=ClusSize(Regf(1,n-1))+1;
%             else
%                 ClusNum=ClusNum+1;
%                 Regf(1,n)=ClusNum;
%                 ClusSize(ClusNum)=ClusSize(ClusNum)+1;
%             end
%         end
%     end
% 
%     for m=2:s1
%         if Regd(m,1)==0
%             if Regf(m-1,1)~=0
%                 Regf(m,1)=Regf(m-1,1);
%                 ClusSize(Regf(m-1,1))=ClusSize(Regf(m-1,1))+1;
%             else
%                 ClusNum=ClusNum+1;
%                 Regf(m,1)=ClusNum;
%                 ClusSize(ClusNum)=ClusSize(ClusNum)+1;
%             end
%         end
%     end
% 
%     for m=2:s1
%         for n=2:s2
%             if Regd(m,n)==0
%                 if Regf(m-1,n)==Regf(m,n-1)
%                     if Regf(m-1,n)+Regf(m,n-1)>0
%                         ClusNum=ClusNum+1;
%                         Regf(m,n)=ClusNum;
%                         ClusSize(ClusNum)=ClusSize(ClusNum)+1;
%                     else
%                         Regf(m,n)=Regf(m,n-1);
%                         ClusSize(Regf(m,n-1))=ClusSize(Regf(m,n-1))+1; 
%                     end
%                 else
%                     if Regf(m-1,n)*Regf(m,n-1)~=0
%                         mmm=max(Regf(m-1,n),Regf(m,n-1));
%                         Regf(m,n)=mmm;
%                         ClusSize(mmm)=ClusSize(mmm)+1;
%                     else
%                         mmm=max(Regf(m-1,n),Regf(m,n-1));
%                         mm=min(Regf(m-1,n),Regf(m,n-1));
%                         for i=1:s1
%                             for j=1:s2
%                                 if Regf(i,j)==mmm
%                                     Regf(i,j)=mm;
%                                 else
%                                     if Regf(i,j)>mmm
%                                         Regf(i,j)=Reg(i,j)-1;
%                                     end
%                                 end
%                             end
%                         end
% 
%                         ClusSize(mm)=ClusSize(mm)+ClusSize(mmm);
%                         Reg(m,n)=mm;
%                         ClusSize(mm)=ClusSize(mm)+1;
% 
%                         for s=mmm:ClusNum
%                             ClusSize(s)=ClusSize(s+1);
%                         end
% 
%                         ClusNum=ClusNum-1;
%                     end
%                 end
%             end
%         end
%     end
% ClusNum=ClusNum;
% CluRegNum=max(max(Reg));
% 
% %for i=1:s1
% %    for j=1:s2
% %        if Reg(i,j)==0
% %            Reg(i,j)=-1000;
% %        end
% %    end
% %end
% 
% for t=ClusNum:-1:1
%     if ClusSize(t)<=thr
%         for i=1:s1
%             for j=1:s2
%                 if Reg(i,j)==t
%                     Reg(i,j)=0;
%                 end
%             end
%         end
%     end
% end
% f=unique(Reg);
% DomId=f(2:length(f));
% FinSize=ClusSize(DomId);
%     
%     
% 
