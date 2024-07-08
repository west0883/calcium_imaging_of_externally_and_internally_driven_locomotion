function [X,Y,mX,sX,mY,sY,It,diff,Xrec,Yrec,T,P,C,W,Wstar]=plsmbtsr2(X,Y,A)
% TSR 2 for PLS-MB (adapted from PLS-ME)
%
% INPUTS:
% X: predictor data matrix with NaN for the missing data
% Y: response data matrix with NaN for the missing data
% A: number of PLS components
% 
% OUTPUTS:
% X: original predictor data matrix with imputed values
% Y: original response data matrix with imputed values
% mX: estimated X mean
% mY: estimated Y mean
% sX: estimated X covariance matrix
% sY: estimated Y covariance matrix
% It: number of iterations needed
% diff: tolerance reached for convergence
% Xrec: reconstructed X using the PLS scores and loadings
% Yrec: reconstructed Y using the PLS scores and loadings
% T: scores of the final PLS model
% P: X loadings of the final PLS model
% C: Y loadings of the final PLS model
% W: weights of the final PLS model
% Wstar: normalized weights of the final PLS model 
%
% Copyright (C) 2016 A. Folch-Fortuny and F. Arteaga
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
[n,px]=size(X);
py=size(Y,2);
misX=isnan(X);
misY=isnan(Y);
for i=n:-1:1,
  r=~isnan(X(i,:));
  patX(i).O=find(r==1); % observed variables
  patX(i).M=find(r==0); % missing variables
  patX(i).nO=size(patX(i).O,2); % number of observed variables
  patX(i).nM=size(patX(i).M,2); % number of missing variables
  r=~isnan(Y(i,:));
  patY(i).O=find(r==1); % observed variables
  patY(i).M=find(r==0); % missing variables
  patY(i).nO=size(patX(i).O,2); % number of observed variables
  patY(i).nM=size(patX(i).M,2); % number of missing variables
end
mX=zeros(1,px);
[rX cX]=find(isnan(X));
for i=1:px,
    x=X(:,i);
    mX(i)=mean(x(~isnan(x)));
end
for i=1:size(rX,1),
    X(rX(i),cX(i))=mX(cX(i));
end
mY=zeros(1,py);
[rY cY]=find(isnan(Y));
for i=1:py,
    y=Y(:,i);
    mY(i)=mean(y(~isnan(y)));
end
for i=1:size(rY,1),
    Y(rY(i),cY(i))=mY(cY(i));
end
diff=100;
It=0;
while It<5000 & diff>10e-10,
  It=It+1;
  Xmis=X(misX);
  Ymis=Y(misY);
  [mXini,sXini,Xc]=MD_mean_std(X);
  [mYini,sYini,Yc]=MD_mean_std(Y);
  SX=cov(Xc);
  SY=cov(Yc);
  [T,W,P,C]=plscode(Xc,Yc,A);
  H=cov(T);
  for i=1:n, 
      if patX(i).nM>0,
          S11=SX(patX(i).O,patX(i).O);
          P1=P(patX(i).O,1:min(A,patX(i).nO));
          R1=W(patX(i).O,1:min(A,patX(i).nO))*pinv(P(:,1:min(A,patX(i).nO))'*W(:,1:min(A,patX(i).nO)));
          B=pinv(R1'*S11*R1)*R1'*P1(:,1:min(A,patX(i).nO))*H(1:min(A,patX(i).nO),1:min(A,patX(i).nO));
          t=B'*R1'*Xc(i,patX(i).O)';
          Xc(i,patX(i).M)=t'*P(patX(i).M,1:min(A,patX(i).nO))';
          if patY(i).nM>0
              Yc(i,patY(i).M)=t'*C(patY(i).M,1:min(A,patX(i).nO))';
          end
      else
          if patY(i).nM>0
              S11=SX(patX(i).O,patX(i).O);
              P1=P(patX(i).O,1:min(A,patX(i).nO));
              R1=W(patX(i).O,1:min(A,patX(i).nO))*pinv(P(:,1:min(A,patX(i).nO))'*W(:,1:min(A,patX(i).nO)));
              B=pinv(R1'*S11*R1)*R1'*P1(:,1:min(A,patX(i).nO))*H(1:min(A,patX(i).nO),1:min(A,patX(i).nO));
              t=B'*R1'*Xc(i,patX(i).O)';
              Yc(i,patY(i).M)=t'*C(patY(i).M,:)';
          end
      end
  end
  X=Xc*diag(sXini)+ones(n,1)*mXini;
  Y=Yc*diag(sYini)+ones(n,1)*mYini;
  d=[(X(misX)-Xmis).^2;(Y(misY)-Ymis).^2];
  diff=mean(d);  
end
[mX,sX,Xc]=MD_mean_std(X);
[mY,sY,Yc]=MD_mean_std(Y);
[T,W,P,C,Wstar]=plscode(Xc,Yc,A);
Xrec=T*P'*diag(sX)+ones(n,1)*mX;
Yrec=T*C'*diag(sY)+ones(n,1)*mY;