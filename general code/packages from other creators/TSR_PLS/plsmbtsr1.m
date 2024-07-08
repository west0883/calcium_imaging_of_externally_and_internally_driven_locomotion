function [X,Y,mX,sX,mY,sY,It,diff,Xrec,Yrec,T,P,C,W,Wstar]=plsmbtsr1(X,Y,A)
% TSR 1 for PLS-MB (adapted from PCA-MB)
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
X=[X,Y];
p=px+py;
mis=isnan(X);
for i=n:-1:1
  r=~isnan(X(i,:));
  pat(i).O=find(r==1); % observed variables
  pat(i).M=find(r==0); % missing variables
  pat(i).nO=size(pat(i).O,2); % number of observed variables
  pat(i).nM=size(pat(i).M,2); % number of missing variables
end
mX=zeros(1,p);
[rX, cX]=find(isnan(X));
for i=1:p
    x=X(:,i);
    mX(i)=mean(x(~isnan(x)));
end
for i=1:size(rX,1)
    X(rX(i),cX(i))=mX(cX(i));
end
maxiter=5000;
conv=1.0e-10;
diff=100;
It=0;
while It < maxiter && diff > conv
  It=It+1;
  Xmis=X(mis);
  [mXini,sXini,Xc]=MD_mean_std(X);
  S=cov(Xc);
  if n>p
      [~, ~, V] = svd(Xc,0); 
  else
      [V, ~, ~] = svd(Xc',0);
  end
  V = V(:,1:A);
  for i=1:n          
    if pat(i).nM>0    
      L=V(pat(i).O,1:min(A,pat(i).nO));
      S11=S(pat(i).O,pat(i).O);
      S21=S(pat(i).M,pat(i).O);
      z1=Xc(i,pat(i).O)';
      z2=S21*L*pinv(L'*S11*L)*L'*z1;
      Xc(i,pat(i).M)=z2';
    end
  end
  X=Xc*diag(sXini)+ones(n,1)*mXini;
  d=(X(mis)-Xmis).^2;
  diff=mean(d);
end
Y=X(:,px+1:end);
X=X(:,1:px);
[mX,sX,Xc]=MD_mean_std(X);
[mY,sY,Yc]=MD_mean_std(Y);
[T,W,P,C,Wstar]=plscode(Xc,Yc,A);
Xrec=T*P'*diag(sX)+ones(n,1)*mX;
Yrec=T*C'*diag(sY)+ones(n,1)*mY;
