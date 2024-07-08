function [T,W,P,C,Wstar]=plscode(X,Y,A)
% PLS
% 
% INPUTS
% X: predictor data matrix
% Y: response data matrix
% A: number of PLS components
% 
% OUTPUTS
% T: scores
% P: X loadings
% C: Y loadings
% W: weights
% Wstar: normalized weights
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
[N,K]=size(X);M=size(Y,2);
T=zeros(N,A);U=zeros(N,A);
W=zeros(K,A);C=zeros(M,A);
P=zeros(K,A);
nIt=zeros(1,A+1);
Tol=1e-10;
Xc=X;Yc=Y;
mX=mean(Xc);mY=mean(Yc);
E=Xc-ones(N,1)*mX;
F=Yc-ones(N,1)*mY;
for a=1:A,
  u=F(:,1);
  dif=100;
  while dif>Tol,
    nIt(a)=nIt(a)+1;
    w=E'*u/(u'*u);
    w=w/norm(w);
    t=E*w;
    c=F'*t/(t'*t);
    unew=F*c/(c'*c);
    dif=norm(unew-u);
    u=unew;
  end
  p=E'*t/(t'*t);
  E=E-t*p';F=F-t*c';
  U(:,a)=u(:,1);T(:,a)=t(:,1);
  W(:,a)=w(:,1);C(:,a)=c(:,1);
  P(:,a)=p(:,1);
end    
Wstar=W*pinv(P'*W);
