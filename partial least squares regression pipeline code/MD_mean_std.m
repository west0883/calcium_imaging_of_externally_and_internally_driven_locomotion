function [m,s,Y]=MD_mean_std(X)
% Autoscaling function for MD
%
% INPUTS
% X: data matrix
%
% OUTPUTS
% m: mean
% s: covariance matrix
% Y: autoscaled data matrix
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
[N,K]=size(X);
m=zeros(1,K);
s=zeros(1,K);
for j=1:K,
    x=X(:,j);
    ind = ~isnan(x);
    n=sum(ind);
    m(j)=mean(x(ind));
    s(j)=sqrt((x(ind)-m(j))'*(x(ind)-m(j))/(n-1));
end
Y=X;
Y=Y-ones(N,1)*m;
for i=1:N,
    for j=1:K,
        Y(i,j)=Y(i,j)/s(j);
    end
end