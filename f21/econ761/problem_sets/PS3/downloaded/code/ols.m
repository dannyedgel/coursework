%This function estimates a simple ols regression.

%Inputs:
%Y - the dependent variable
%X - exogenous variables
%fid - a dummy variable =1 to give output =0 otherwise
%Important Outputs:
%f.b - estimated coefficients
%f.seb - standard errors
% Written by Abe Dunn, 6/14/03

function f = ols(y,X,fid)

f.b=(X'*X)^(-1)*X'*y;
f.e=y-X*f.b;
f.N=size(X,1);
f.K=size(X,2);
f.se2=(f.e'*f.e)/(f.N-f.K);
f.varb=f.se2*(X'*X)^(-1);
f.seb= full(diag(f.varb).^(1/2)); % DANNY EDIT: only sqrt-ing diagonal and 
                                  % changed from sparse to full
f.R2=1-f.e'*f.e/(y'*y-f.N*mean(y)^2);
f.adjR2=1-(f.N-1)/(f.N-f.K)*(1-f.R2);
%start of a comment and end of a comment%


if fid
   fprintf(fid,'Estimates: \n');
   for i=1:f.K
   	fprintf(fid,'\n  b%d:  %4.3f',i,f.b(i)); 
      fprintf(fid,'\n   (%4.3f) \n', f.seb(i,1)); % DANNY EDIT: sig digs
   end
   fprintf('\n\n R-squared: %3.3f  \n adjR-squared: %3.3f \n\n',f.R2, f.adjR2);
end

