function X = TDMA(A,B)
%
%  Function to solve tridiagonal matrix
%  A*X = B
%  A - Co-efficient square matrix / INPUT
%  X - Unknown column matrix  / FINAL OUTPUT
%  B - Constants column matrix / INPUT
%
m = size(A);
if length(m) <= 1
    error('The input should be a square matrix\n');
elseif m(1)~= m(2)
    error('The input should be a square matrix\n');
end    
%
X = zeros(length(B),1); %preaalocating array for X matrix
c = B;
n = length(B);
% getting the elements
% d - diagonal elements $ b - below diagonal $ a - above diagonal
d = diag(A);
b = zeros(n,1);
b(2:end) = diag(A,-1);
a = zeros(n,1);
a(1:end-1) = diag(A,+1);
% Forward Elimination
dummy = zeros(n,1); % Temperarory Variable                                                       
dummy(2:end) = b(2:end)./d(1:end-1);
d(2:end) = d(2:end) - dummy(2:end).*a(1:end-1);
c(2:end) = c(2:end) - dummy(2:end).*c(1:end-1);
% Backward Substitution
% Getting the last element
X(end) = c(end)/d(end);
% for remaining elements
for i = n-1:-1:1
    X(i) = (c(i) - a(i)*X(i+1))/d(i);
end
end