function [X,V]=histogram(A,ab,n)
% Written by Michael Wester and Stanly Steinberg in 2008.
% Sort A into bins whose edges are given by ab(1)+i*(ab(2)-ab(1))/n, i=0:n.
L = ab(2)-ab(1);
X = linspace(L/(2*n),1-L/(2*n),n);
V=zeros(1,n);
for i=1:length(A)
	k=floor(n*(A(i)-ab(1))/L)+1;
	V(k)=V(k) + 1;
end

end
