

np = 8;
x = linspace(0,pi,np).';
y = linspace(0,pi,np).';


a_c = [0,1];
a = cos(x.*a_c);

b_c = [0,1,2];
b = cos(y.*b_c);

z = zeros(length(x), length(a_c)*length(b_c));
item = 0;
for ii = 1:size(a,2)
    for jj = 1:size(b,2)
        item = item +1;
        z(:,item) = a(:,ii).*b(:,jj);
    end
end
