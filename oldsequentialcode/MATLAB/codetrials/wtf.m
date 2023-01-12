clear all
t = []
p1 = [1,2,3,4,5]
p2 = [6,7,8,9,10]
for i = 1:length(p1)
exp_vector = p1(i)*ones(size(p2))';
t = [t;exp_vector,p2']; %#ok<*AGROW>
end
tt = t(:,1).*t(:,2);
pl = reshape(tt,length(p2),[])'

%%
[X,Y,Z] = peaks(5)
surf(X,Y,Z)

