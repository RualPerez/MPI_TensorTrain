% Matlab test script: Runs the check.m which is equivalent to our
% algorithm1
% Remember to clone the following repository before running it:
% https://github.com/kbatseli/TTClassifier

d = 10;
general_poly_order = 2;
general_rank = 3;


Input = ones(1,d);
x=cell(1,d);
n = general_poly_order *ones(1,d);
r = general_rank *ones(1,d+1); r(1) = 1;    r(d+1)=1;
for i=1:d
    x{i}=ones(r(i),n(i),r(i+1));
end
check(Input,x)