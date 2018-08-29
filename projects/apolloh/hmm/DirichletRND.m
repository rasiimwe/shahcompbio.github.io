function X = DirichletRND(alpha)

n = length(alpha);
Y = zeros(1,n);
for i = 1:n
    Y(i) = gamrnd(alpha(i),1);
end

X = Y/sum(Y);