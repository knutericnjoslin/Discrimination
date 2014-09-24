
n=1:1:40;
m=1:1:40;
homogenous_wage = zeros(40,40);
for i=1:1:40
homogenous_wage(i,:) = n(i)./(m.*(m./(m-1)).^n(i) -m -(n(i)./(m-1)));
end;

homogenous_wage_profit = zeros(40,40);
for i=1:1:40
homogenous_wage_profit(i,:) = homogenous_wage(i,:).*(1 - (1 - 1./m).^n);
end;
