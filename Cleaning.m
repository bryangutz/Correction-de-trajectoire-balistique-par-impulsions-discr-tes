N = 100; % can be increased

t_table_2 = linspace(0,1,N+1).';

X_table_2 = deval(sol, t_table_2).'; % row vector

A_table_2 = zeros(n,n,N+1);

for k = 1:N+1
    xk = X_table_2(k,:).';
    A_table_2(:,:,k) = compute_A_numeric(xk, f_dyn);
end

%%

f_1 = @(t,y) A_of_t(t,t_table_2,A_table_2)*y;
g_1 = @(t,y) reshape(f_1(t,reshape(y,n,n)),n^2,1);

f_2 = @(t,y) compute_A_numeric(deval(sol,t),f_dyn)*y;
g_2 = @(t,y) reshape(f_2(t,reshape(y,n,n)),n^2,1);

soly_1 = ode45(g_1,0:1,reshape(eye(n),n^2,1));
soly_2 = ode45(g_2,0:1,reshape(eye(n),n^2,1));

y_1 = deval(soly_1,t_table_2).';
y_2 = deval(soly_2,t_table_2).';

y_1_rs = zeros(n,n,length(t_table_2));
y_2_rs = zeros(n,n,length(t_table_2));

for k=1:length(t_table_2)
    y_1_rs(:,:,k) = reshape(y_1(k,:),n,n);
    y_2_rs(:,:,k) = reshape(y_2(k,:),n,n);
end




