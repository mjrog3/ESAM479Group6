function data = coupled_system(init_cond, N, rx, ry, betaxy, betayx)
% general function for coupled system, based on Sugihara equation (1)
% for asymmetric case: 
% rx = 3.7; ry = 3.7; betayx = .32; betaxy=0
% init_cond = (0.2,0.4);

data = NaN(N,2);
data(1,:) = init_cond;
for ii=1:(N-1)
    data(ii+1,1) = data(ii,1)*(rx-rx*data(ii,1)-betaxy*data(ii,2));
    data(ii+1,2) = data(ii,2)*(ry-ry*data(ii,2)-betayx*data(ii,1));
end
end
