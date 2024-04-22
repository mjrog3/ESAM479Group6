function data = coupled_system(init_cond, N, rx, ry, betaxy, betayx)
% general function for coupled system, based on Sugihara equation (1)
% for asymmetric case, betaxy=0
    data = NaN(N,2);
    data(1,:) = init_cond;

    for ii=1:(N-1)
        data(ii+1,:) = step(data(ii,:),rx,ry,betaxy,betayx);
    end
end

function Y_out = step(data, rx, ry, betaxy, betayx)
    X = data(1);
    Y = data(2);

    Y_out = NaN(size(data));
    Y_out(1) = X*(rx-rx*X-betaxy*Y);
    Y_out(2) = Y*(ry-ry*Y-betayx*X);
end