function data = tentMapData(numRows, steps, epsilon)
    IC = rand(numRows, 1);
    data = zeros(numRows + 1, steps);
    data(:,1) = [IC(end) ; IC]; % periodic boundary conditions
    
    for i = 2:steps
        data(2:end, i) = tentMap(data(1:end-1, i-1), data(2:end, i-1), epsilon);
        data(1,i) = data(end, i);
    end
    
    data = (data > 0.5);
end

function out = tentMap(x_alt, x_self, epsilon)
    in = epsilon * x_alt + (1-epsilon) * x_self;
    out = (in < 0.5)*2.*in + (in >= 0.5).*(2-2*in);
end
