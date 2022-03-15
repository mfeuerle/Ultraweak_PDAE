function [ M ] = assemble_matrix(MS, MT)

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

M.inner_product =@(x1,x2) eval_bothsides(x1',x2,MS,MT);
M.norm  =@(x) sqrt(eval_bothsides(x',x,MS,MT));
M.eval  =@(x) eval_matvec(x,MS,MT);
M.evalT =@(x) eval_transposed(x,MS,MT);
M.eval2 =@(x1,x2) eval_bothsides(x1,x2,MS,MT);
M.full  =@() full_matrix(MS,MT);

end

function Bx = eval_matvec(x,MS,MT)
[K2, K1] = size(MT{1,1}{1});
n = size(MS{1,1}{1},1);
d = size(MS{2,2}{1},1);

Bx = zeros(n*K2+d,1);

xk = cell(2,1);
Bxk = cell(2,1);
for k = 1:size(x,2)
    xk{1} = reshape(x(1:n*K1,k), [n,K1]);
    xk{2} = x(n*K1+1:end,k);
    Bxk{1} = zeros(n,K2);
    Bxk{2} = zeros(d,1);
    
    for i = 1:2
        for j = 1:2
            for l = 1:length(MS{i,j})
                Bxk{i} = Bxk{i} + MS{i,j}{l} * xk{j} * MT{i,j}{l}';
            end
        end
    end
    Bx(:,k) = [reshape(Bxk{1}, [n*K2, 1]); Bxk{2}];
end
end

function Bx = eval_transposed(x,MS,MT)
K = size(MT{1,1}{1},1);
n = size(MS{1,1}{1},1);
d = size(MS{2,2}{1},1);

Bx = zeros(size(x));

xk = cell(2,1);
Bxk = cell(2,1);
for k = 1:size(x,2)
    xk{1} = reshape(x(1:n*K,k), [n,K]);
    xk{2} = x(n*K+1:end,k);
    Bxk{1} = zeros(n,K);
    Bxk{2} = zeros(d,1);
    
    for i = 1:2
        for j = 1:2
            for l = 1:length(MS{i,j})
                Bxk{i} = Bxk{i} + MS{j,i}{l}' * xk{j} * MT{j,i}{l};
            end
        end
    end
    Bx(:,k) = [reshape(Bxk{1}, [n*K, 1]); Bxk{2}];
end
end

function x1Bx2 = eval_bothsides(x1,x2,MS,MT)
K = size(MT{1,1}{1},1);
n = size(MS{1,1}{1},1);
d = size(MS{2,2}{1},1);

Bx2 = zeros(size(x2));

xk = cell(2,1);
Bxk = cell(2,1);
for k = 1:size(x2,2)
    xk{1} = reshape(x2(1:n*K,k), [n,K]);
    xk{2} = x2(n*K+1:end,k);
    Bxk{1} = zeros(n,K);
    Bxk{2} = zeros(d,1);
    
    for i = 1:2
        for j = 1:2
            for l = 1:length(MS{i,j})                
                Bxk{i} = Bxk{i} + MS{i,j}{l} * xk{j} * MT{i,j}{l}';
            end
        end
    end
    Bx2(:,k) = [reshape(Bxk{1}, [n*K, 1]); Bxk{2}];
end
x1Bx2 = x1*Bx2;
end

function B_full = full_matrix(MS,MT)
B_full = cell(2,2);
for i = 1:2
    for j = 1:2
        B_full{i,j} = kron(MT{i,j}{1},MS{i,j}{1});
        for l = 2:length(MS{i,j})
            B_full{i,j} = B_full{i,j} + kron(MT{i,j}{l},MS{i,j}{l});
        end
    end
end
B_full = [B_full{1,1}, B_full{1,2};
          B_full{2,1}, B_full{2,2}];
end