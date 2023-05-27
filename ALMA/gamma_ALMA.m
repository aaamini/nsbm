Z_all = zeros(1100, 20);
xi_all = zeros(1100, 3, 200);

for file = 1:4

    for rep = 1:275

        if rep == 1
            load("./gamma_networks" + file + ".mat")

            nreps = size(A, 1);
            n = size(A, 3);
            L = size(A, 2);

            params.K = [2, 3, 5];
            params.L = L;
            params.M = 3;
            params.n = n;

            Z = zeros(nreps, L);
            xi = zeros(nreps, params.M, n);
        end
    
        [Q1, W1] = AltMin(squeeze(A(rep,:,:,:)), params);
        labelhat = kmeans(W1, params.M, "Replicates", 100);

        Z(rep,:) = labelhat;

        label_comhat = cell(1, params.M);

        for i = 1:params.M

            Qslice = squeeze(-Q1(i,:,:));
            k = params.K(i);
            [U, S, V] = svd(Qslice);
            label_comhat{i} = kmeans(U(:,1:k), k, "Replicates", 100);

            xi(rep,i,:) = label_comhat{i};
        end

    end

    Z_all(275*(file-1)+1:275*file,:) = Z;
    xi_all(275*(file-1)+1:275*file,:,:) = xi;

end

save("./gamma_ALMA.mat", "Z", "xi")
