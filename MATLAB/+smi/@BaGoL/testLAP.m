% Test to see if LAP/Cost Matrix can work for MAP estimates


N = 100
ChainLength = 100;
ClusterSize = 200 % Units are nm 

TruePositions = ClusterSize * rand(N,2);
PositionSigmas = gamrnd(4,1.5, N ,1);

clear Chain 
for ii = 1:ChainLength
    Chain(ii).X = randn(N,1) .* PositionSigmas + TruePositions(:,1);
    Chain(ii).Y = randn(N,1) .* PositionSigmas + TruePositions(:,2);
end

% Plot True
figure; plot([Chain(:).X]', [Chain(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([Chain(:).X],2), mean([Chain(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
title("True Mapping")

% Assignment when correctly ordered
ChainUnmixed = Chain;
CostMatrix = zeros(N,N);
% Do a first pass with aligning using MSE to everything before.
for ii = 2:ChainLength
    X = [ChainUnmixed(1:ii-1).X];
    Y = [ChainUnmixed(1:ii-1).Y];
    for jj = 1:N
            CostMatrix(:,jj) = mean( (X - ChainUnmixed(ii).X(jj)).^2 ,2) + ... 
                mean( (Y - ChainUnmixed(ii).Y(jj)).^2 ,2);
    end
    [AssignIDs, Cost] = smi.SPT.solveLAP(CostMatrix);
    sum(Cost)
    AssignIDs';
    ChainUnmixed(ii).X(AssignIDs) = ChainUnmixed(ii).X;
    ChainUnmixed(ii).Y(AssignIDs) = ChainUnmixed(ii).Y;
end

figure; plot([ChainUnmixed(:).X]', [ChainUnmixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainUnmixed(:).X],2), mean([ChainUnmixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
title('Unmixed from Correct Mapping')

% Now try mixing them up
ChainMixed = Chain;
NSwaps = ceil(N/10)
for ii = 2: ChainLength
    % draw 2 ID to swap 
    for jj=1:NSwaps
    ID1 = randi(N)
    ID2 = randi(N)
    IDVec = 1:N 
    IDVec(ID1) = ID2
    IDVec(ID2)  = ID1
    ChainMixed(ii).X = ChainMixed(ii).X(IDVec)
    ChainMixed(ii).Y = ChainMixed(ii).Y(IDVec)
    end
end
% Check 
figure; plot([ChainMixed(:).X]', [ChainMixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainMixed(:).X],2), mean([ChainMixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
title('Mixed Assignments')

% Now unmix 
ChainUnmixed = ChainMixed;
CostMatrix = zeros(N,N);
% Do a first pass with aligning using MSE to everything before.
for ii = 2:ChainLength
    X = [ChainUnmixed(1:ii-1).X];
    Y = [ChainUnmixed(1:ii-1).Y];
    for jj = 1:N
            CostMatrix(:,jj) = mean( (X - ChainUnmixed(ii).X(jj)).^2 ,2) + ... 
                mean( (Y - ChainUnmixed(ii).Y(jj)).^2 ,2);
    end
    [AssignIDs, Cost] = smi.SPT.solveLAP(CostMatrix);
    sum(Cost)
    ChainUnmixed(ii).X(AssignIDs) = ChainUnmixed(ii).X;
    ChainUnmixed(ii).Y(AssignIDs) = ChainUnmixed(ii).Y;
end
figure; plot([ChainUnmixed(:).X]', [ChainUnmixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainUnmixed(:).X],2), mean([ChainUnmixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
title('Unmixed from Mixed')

% kmeans on mixed: 
A = [ChainMixed(:).X]
B = [ChainMixed(:).Y]
KIDs = kmeans([A(:) B(:)], N)
clear MAPN
for ii = 1:N 
    MAPN.X(ii) = mean(A(KIDs==ii))
    MAPN.Y(ii) = mean(B(KIDs==ii))
end
% Check 
figure; scatter(A(:),B(:),1, KIDs)
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
plot(MAPN.X, MAPN.Y, 'rx','MarkerSize',16,'LineWidth',3)
title('KMeans from Mixed')

