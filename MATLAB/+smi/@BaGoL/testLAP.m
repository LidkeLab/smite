% Test to see if LAP/Cost Matrix can work for MAP estimates


N = 100 
ChainLength = 100
ClusterSize = 100 % Units are nm 

TruePositions = ClusterSize * rand(N,2)
PositionSigmas = gamrnd(4,2, N ,1)

clear Chain 
for ii = 1:ChainLength
    Chain(ii).X = randn(N,1) .* PositionSigmas + TruePositions(:,1)
    Chain(ii).Y = randn(N,1) .* PositionSigmas + TruePositions(:,2)
end

% Check 
figure; plot([Chain(:).X]', [Chain(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([Chain(:).X],2), mean([Chain(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)


% Indexes are not mixed up yet. 

% Try LAP 

ChainUnmixed = Chain
CoordsRef = [ChainUnmixed(1).X ChainUnmixed(1).Y]
Coords2Match = [ChainUnmixed(2).X ChainUnmixed(2).Y]
CostMatrix = pdist2(CoordsRef,Coords2Match).^2
AssignIDs = smi.SPT.solveLAP(CostMatrix)

ChainUnmixedUnmixed(2).X = ChainUnmixed(2).X(AssignIDs)
ChainUnmixed(2).Y = ChainUnmixed(2).Y(AssignIDs)


ChainUnmixed = Chain
% Do a first pass with aligning to average of everything before: 
for ii = 2:ChainLength
    CoordsRef = [mean([ChainUnmixed(1:ii-1).X],2) mean([ChainUnmixed(1:ii-1).Y],2)]
    Coords2Match = [ChainUnmixed(ii).X ChainUnmixed(ii).Y]
    CostMatrix = pdist2(CoordsRef,Coords2Match).^2
    AssignIDs = smi.SPT.solveLAP(CostMatrix)
    ChainUnmixed(ii).X = ChainUnmixed(ii).X(AssignIDs)
    ChainUnmixed(ii).Y = ChainUnmixed(ii).Y(AssignIDs)
end


% Check 
figure; plot([ChainUnmixed(:).X]', [ChainUnmixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainUnmixed(:).X],2), mean([ChainUnmixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)





% Now try mixing them up, then unmixing:
ChainMixed = Chain;
for ii = 2: ChainLength
    % draw 2 ID to swap 
    ID1 = randi(N)
    ID2 = randi(N)
    IDVec = 1:N 
    IDVec(ID1) = ID2
    IDVec(ID2)  = ID1
    ChainMixed(ii).X = ChainMixed(ii).X(IDVec)
    ChainMixed(ii).Y = ChainMixed(ii).Y(IDVec)
end
% Check 
figure; plot([ChainMixed(:).X]', [ChainMixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainMixed(:).X],2), mean([ChainMixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)


% Now unmix 


ChainUnmixed = ChainMixed;
% Do a first pass with aligning to average of everything before: 
for ii = 2:ChainLength
    CoordsRef = [mean([ChainUnmixed(1:ii-1).X],2) mean([ChainUnmixed(1:ii-1).Y],2)];
    Coords2Match = [ChainUnmixed(ii).X ChainUnmixed(ii).Y];
    CostMatrix = pdist2(CoordsRef,Coords2Match).^2;
    AssignIDs = smi.SPT.solveLAP(CostMatrix);
    ChainUnmixed(ii).X = ChainUnmixed(ii).X(AssignIDs);
    ChainUnmixed(ii).Y = ChainUnmixed(ii).Y(AssignIDs);
end

% Update Pass using mean as reference (this can oscillate)
CoordsRef = [mean([ChainUnmixed(:).X],2) mean([ChainUnmixed(:).Y],2)];
for ii = 2:ChainLength
    Coords2Match = [ChainUnmixed(ii).X ChainUnmixed(ii).Y];
    CostMatrix = pdist2(CoordsRef,Coords2Match).^2;
    AssignIDs = smi.SPT.solveLAP(CostMatrix);
    ChainUnmixed(ii).X = ChainUnmixed(ii).X(AssignIDs);
    ChainUnmixed(ii).Y = ChainUnmixed(ii).Y(AssignIDs);
end

% Check 
figure; plot([ChainUnmixed(:).X]', [ChainUnmixed(:).Y]','+')
hold on 
plot(TruePositions(:,1),TruePositions(:,2),'ko','MarkerSize',16,'LineWidth',3)
MAPN =  [mean([ChainUnmixed(:).X],2), mean([ChainUnmixed(:).Y],2)];
plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
[MAPN(:,1), MAPN(:,2)]


% Try a kmeans on mixed: 

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













