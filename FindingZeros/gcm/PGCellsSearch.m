function [G, Ng] = PGCellsSearch(N, GCM)
% Identifiy persistent cells and transient cells. The codes
% in this folder are based on section 11.4 and 11.5
%
Nc = prod(N);
G = zeros(Nc+1,1); % GCM group No.
% Ng = []; % Memeber No. of each GCM group
%
SCM = SelectedSCM(N, GCM);
[Gr, ~, St] = search(SCM, N);
Epc = find(St==0); % periodic cells of SCM
%
g = 1; % persistent group counter
%
% Sink absorbing cells
SinkCells = find(Gr==1 & St==0);
G(Gr==1 & St==0) = g;
G(Nc+1) = g;
%
% Other absorbing cells
AbsorbCells = [];
NonSinkEpc = setdiff(Epc,SinkCells);
for i = 1:length(NonSinkEpc)
    cell = NonSinkEpc(i);
    if GCM{cell,1} == 1 && GCM{cell,2} == cell % Eqn (11.4.1)
        g = g+1;
        G(cell) = g;
        AbsorbCells = [AbsorbCells; cell];
    end
end
%
AbsorbCells = [AbsorbCells;SinkCells];
%
% Transient cells lead to absorbing cells (Sec. 11.5.2)
for i = 1:length(AbsorbCells)
    zprime = (Gr == Gr(AbsorbCells(i)) & St>0);
    G(zprime) = -2;
end
%
% Generate processing sequence in Epc
EpcVirgin = setdiff(Epc, AbsorbCells);
% PotentialPG = [];
%
for i = 1:length(EpcVirgin)
    %
    cellStart = EpcVirgin(i);
    flag = true;
    cellOld = cellStart;
    %
    M = 1; % Shooting sequence index
    ShootSeq = [];
    %
    while flag % Control iterations of shooting
        ImgCellNo = 1; % No. of image cell use in selected SCM
        ShootSeq = [ShootSeq; cellOld]; % store shooting sequence
        %
        % Various possibilities of cellEvent/C(B(M),1) discussed in 11.5.3
        switch G(ShootSeq(end))
            case 0 % normal
                cellOld = ShootSeq(end);
                G(cellOld) = -1; % processing
                C = GCM{cellOld,2};
                cellNew = C(ImgCellNo);
                cellOld = cellNew;
                M = M+1;
                continue
            case -1 % case (B)
                BM = ShootSeq(end-1);
                I = GCM{BM,1}; % No. of image cells of BM
                C = GCM{BM,2}; % Image cell array of BM
                %
                if I > 1 && ~all(G(C) == -1) % case (B-I)
                    %
                    ImgCellNo = 2; % use the second image cell of BM
                    opt = true;
                    %
                    while opt % repeat checking subcases in (B-I)
                        ShootSeq(end) = C(ImgCellNo);
                        %
                        switch G(ShootSeq(end))
                            case 0 % case (B-I-2)
                                G(ShootSeq(end)) = -1; % processing modified image cell
                                cellOld = ShootSeq(end); % continue the sequence
                                break
                            case -1 % case (B-I-3)
                                if length(C) >= ImgCellNo
                                    ImgCellNo = ImgCellNo + 1; % rise up to the next image cell
                                    ShootSeq(end) = C(ImgCellNo);
                                    continue % go on checking three possiblities in case (B-I)
                                end
                            otherwise % case (B-I-1)
                                % apply case (A) and stop the shooting
                                G(ShootSeq(1:end-1)) = -2;
                                flag = false;
                                break
                        end
                    end
                    %
                elseif all(G(C) == -1) % case (B-II)
                    PotentialPG = [BM; C]; % potential PG at this time
                    %
                    % Testing for a potential PG as described in
                    % 11.5.4, four cases are considered
                    G(PotentialPG) = -3;
                    j = 1;
                    check = true;
                    while ismember(-3,G) && check % Loops over potential PG to check persistent cells
                        Testcell = PotentialPG(j);
                        CTestcell = GCM{Testcell,2};
                        survive = true;
                        %
                        % check every image cell
                        for k = 1:length(CTestcell)
                            switch G(CTestcell(k))
                                case -3 % candidate cell
                                    continue
                                case -1 % under processing, but not in the candidate set
                                    G(CTestcell(k)) = -3;
                                    PotentialPG = [PotentialPG; CTestcell(k)];
                                    continue
                                case 0 % virgin cell
                                    cellOld = CTestcell(k);
                                    G(PotentialPG) = -1; % change back as under processing
                                    %
                                    % terminate testing processing and continue the sequence
                                    check = false;
                                    survive = false;
                                    break
                                otherwise % case (A) prevails, terminate current sequence
                                    G([ShootSeq(1:end-1);Testcell]) = -2;
                                    check = false;
                                    flag = false;
                                    survive = false;
                                    break
                            end
                        end
                        j = j + 1;
                    end
                    %
                    % Those who can survive the above loops are
                    % considered as a new persistent group
                    if survive
                        g = g+1;
                        G(PotentialPG) = g;
                    end
                end
                %
            otherwise % case (A)
                G(ShootSeq(1:end-1)) = -2;
                flag = false;
        end
    end
    %
end
%
% Assign member number array
G = G(1:Nc);
Gc = unique(G);
Ng = zeros(length(Gc),1);
for i = 1:length(Ng)
    members = find(Gc(i)==G);
    Ng(i) = length(members);
end