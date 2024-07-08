% manual_region_assignment.m
% Sarah West
% 4/4/22

% Manual assignment of ICs for each mouse based on the regions plotted from
% the LocaNMF atlas. Run after the "FindAtlasRegions" function in the IC
% analysis pipeline.

%% Mouse 1087
mouse = '1087';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal WITH high fine tuning. Second column is the node number
% (common across mice).

region_assignments = [
    
    14, 1
    20, 3
    8, 7
    24, 5
    33, 11
    31, 9
    5, 15
    28, 13
    12, 19
    10, 21
    1, 17
    15, 25
    25, 23
    22, 27
    16, 29
    17, 31 
    % start of right side
    13, 2
    30, 4
    26, 6
    7, 8
    3, 10
    4, 10
    6, 12
    2, 16
    32, 18
    9, 20
    11, 22
    21, 26
    23, 14
    29, 24
    19, 28
    27, 30
    18, 32
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{1} = region_assignments;

node_vector = 1:32; 
IC_vector = 1:33; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(a)
    warning(['Nodes ' num2str(a) ' are missing.']);
end 
if ~isempty(b)
    error(['ICs ' num2str(b) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')


% [From pre- conditional thresholding ICs]

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal. Second column is the index from the
% atlas.

% region_assignments = [ 
% 
% 1, 19;
% 2, 30;
% 3, 16;
% 4, 29;
% 5, 24;
% 6, 18;
% 7, 17;
% 8, 28;
% 9, 33;
% 10, 34;
% 11, 27;
% 12, 12;
% 13, 11;
% 14, 41;
% 15, 49;
% 16, 61;
% 17, 62;
% 18, 46;
% 19, 13;
% 20, 42;
% 21, 45;
% 22, 26;
% 23, 21;
% 24, 35;
% 25, 22;
% 26, 50;
% 27, 25;
% 28, 36;
% 29, 14;
% ];

%% Mouse 1088
mouse = '1088';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal WITH high fine tuning. Second column is the node number
% (common across mice).

region_assignments = [
    % Left side
    13, 1;
    16, 3;
    17, 3;
    19, 5;
    23, 7;
    34, 7;
    6, 11;
    36, 13;
    11, 7;
    26, 9; 
    35, 15;
    1, 17; 
    7, 19;
    9, 21;
    24, 23;
    31, 25;
    28, 27;
    4, 29;
    38, 31; 

    % Right side
    14, 2;
    20, 4;
    29, 6; 
    3, 8;
    33, 10;
    5, 12;
    30, 14;
    2, 16;
    37, 18;
    10, 20;
    8, 22;
    25, 24; 
    12, 26; 
    32, 26;
    27, 28
    21, 30;
    22, 32;
    
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{2} = region_assignments;
node_vector = 1:32; 
IC_vector = 1:38; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')

%% Mouse 1096
mouse = '1096';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\' mouse '\sources.mat'])
figure; imagesc(sources.overlay); colorbar;
title(mouse);

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal WITH high fine tuning. Second column is the node number
% (common across mice).

region_assignments = [
   16	1
38	2
24	3
40	3
23	4
33	5
21	6
15	7
25	8
26	8
27	6
6	9
9	9
7	10
8	10
10	11
12	12
41	13
29	14
5	15
4	16
2	17
3	17
1	18
18	19
19	19
14	20
13	21
39	21
20	22
31	23
37	24
22	25
28	26
36	27
34	28
35	28
30	29
17	30
32	31
11	32   
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{3} = region_assignments;

node_vector = 1:32; 
IC_vector = 1:41; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')

%%  Mouse 1100
mouse = '1100';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)
region_assignments_all{4} = [NaN, NaN];

load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\' mouse '\sources.mat'])
figure; imagesc(sources.overlay); colorbar;
title(mouse);

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal post add back. Second column is the node number
% (common across mice).
region_assignments = [
    % Left side
    17, 1;
    15, 3;
    22, 5;
    28, 7;
    29, 9;
    11, 11;
    18, 13;
    6, 15;
    1, 17;
    31, 19;
    10, 21;
    33, 23;
    32, 25;
    12, 27;
    19, 29;
    5, 31;

    % Right side
    24, 2;
    27, 4;
    25, 6;
    2, 8;
    14, 10;
    7, 12;
    20, 14;
    3, 16;
    30, 18;
    13, 20;
    9, 22;
    23, 24;
    16, 26;
    21, 28;
    26, 30;
    4, 32;
    8, 32;
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{4} = region_assignments;

node_vector = 1:32; 
IC_vector = 1:33; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')

%% Mouse 1106
mouse = '1106';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\' mouse '\sources.mat'])
figure; imagesc(sources.overlay); colorbar;
title(mouse);

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal WITH high fine tuning. Second column is the node number
% (common across mice).

region_assignments = [
    % Left side
    20, 1;
    12, 3;
    22, 5;
    16, 7;
    35, 9;
    7, 11;
    8, 11;
    18, 13;
    1, 15; 
    29, 17;
    13, 19;
    10, 21;
    28, 23;
    15, 25;
    33, 25;
    14, 27;
    31, 29;
    27, 29;
    30, 31;

    % Right side
    32, 2;
    17, 4; 
    26, 6;
    6, 8;
    36, 10;
    3, 12;
    24, 14;
    5, 16; 
    2, 18;
    9, 20;
    11, 22;
    25, 24;
    21, 26; 
    23, 28;
    19, 30;
    4, 32; 
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{5} = region_assignments;
node_vector = 1:32; 
IC_vector = 1:36; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')

%% Mouse 1107
mouse = '1107';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\' mouse '\sources.mat'])
figure; imagesc(sources.overlay); colorbar;
title(mouse);

% Numeric array, because everything's easier that way. First column is IC
% number from after artifact removal WITH high fine tuning. Second column is the node number
% (common across mice).

region_assignments = [
    % Left side
     16, 1;
     30, 3;
     25, 5;
     7, 7;
     8, 7;
     36, 9;
     9, 11;
     17, 13;       % No 13.
     6, 15;
     1, 17;
     2, 17;
     13, 19;
     10, 21;
     34, 23;
     15, 25;
     21, 25;
     32, 27;
     27, 29;
     29, 31;

    % Right side
    31, 2;
    19, 4;
    24, 6;
    4, 8;
    28, 10;
    12, 12;
    33, 14;
    3, 16;
    35, 18;
    14, 20;
    11, 22;
    26, 24;
    20, 26;
    22, 26;
    18, 28;
    5, 30;
    23, 32;  
];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{6} = region_assignments;

node_vector = 1:32; 
IC_vector = 1:36; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')


%% Mouse 1099
mouse = '1099';
dir_out = ['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\'];
mkdir(dir_out)

load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\' mouse '\sources.mat'])
figure; imagesc(sources.overlay); colorbar;
title(mouse);
region_assignments = [
  
    % Left side
    24, 1;
    9, 3;
    31, 5;
    37, 5;
    7, 7;
    44, 9;
    5, 11;
    32, 13;
    42, 13;
    1, 15;
    12, 15;
    4, 17;
    28, 19;
    29, 19;
    19, 21;
    36, 23;
    22, 25;
    30, 27;
    25, 29;
    23, 31;
    33, 31;
    23, 31;

    % Right side
    21, 2;
    8, 4;
    35, 6;
    16, 8;
    17, 8;
    11, 8;
    13, 10;
    14, 10;
    6, 12;
    20, 14;
    3, 16;
    40, 16;
    39, 16;
    2, 18; 
    15, 20;
    18, 22;
    34, 24;
    38, 24;
    27, 26;
    26, 28;
    10, 30;
    41, 32;
    43, 32;
    ];

[~,ordered] = sort(region_assignments(:,2));
region_assignments = region_assignments(ordered,:);
region_assignments_all{8} = region_assignments;

node_vector = 1:32; 
IC_vector = 1:44; 

a = setdiff(IC_vector, region_assignments(:,1)); 
b = setdiff(node_vector, region_assignments(:,2));
if ~isempty(b)
    warning(['Nodes ' num2str(b) ' are missing.']);
end 
if ~isempty(a)
    warning(['ICs ' num2str(a) ' are missing.']);
end
save([dir_out 'region_assignments.mat'], 'region_assignments')


%% count all node occurances
holder = cellfun(@(x) x(:,2), region_assignments_all,'UniformOutput',false);
holder =  cellfun(@unique, holder, 'UniformOutput', false);
holder2 = vertcat(holder{:});

count = NaN(32,1);
for i = 1:32
   
    count(i) = sum(holder2 == i);

end 