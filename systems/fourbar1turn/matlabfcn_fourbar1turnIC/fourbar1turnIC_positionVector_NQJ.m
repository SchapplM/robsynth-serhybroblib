% Calculate position vector for OL to TE/DE
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[6x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = fourbar1turnIC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:12
% EndTime: 2020-05-07 11:33:12
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 0;];
posNQJ = t1(:);
