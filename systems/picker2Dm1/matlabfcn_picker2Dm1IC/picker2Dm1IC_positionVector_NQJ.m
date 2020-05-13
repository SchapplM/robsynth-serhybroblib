% Calculate position vector for OL to TE/DE
% picker2Dm1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[15x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = picker2Dm1IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-05-11 05:54:43
% EndTime: 2020-05-11 05:54:43
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0;];
posNQJ = t1(:);
