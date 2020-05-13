% Calculate position vector for OL to TE/DE
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[6x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = fivebar1IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 0;];
posNQJ = t1(:);
