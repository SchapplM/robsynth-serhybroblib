% Calculate position vector for OL to TE/DE
% palh3m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[12x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = palh3m1IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 1; 1;];
posNQJ = t1(:);
