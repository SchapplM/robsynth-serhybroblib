% Calculate position vector for OL to TE/DE
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[16x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = palh1m2IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-05-02 23:42:14
% EndTime: 2020-05-02 23:42:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0;];
posNQJ = t1(:);
