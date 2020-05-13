% Calculate position vector for OL to TE/DE
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[12x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = palh3m2IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 2
% StartTime: 2020-05-07 04:55:19
% EndTime: 2020-05-07 04:55:19
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0;];
posNQJ = t1(:);
