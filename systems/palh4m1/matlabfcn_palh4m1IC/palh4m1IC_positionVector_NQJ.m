% Calculate position vector for OL to TE/DE
% palh4m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% 
% Output:
% posNQJ[9x1]
%   position vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function posNQJ = palh4m1IC_positionVector_NQJ()
%% Coder Information
%#codegen
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->9)
unknown=NaN(9,1);
unknown(1,1) = 1;
unknown(2,1) = 1;
unknown(3,1) = 1;
unknown(4,1) = 1;
unknown(5,1) = 1;
unknown(6,1) = 1;
unknown(7,1) = 1;
unknown(8,1) = 1;
unknown(9,1) = 0;
posNQJ  = unknown(:);
