% Implicit kinematic constraints of
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = palh2m2IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:48:58
% EndTime: 2020-05-03 06:48:58
% DurationCPUTime: 0.03s
% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [qJ(2) + qJ(3); qJ(4) + qJ(5);];
h = t1;
