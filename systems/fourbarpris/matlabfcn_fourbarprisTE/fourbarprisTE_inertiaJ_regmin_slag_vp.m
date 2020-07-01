% Calculate minimal parameter regressor of joint inertia matrix for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% MM_reg [((1+1)*1/2)x9]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbarprisTE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_inertiaJ_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_inertiaJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:03
% EndTime: 2020-06-27 17:07:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (84->14), mult. (42->15), div. (20->7), fcn. (2->2), ass. (0->9)
t13 = (qJ(1) + pkin(3));
t18 = -pkin(2) + t13;
t19 = -pkin(2) - t13;
t17 = (0.1e1 / (pkin(1) - t19) / (pkin(1) + t18) / (pkin(1) + t19) / (pkin(1) - t18));
t15 = (t13 ^ 2);
t14 = (qJ(1) ^ 2);
t2 = pkin(1) ^ 2 - pkin(2) ^ 2 - t14 - (2 * qJ(1) + pkin(3)) * pkin(3);
t16 = t17 * t2 ^ 2 / t15;
t1 = [-t16, 0, 0, 2 * t2 / t13 * (-1 / t17) ^ (-0.1e1 / 0.2e1), -2 * qJ(1) * t16, -t14 * t16 + 1, -4 * t15 * t17, 0, 0;];
MM_reg = t1;
