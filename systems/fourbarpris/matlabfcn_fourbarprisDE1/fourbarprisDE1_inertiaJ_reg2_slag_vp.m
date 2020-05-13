% Calculate inertial parameters regressor of joint inertia matrix for
% fourbarprisDE1
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
% MM_reg [((1+1)*1/2)x(1*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbarprisDE1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaJ_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaJ_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t13 = (qJ(1) + pkin(3));
t18 = -pkin(2) + t13;
t19 = -pkin(2) - t13;
t17 = (0.1e1 / (pkin(1) - t19) / (pkin(1) + t18) / (pkin(1) + t19) / (pkin(1) - t18));
t15 = (t13 ^ 2);
t14 = (qJ(1) ^ 2);
t2 = pkin(1) ^ 2 - pkin(2) ^ 2 - t14 - (2 * qJ(1) + pkin(3)) * pkin(3);
t16 = t17 * t2 ^ 2 / t15;
t1 = [0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 2 * t2 / t13 * (-1 / t17) ^ (-0.1e1 / 0.2e1), 0, -2 * qJ(1) * t16, -t14 * t16 + 1, 0, 0, 0, 0, 0, -4 * t15 * t17, 0, 0, 0, 0;];
MM_reg = t1;
