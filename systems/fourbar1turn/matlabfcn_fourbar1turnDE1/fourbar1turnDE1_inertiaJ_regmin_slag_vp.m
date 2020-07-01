% Calculate minimal parameter regressor of joint inertia matrix for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% MM_reg [((2+1)*2/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbar1turnDE1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:24
% EndTime: 2020-06-27 16:35:27
% DurationCPUTime: 0.42s
% Computational Cost: add. (3825->47), mult. (5333->122), div. (258->14), fcn. (1498->8), ass. (0->59)
t73 = pkin(4) ^ 2;
t72 = pkin(3) ^ 2;
t71 = -2 * pkin(2);
t43 = pkin(2) ^ 2;
t44 = pkin(1) ^ 2;
t36 = cos(qJ(2));
t66 = pkin(2) * t36;
t58 = -0.2e1 * t66;
t60 = pkin(1) * t58 + t44;
t30 = t43 + t60;
t59 = t72 - t73;
t26 = t30 + t59;
t32 = pkin(1) * t36 - pkin(2);
t35 = sin(qJ(2));
t69 = (-pkin(3) - pkin(4));
t24 = ((pkin(2) - t69) * (pkin(2) + t69)) + t60;
t68 = (-pkin(3) + pkin(4));
t25 = ((pkin(2) - t68) * (pkin(2) + t68)) + t60;
t45 = sqrt(-t24 * t25);
t62 = t35 * t45;
t15 = -pkin(1) * t62 - t32 * t26;
t65 = t35 * pkin(1);
t23 = t26 * t65;
t18 = -t32 * t45 + t23;
t14 = t18 ^ 2;
t28 = 0.1e1 / t30;
t29 = 0.1e1 / t30 ^ 2;
t41 = 0.1e1 / pkin(3);
t46 = t15 ^ 2;
t54 = t28 * t41 * ((t14 + t46) / t72 * t29) ^ (-0.1e1 / 0.2e1);
t4 = (-t15 * t36 + t18 * t35) * t54;
t70 = 0.2e1 * t4;
t67 = pkin(2) * t35;
t57 = pkin(2) * t65;
t64 = 0.1e1 / t45 * (-t24 - t25) * t57;
t63 = t29 / t73;
t61 = t36 * t45;
t56 = t29 * t67;
t38 = 0.1e1 / pkin(4);
t27 = t30 - t59;
t22 = t27 * t67;
t31 = pkin(1) - t66;
t17 = t31 * t45 + t22;
t13 = t17 ^ 2;
t16 = -pkin(2) * t62 + t31 * t27;
t47 = t16 ^ 2;
t8 = (t13 + t47) * t63;
t55 = t28 * t38 * t8 ^ (-0.1e1 / 0.2e1);
t53 = 0.1e1 / t8 * t63;
t12 = 0.1e1 / t47;
t34 = t35 ^ 2;
t1 = 0.2e1 * (-((t31 * t64 + (t36 * t27 + t62) * pkin(2)) * t28 / 0.2e1 + (t43 * t34 * t28 - t17 * t56) * pkin(1)) / t16 - (-(t22 + (-t35 * t64 - t61) * pkin(2)) * t28 / 0.2e1 + (t16 * t29 - t28 * t31) * t57) * t17 * t12) * pkin(4) * t30 * t38 / (t13 * t12 + 0.1e1);
t51 = t1 * t55;
t49 = -0.2e1 * pkin(1) * t55;
t11 = 0.1e1 / t46;
t2 = 0.1e1 + (((0.2e1 * t44 * t34 * pkin(2) - t32 * t64) * t28 + ((t36 * t26 + t62) * t28 - 0.2e1 * t18 * t56) * pkin(1)) / t15 - (t23 * t28 + (-t28 * t61 + ((t32 * t71 - t64) * t28 + t15 * t29 * t71) * t35) * pkin(1)) * t18 * t11) * pkin(3) / (t14 * t11 + 0.1e1) * t30 * t41;
t48 = pkin(2) * t2 * t54;
t3 = (-t15 * t35 - t18 * t36) * t54;
t5 = [1, 0, 0, t34, 0.2e1 * t35 * t36, 0, 0, 0, 0, 0, t3 ^ 2, t3 * t70, 0, t4 ^ 2, 0, 0, t66 * t70, t3 * t58, t13 * t53, -0.2e1 * t17 * t16 * t53, 0, 0, 0, t16 * t49, t17 * t49; 0, 0, 0, 0, 0, t35, t36, 0, 0, 0, 0, 0, t3 * t2, 0, t4 * t2, 0, 0, 0, 0, 0, t17 * t51, -t16 * t51, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2, -0.2e1 * t15 * t48, 0.2e1 * t18 * t48, 0, 0, 0, 0, t1 ^ 2, 0, 0;];
MM_reg = t5;
