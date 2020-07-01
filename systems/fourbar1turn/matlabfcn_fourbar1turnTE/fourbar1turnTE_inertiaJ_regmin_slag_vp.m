% Calculate minimal parameter regressor of joint inertia matrix for
% fourbar1turnTE
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
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = fourbar1turnTE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:33
% EndTime: 2020-06-27 16:22:37
% DurationCPUTime: 0.38s
% Computational Cost: add. (2796->47), mult. (3816->120), div. (160->12), fcn. (1070->4), ass. (0->54)
t64 = pkin(4) ^ 2;
t63 = -2 * pkin(2);
t37 = pkin(2) ^ 2;
t38 = pkin(1) ^ 2;
t31 = cos(qJ(2));
t57 = pkin(2) * t31;
t47 = -0.2e1 * t57;
t49 = pkin(1) * t47 + t38;
t25 = t37 + t49;
t48 = pkin(3) ^ 2 - t64;
t21 = t25 + t48;
t27 = pkin(1) * t31 - pkin(2);
t30 = sin(qJ(2));
t60 = (-pkin(3) - pkin(4));
t19 = ((pkin(2) - t60) * (pkin(2) + t60)) + t49;
t59 = (-pkin(3) + pkin(4));
t20 = ((pkin(2) - t59) * (pkin(2) + t59)) + t49;
t39 = sqrt(-t19 * t20);
t51 = t30 * t39;
t10 = -pkin(1) * t51 - t27 * t21;
t56 = t30 * pkin(1);
t18 = t21 * t56;
t13 = -t27 * t39 + t18;
t23 = 0.1e1 / t25;
t36 = 0.1e1 / pkin(3);
t53 = t23 * t36;
t61 = -t31 / 0.2e1;
t4 = (t30 * t13 / 0.2e1 + t10 * t61) * t53;
t62 = 0.2e1 * t4;
t58 = pkin(2) * t30;
t46 = pkin(2) * t56;
t55 = 0.1e1 / t39 * (-t19 - t20) * t46;
t33 = 0.1e1 / pkin(4);
t54 = t23 * t33;
t24 = 0.1e1 / t25 ^ 2;
t52 = t24 / t64;
t50 = t31 * t39;
t45 = t24 * t58;
t22 = t25 - t48;
t26 = pkin(1) - t57;
t11 = -pkin(2) * t51 + t26 * t22;
t44 = t11 * t54;
t17 = t22 * t58;
t12 = t26 * t39 + t17;
t43 = t12 * t54;
t29 = t30 ^ 2;
t7 = 0.1e1 / t10 ^ 2;
t2 = 0.1e1 + (((0.2e1 * t38 * t29 * pkin(2) - t27 * t55) * t23 + ((t31 * t21 + t51) * t23 - 0.2e1 * t13 * t45) * pkin(1)) / t10 - (t18 * t23 + (-t23 * t50 + ((t27 * t63 - t55) * t23 + t10 * t24 * t63) * t30) * pkin(1)) * t13 * t7) * pkin(3) * t25 * t36 / (t13 ^ 2 * t7 + 0.1e1);
t42 = pkin(2) * t2 * t53;
t9 = t12 ^ 2;
t8 = 0.1e1 / t11 ^ 2;
t3 = (-t30 * t10 / 0.2e1 + t13 * t61) * t53;
t1 = 0.2e1 * (-((t26 * t55 + (t31 * t22 + t51) * pkin(2)) * t23 / 0.2e1 + (t37 * t29 * t23 - t12 * t45) * pkin(1)) / t11 - (-(t17 + (-t30 * t55 - t50) * pkin(2)) * t23 / 0.2e1 + (t11 * t24 - t23 * t26) * t46) * t12 * t8) * pkin(4) * t25 * t33 / (t9 * t8 + 0.1e1);
t5 = [1, 0, 0, t29, 0.2e1 * t30 * t31, 0, 0, 0, 0, 0, t3 ^ 2, t3 * t62, 0, t4 ^ 2, 0, 0, t57 * t62, t3 * t47, t9 * t52 / 0.4e1, -t12 * t11 * t52 / 0.2e1, 0, 0, 0, -pkin(1) * t44, -pkin(1) * t43; 0, 0, 0, 0, 0, t30, t31, 0, 0, 0, 0, 0, t3 * t2, 0, t4 * t2, 0, 0, 0, 0, 0, t1 * t43 / 0.2e1, -t1 * t44 / 0.2e1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2, -t10 * t42, t13 * t42, 0, 0, 0, 0, t1 ^ 2, 0, 0;];
MM_reg = t5;
