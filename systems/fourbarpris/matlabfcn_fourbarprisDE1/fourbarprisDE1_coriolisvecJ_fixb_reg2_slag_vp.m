% Calculate inertial parameters regressor of coriolis joint torque vector for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% tauc_reg [1x(1*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbarprisDE1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:10:10
% EndTime: 2020-05-07 09:10:13
% DurationCPUTime: 0.24s
% Computational Cost: add. (779->43), mult. (485->108), div. (175->12), fcn. (10->2), ass. (0->52)
t31 = (qJD(1) ^ 2);
t32 = (qJ(1) ^ 2);
t6 = pkin(1) ^ 2 - pkin(2) ^ 2 - t32 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t67 = t6 ^ 2;
t61 = t31 * t67;
t30 = qJ(1) + pkin(3);
t50 = -pkin(2) + t30;
t20 = pkin(1) - t50;
t66 = 1 / t20;
t51 = -pkin(2) - t30;
t22 = pkin(1) + t51;
t65 = 1 / t22;
t27 = 1 / t30;
t19 = pkin(1) - t51;
t21 = pkin(1) + t50;
t60 = t21 * t22;
t49 = t20 * t60;
t46 = t19 * t49;
t2 = (-t46) ^ (-0.1e1 / 0.2e1);
t12 = 1 / t20 ^ 2;
t16 = 1 / t22 ^ 2;
t9 = 1 / t19;
t64 = t6 * t9;
t26 = t30 ^ 2;
t63 = t26 * t9;
t28 = 1 / t30 ^ 2;
t62 = t28 * t9;
t59 = t30 * t31;
t8 = t30 * qJD(1);
t58 = qJD(1) * t8;
t57 = 2 * t64;
t56 = t9 * t61;
t55 = t65 * t63;
t54 = t32 * t62;
t10 = 1 / t19 ^ 2;
t53 = t10 * t61;
t52 = 4 * t58;
t48 = t28 * t56;
t29 = t27 / t26;
t47 = t29 * t56;
t45 = t65 * t48;
t44 = t54 * t61;
t43 = t16 * t48;
t13 = 1 / t21;
t42 = t13 * t45;
t14 = 1 / t21 ^ 2;
t41 = t14 * t45;
t40 = t44 / 0.2e1;
t39 = t12 * t42;
t38 = t2 / t46 * t6 * (-t49 + (t60 + (t21 - t22) * t20) * t19);
t1 = -t39 / 0.2e1 + (t41 / 0.2e1 + (-t43 / 0.2e1 + (t47 + (t53 / 0.2e1 + ((2 * t58 - t59) * t57)) * t28) * t65) * t13) * t66;
t3 = [0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, (-(t27 * t38) / 0.2e1 + ((2 * t27 * t30 - t6 * t28) * t2)) * t31 + ((qJD(1) * t38 - 4 * t8 * t2) * t27 * qJD(1)), 0, -t66 * t42 + (-t39 + (t41 + (-t43 + (2 * t47 + (t53 + (t52 - 2 * t59) * t57) * t28) * t65) * t13) * t66) * qJ(1), t65 * t14 * t66 * t40 + (-(t16 * t66 * t44) / 0.2e1 - (t12 * t40 - ((t6 * t52 * t54) + (-(t67 * qJ(1) * t62) + (-(2 * t27 * t64) + ((t29 * t9) + (t28 * t10) / 0.2e1) * t67) * t32) * t31) * t66) * t65) * t13, 0, 0, 0, 0, 0, 2 * (-t12 * t13 * t55 + (t14 * t55 + (-t16 * t63 + (t10 * t26 - 2 * t30 * t9) * t65) * t13) * t66) * t31, 0, 0, 0, 0;];
tauc_reg = t3;
