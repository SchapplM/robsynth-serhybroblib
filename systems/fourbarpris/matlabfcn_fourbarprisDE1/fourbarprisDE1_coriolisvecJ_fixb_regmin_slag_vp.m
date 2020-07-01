% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [1x9]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbarprisDE1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:54
% EndTime: 2020-06-27 17:14:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (652->43), mult. (400->108), div. (140->12), fcn. (10->2), ass. (0->51)
t30 = (qJD(1) ^ 2);
t31 = (qJ(1) ^ 2);
t5 = pkin(1) ^ 2 - pkin(2) ^ 2 - t31 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t66 = t5 ^ 2;
t60 = t66 * t30;
t29 = qJ(1) + pkin(3);
t49 = -pkin(2) + t29;
t19 = pkin(1) - t49;
t65 = 1 / t19;
t50 = -pkin(2) - t29;
t21 = pkin(1) + t50;
t64 = 1 / t21;
t26 = 1 / t29;
t18 = pkin(1) - t50;
t20 = pkin(1) + t49;
t59 = t20 * t21;
t48 = t19 * t59;
t45 = t18 * t48;
t1 = (-t45) ^ (-0.1e1 / 0.2e1);
t11 = 1 / t19 ^ 2;
t15 = 1 / t21 ^ 2;
t8 = 1 / t18;
t63 = t5 * t8;
t25 = t29 ^ 2;
t62 = t25 * t8;
t27 = 1 / t29 ^ 2;
t61 = t27 * t8;
t58 = t29 * t30;
t7 = t29 * qJD(1);
t57 = qJD(1) * t7;
t56 = 2 * t63;
t55 = t8 * t60;
t9 = 1 / t18 ^ 2;
t54 = t9 * t60;
t53 = t64 * t62;
t52 = t31 * t61;
t51 = 4 * t57;
t47 = t27 * t55;
t28 = t26 / t25;
t46 = t28 * t55;
t44 = t64 * t47;
t43 = t52 * t60;
t42 = t15 * t47;
t12 = 1 / t20;
t41 = t12 * t44;
t13 = 1 / t20 ^ 2;
t40 = t13 * t44;
t39 = t43 / 0.2e1;
t38 = t11 * t41;
t37 = t1 / t45 * t5 * (-t48 + (t59 + (t20 - t21) * t19) * t18);
t2 = [-t38 / 0.2e1 + (t40 / 0.2e1 + (-t42 / 0.2e1 + (t46 + (t54 / 0.2e1 + ((2 * t57 - t58) * t56)) * t27) * t64) * t12) * t65, 0, 0, (-(t26 * t37) / 0.2e1 + ((2 * t26 * t29 - t5 * t27) * t1)) * t30 + ((qJD(1) * t37 - 4 * t7 * t1) * t26 * qJD(1)), -t65 * t41 + (-t38 + (t40 + (-t42 + (2 * t46 + (t54 + (t51 - 2 * t58) * t56) * t27) * t64) * t12) * t65) * qJ(1), t64 * t13 * t65 * t39 + (-(t15 * t65 * t43) / 0.2e1 - (t11 * t39 - ((t5 * t51 * t52) + (-(t66 * qJ(1) * t61) + (-(2 * t26 * t63) + ((t28 * t8) + (t27 * t9) / 0.2e1) * t66) * t31) * t30) * t65) * t64) * t12, 2 * (-t11 * t12 * t53 + (t13 * t53 + (-t15 * t62 + (t25 * t9 - 2 * t29 * t8) * t64) * t12) * t65) * t30, 0, 0;];
tauc_reg = t2;
