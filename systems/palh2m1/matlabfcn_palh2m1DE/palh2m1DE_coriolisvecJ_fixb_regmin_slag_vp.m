% Calculate minimal parameter regressor of coriolis joint torque vector for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m1DE_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:35
% EndTime: 2020-06-30 17:39:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (312->66), mult. (768->140), div. (0->0), fcn. (552->6), ass. (0->69)
t31 = qJD(2) + qJD(3);
t39 = sin(qJ(2));
t41 = cos(qJ(3));
t38 = sin(qJ(3));
t42 = cos(qJ(2));
t73 = t42 * t38;
t48 = t41 * t39 + t73;
t14 = t31 * t48;
t80 = 2 * qJD(1);
t32 = qJD(1) + qJD(4);
t40 = cos(qJ(4));
t84 = t32 * t40;
t72 = t42 * t41;
t75 = t39 * t38;
t23 = -0.2e1 * t72 * t75;
t35 = t41 ^ 2;
t36 = t42 ^ 2;
t83 = t23 + (-t38 ^ 2 + t35) * (t36 - 0.1e1 / 0.2e1);
t70 = t39 ^ 2 - t36;
t82 = t23 - t70 * (t35 - 0.1e1 / 0.2e1);
t44 = qJD(1) ^ 2;
t81 = -2 * t44;
t22 = t72 - t75;
t69 = pkin(3) * qJD(3);
t19 = t22 * t69;
t28 = t41 * pkin(3) + pkin(2);
t50 = -pkin(3) * t75 + t28 * t42;
t46 = t22 * qJD(3);
t65 = qJD(2) * t42;
t66 = qJD(2) * t39;
t9 = t28 * t65 + (-t38 * t66 + t46) * pkin(3);
t79 = -t50 * qJD(2) - t19 + t9;
t17 = pkin(1) + pkin(4) + t50;
t78 = t17 * t40;
t29 = t42 * pkin(2) + pkin(1);
t77 = t29 * t44;
t18 = pkin(3) * t73 + t39 * t28;
t11 = t18 * qJD(2) + t48 * t69;
t37 = sin(qJ(4));
t76 = t37 * t11;
t74 = t39 * t44;
t71 = t42 * t44;
t68 = qJD(1) * t37;
t67 = qJD(1) * t40;
t64 = qJD(4) * t40;
t63 = -qJD(1) - t32;
t62 = -qJD(4) + t32;
t47 = t22 * qJD(2);
t13 = t46 + t47;
t3 = t9 * qJD(2) + t13 * t69;
t8 = -t28 * t66 + (-t48 * qJD(3) - t38 * t65) * pkin(3);
t60 = t11 * t64 + t37 * t3 + t8 * t67;
t59 = qJD(1) * qJD(2);
t58 = t44 * t48 * t22;
t57 = pkin(2) * t66;
t56 = t17 * t68;
t55 = 0.2e1 * t59;
t53 = t63 * t37;
t52 = -0.2e1 * pkin(1) * t59;
t51 = qJD(3) * (-qJD(2) - t31);
t49 = pkin(2) * qJD(2) * (-qJD(3) + t31);
t43 = qJD(2) ^ 2;
t20 = t48 * pkin(3);
t16 = t22 * t77;
t15 = t48 * t77;
t12 = pkin(3) * t47 + t19;
t4 = (t22 * t31 - t13) * qJD(1);
t2 = t40 * t3;
t1 = [0, 0, 0, t39 * t42 * t55, -t70 * t55, -t43 * t42, t43 * t39, 0, t39 * t52, t42 * t52, t48 * t13 * t80, 0.4e1 * (t82 * qJD(2) + t83 * qJD(3)) * qJD(1), -t31 * t13, t31 * t14, 0, (-t14 * t29 - t22 * t57) * t80, (-t13 * t29 + t48 * t57) * t80, t8 * t80, 0, t17 * qJD(4) * t53 + t8 * t84 + t60, t2 + t8 * t53 + (t63 * t78 - t76) * qJD(4); 0, 0, 0, -t39 * t71, t70 * t44, 0, 0, 0, pkin(1) * t74, pkin(1) * t71, -t58, t82 * t81, t4, 0, 0, t15 + (t22 * t74 + t38 * t51) * pkin(2), t16 + (t41 * t51 - t48 * t74) * pkin(2), t18 * t44, 0, (t18 * t84 + t79 * t37) * t32, (-t32 * t37 * t18 + t79 * t40) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t83 * t81, t4, 0, 0, t38 * t49 + t15, t41 * t49 + t16, t20 * t44, 0, (t20 * t67 - t37 * t12 + (t13 * t37 + t48 * t64) * pkin(3)) * t32, (-t20 * t68 - t40 * t12 + (-qJD(4) * t37 * t48 + t13 * t40) * pkin(3)) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t56 - (t40 * t11 - t56) * t32 + t60, t2 + t62 * t76 + (-t37 * t8 + t62 * t78) * qJD(1);];
tauc_reg = t1;
