% Calculate inertial parameters regressor of coriolis joint torque vector for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m2DE_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:45
% EndTime: 2020-05-03 01:06:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (197->54), mult. (535->127), div. (0->0), fcn. (286->6), ass. (0->68)
t21 = (qJD(1) + qJD(4));
t79 = t21 ^ 2;
t55 = (qJD(1) * qJD(2));
t78 = -2 * t55;
t77 = qJD(2) - qJD(3);
t76 = 2 * qJD(1);
t28 = sin(qJ(2));
t75 = pkin(4) * t28;
t27 = sin(qJ(3));
t74 = pkin(5) * t27;
t64 = pkin(4) * qJD(2);
t49 = t28 * t64;
t61 = qJD(3) * t27;
t14 = pkin(5) * t61 + t49;
t26 = sin(qJ(4));
t73 = t26 * t14;
t34 = qJD(1) ^ 2;
t72 = t27 * t34;
t71 = t28 * t34;
t29 = cos(qJ(4));
t70 = t29 * t14;
t30 = cos(qJ(3));
t69 = t30 * t34;
t32 = qJD(3) ^ 2;
t68 = t32 * t30;
t31 = cos(qJ(2));
t33 = qJD(2) ^ 2;
t67 = t33 * t31;
t66 = t27 ^ 2 - t30 ^ 2;
t65 = t28 ^ 2 - t31 ^ 2;
t63 = pkin(1) * t34;
t36 = t27 * t31 - t28 * t30;
t62 = t36 * qJD(3);
t60 = t14 * qJD(1);
t59 = qJD(1) + t21;
t58 = -qJD(4) + t21;
t16 = t31 * pkin(4) + pkin(1);
t4 = t77 * t36;
t57 = pkin(4) * pkin(5) * t4;
t53 = pkin(4) * t67;
t13 = pkin(5) * t68 + t53;
t56 = qJD(4) * t70 + t26 * t13 - t29 * t60;
t54 = qJD(1) * qJD(3);
t18 = pkin(4) * t71;
t52 = t27 * t69;
t51 = t31 * t71;
t50 = -0.2e1 * t60;
t15 = pkin(2) + t16;
t12 = t30 * pkin(5) + t15;
t9 = pkin(3) + t12;
t48 = t26 * t9 * qJD(1);
t47 = t28 * t55;
t46 = qJD(2) * t62;
t45 = t59 * t9;
t44 = pkin(5) * t46;
t43 = pkin(4) * t46;
t42 = t27 * t30 * t54;
t41 = t31 * t47;
t40 = pkin(1) * t78;
t39 = -0.2e1 * pkin(4) * t47;
t38 = t79 * t26;
t37 = t79 * t29;
t35 = t27 * t28 + t30 * t31;
t8 = t29 * t13;
t3 = t77 * t35;
t2 = qJD(2) * t57;
t1 = qJD(3) * t57;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t41, t65 * t78, t67, -0.2e1 * t41, -t33 * t28, 0, t28 * t40, t31 * t40, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t53, t16 * t39, 0.2e1 * t42, -0.2e1 * t66 * t54, t68, -0.2e1 * t42, -t32 * t27, 0, (-t15 * t61 - t30 * t49) * t76, (-qJD(3) * t15 * t30 + t27 * t49) * t76, -t53, t15 * t39, 0, 0, 0, 0, 0, 0, t50, 0, -t13, t12 * t50, 0, 0, 0, 0, 0, 0, -qJD(4) * t26 * t45 - t21 * t70 + t56, t8 + t59 * t73 + (-t29 * t45 - t73) * qJD(4), 0, t9 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t65 * t34, 0, t51, 0, 0, t28 * t63, t31 * t63, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, t16 * t18, 0, 0, 0, 0, 0, 0, (t28 * t69 + (-qJD(2) * t36 + t4) * qJD(3)) * pkin(4), (-t27 * t71 + (-qJD(2) * t35 + t3) * qJD(3)) * pkin(4), 0, t15 * t18, 0, 0, 0, 0, 0, 0, t18, 0, 0, t1 + (t12 * t71 - t44) * pkin(4), 0, 0, 0, 0, 0, 0, t37 * t75, -t38 * t75, 0, t1 + (t71 * t9 - t44) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t66 * t34, 0, t52, 0, 0, t15 * t72 + (t4 + t62) * t64, t15 * t69 + (qJD(3) * t35 + t3) * t64, 0, 0, 0, 0, 0, 0, 0, 0, pkin(5) * t72, 0, 0, t2 + (t12 * t72 + t43) * pkin(5), 0, 0, 0, 0, 0, 0, t37 * t74, -t38 * t74, 0, t2 + (t72 * t9 + t43) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t48 - t21 * (-t48 + t70) + t56, t8 + t58 * t73 + (t29 * t58 * t9 + t73) * qJD(1), 0, 0;];
tauc_reg = t5;
