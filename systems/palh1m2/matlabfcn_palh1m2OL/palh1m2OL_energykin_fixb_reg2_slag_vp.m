% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_reg [1x(13*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m2OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energykin_fixb_reg2_slag_vp: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_energykin_fixb_reg2_slag_vp: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energykin_fixb_reg2_slag_vp: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:46:52
% EndTime: 2020-05-02 21:46:52
% DurationCPUTime: 0.37s
% Computational Cost: add. (526->87), mult. (1331->257), div. (0->0), fcn. (1042->20), ass. (0->84)
t73 = qJD(1) ^ 2;
t96 = t73 / 0.2e1;
t70 = cos(qJ(3));
t95 = pkin(1) * t70;
t51 = qJD(2) + qJD(8);
t94 = pkin(2) * t51;
t52 = qJD(2) + qJD(7);
t93 = pkin(4) * t52;
t92 = cos(qJ(5));
t68 = cos(qJ(6));
t91 = t68 * t73;
t71 = cos(qJ(2));
t90 = t71 * t70;
t89 = t71 * t73;
t72 = qJD(2) ^ 2;
t88 = t72 * pkin(1) ^ 2;
t63 = sin(qJ(3));
t44 = t63 * pkin(1) + pkin(5);
t62 = sin(qJ(4));
t69 = cos(qJ(4));
t85 = pkin(5) * qJD(3);
t21 = (t44 * t62 - t69 * t95) * qJD(2) + t62 * t85;
t87 = sin(qJ(10));
t86 = pkin(1) * qJD(2);
t84 = qJD(1) * pkin(15);
t64 = sin(qJ(2));
t26 = -pkin(5) * t90 + (t63 * pkin(5) + pkin(1)) * t64 - pkin(15);
t83 = t26 * qJD(1);
t43 = t64 * pkin(1) - pkin(15);
t82 = t43 * qJD(1);
t81 = qJD(1) * qJD(2);
t80 = qJD(1) * qJD(6);
t53 = qJD(2) + qJD(3);
t47 = qJD(9) + t51;
t79 = t47 * t94;
t78 = t53 * t86;
t59 = sin(qJ(7));
t77 = t59 * t86;
t67 = cos(qJ(7));
t76 = t67 * t86;
t32 = (t63 * t71 + t64 * t70) * qJD(1);
t35 = (-t63 * t64 + t90) * qJD(1);
t14 = t62 * t32 - t69 * t35;
t22 = (t44 * t69 + t62 * t95) * qJD(2) + t69 * t85;
t66 = cos(qJ(8));
t65 = cos(qJ(9));
t61 = sin(qJ(5));
t60 = sin(qJ(6));
t58 = sin(qJ(8));
t57 = sin(qJ(9));
t56 = cos(qJ(10));
t55 = cos(pkin(19));
t54 = sin(pkin(19));
t50 = t51 ^ 2;
t49 = pkin(15) ^ 2 * t96;
t48 = qJD(4) + t53;
t46 = qJD(10) + t52;
t40 = t43 ^ 2 * t96;
t39 = -t54 * t87 - t55 * t56;
t38 = t54 * t56 - t55 * t87;
t37 = t55 * t93 + t76;
t36 = t54 * t93 + t77;
t34 = (-t59 * t64 + t67 * t71) * qJD(1);
t33 = (-t58 * t64 + t66 * t71) * qJD(1);
t30 = (-t59 * t71 - t64 * t67) * qJD(1);
t28 = (-t58 * t71 - t64 * t66) * qJD(1);
t23 = -t28 * pkin(2) - t84;
t20 = -t48 * pkin(9) - t22;
t19 = t48 * pkin(11) + t21;
t18 = t69 * t32 + t62 * t35;
t16 = t57 * t28 + t65 * t33;
t13 = -t65 * t28 + t57 * t33;
t12 = qJD(5) + t14;
t11 = t82 + (-t30 * t55 - t34 * t54) * pkin(4);
t10 = t92 * t18 + t61 * t48;
t8 = t61 * t18 - t92 * t48;
t7 = t39 * t36 + t38 * t37;
t6 = -t38 * t36 + t39 * t37;
t5 = t38 * t30 + t39 * t34;
t4 = t39 * t30 - t38 * t34;
t3 = t14 * pkin(9) - t18 * pkin(11) + t83;
t2 = t92 * t19 + t61 * t3;
t1 = -t61 * t19 + t92 * t3;
t9 = [0, 0, 0, 0, 0, t96, 0, 0, 0, 0, t71 ^ 2 * t96, -t64 * t89, t71 * t81, t64 ^ 2 * t96, -t64 * t81, t72 / 0.2e1, -t73 * pkin(15) * t64, -pkin(15) * t89, 0, t49, t32 ^ 2 / 0.2e1, t32 * t35, t32 * t53, t35 ^ 2 / 0.2e1, t35 * t53, t53 ^ 2 / 0.2e1, -t35 * t82 + t63 * t78, t32 * t82 + t70 * t78, (-t32 * t63 - t35 * t70) * t86, t40 + (t70 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1) * t88, t18 ^ 2 / 0.2e1, -t18 * t14, t18 * t48, t14 ^ 2 / 0.2e1, -t14 * t48, t48 ^ 2 / 0.2e1, t14 * t83 + t22 * t48, t18 * t83 - t21 * t48, -t21 * t14 - t22 * t18, t21 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t26 ^ 2 * t96, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t12, t8 ^ 2 / 0.2e1, -t8 * t12, t12 ^ 2 / 0.2e1, t1 * t12 + t20 * t8, t20 * t10 - t2 * t12, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t60 ^ 2 * t96, t60 * t91, t60 * t80, t68 ^ 2 * t96, t68 * t80, qJD(6) ^ 2 / 0.2e1, -pkin(14) * t91, pkin(14) * t73 * t60, 0, pkin(14) ^ 2 * t96, t34 ^ 2 / 0.2e1, t34 * t30, t34 * t52, t30 ^ 2 / 0.2e1, t30 * t52, t52 ^ 2 / 0.2e1, -t30 * t82 + t52 * t76, t34 * t82 - t52 * t77, (t30 * t59 - t34 * t67) * t86, t40 + (t59 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1) * t88, t33 ^ 2 / 0.2e1, t33 * t28, t33 * t51, t28 ^ 2 / 0.2e1, t28 * t51, t50 / 0.2e1, t28 * t84, -t33 * t84, 0, t49, t16 ^ 2 / 0.2e1, -t16 * t13, -t16 * t47, t13 ^ 2 / 0.2e1, t13 * t47, t47 ^ 2 / 0.2e1, -t23 * t13 - t65 * t79, -t23 * t16 + t57 * t79, (-t13 * t57 - t16 * t65) * t94, t23 ^ 2 / 0.2e1 + (t57 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t50, t5 ^ 2 / 0.2e1, t5 * t4, t46 * t5, t4 ^ 2 / 0.2e1, t46 * t4, t46 ^ 2 / 0.2e1, -t11 * t4 + t6 * t46, t11 * t5 - t7 * t46, t7 * t4 - t6 * t5, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t9;
