% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m1OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energykin_fixb_reg2_slag_vp: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_energykin_fixb_reg2_slag_vp: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energykin_fixb_reg2_slag_vp: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:36:39
% EndTime: 2020-04-15 19:36:40
% DurationCPUTime: 0.33s
% Computational Cost: add. (526->83), mult. (1364->249), div. (0->0), fcn. (1036->20), ass. (0->80)
t72 = qJD(1) ^ 2;
t91 = t72 / 0.2e1;
t50 = qJD(2) + qJD(8);
t90 = pkin(2) * t50;
t51 = qJD(2) + qJD(7);
t89 = pkin(4) * t51;
t88 = cos(qJ(5));
t67 = cos(qJ(6));
t87 = t67 * t72;
t70 = cos(qJ(2));
t86 = t70 * t72;
t71 = qJD(2) ^ 2;
t85 = t71 * pkin(1) ^ 2;
t84 = sin(qJ(10));
t83 = pkin(1) * qJD(2);
t82 = qJD(1) * pkin(15);
t81 = qJD(1) * qJD(2);
t80 = qJD(1) * qJD(6);
t52 = qJD(2) + qJD(3);
t46 = qJD(9) + t50;
t79 = t46 * t90;
t69 = cos(qJ(3));
t78 = t69 * t83;
t58 = sin(qJ(7));
t77 = t58 * t83;
t62 = sin(qJ(3));
t76 = t62 * t83;
t66 = cos(qJ(7));
t75 = t66 * t83;
t63 = sin(qJ(2));
t30 = (t62 * t70 + t63 * t69) * qJD(1);
t33 = (-t62 * t63 + t69 * t70) * qJD(1);
t61 = sin(qJ(4));
t68 = cos(qJ(4));
t14 = t30 * t61 - t68 * t33;
t41 = t63 * qJD(1) * pkin(1) - t82;
t40 = pkin(5) * t52 + t76;
t23 = t61 * t40 - t68 * t78;
t21 = -pkin(5) * t33 + t41;
t22 = t40 * t68 + t61 * t78;
t65 = cos(qJ(8));
t64 = cos(qJ(9));
t60 = sin(qJ(5));
t59 = sin(qJ(6));
t57 = sin(qJ(8));
t56 = sin(qJ(9));
t55 = cos(qJ(10));
t54 = cos(pkin(19));
t53 = sin(pkin(19));
t49 = t50 ^ 2;
t48 = pkin(15) ^ 2 * t91;
t47 = qJD(4) + t52;
t45 = qJD(10) + t51;
t39 = t41 ^ 2 / 0.2e1;
t38 = -t53 * t84 - t54 * t55;
t37 = t53 * t55 - t54 * t84;
t36 = t54 * t89 + t75;
t35 = t53 * t89 + t77;
t32 = (-t58 * t63 + t66 * t70) * qJD(1);
t31 = (-t57 * t63 + t65 * t70) * qJD(1);
t28 = (-t58 * t70 - t63 * t66) * qJD(1);
t26 = (-t57 * t70 - t63 * t65) * qJD(1);
t24 = -pkin(2) * t26 - t82;
t20 = pkin(11) * t47 + t23;
t19 = -pkin(9) * t47 - t22;
t18 = t30 * t68 + t33 * t61;
t16 = t26 * t56 + t31 * t64;
t13 = -t26 * t64 + t31 * t56;
t12 = qJD(5) + t14;
t11 = (-t28 * t54 - t32 * t53) * pkin(4) + t41;
t10 = t18 * t88 + t60 * t47;
t8 = t60 * t18 - t47 * t88;
t7 = t35 * t38 + t36 * t37;
t6 = -t35 * t37 + t36 * t38;
t5 = t28 * t37 + t32 * t38;
t4 = t28 * t38 - t32 * t37;
t3 = pkin(9) * t14 - pkin(11) * t18 + t21;
t2 = t20 * t88 + t60 * t3;
t1 = -t60 * t20 + t3 * t88;
t9 = [0, 0, 0, 0, 0, t91, 0, 0, 0, 0, t70 ^ 2 * t91, -t63 * t86, t70 * t81, t63 ^ 2 * t91, -t63 * t81, t71 / 0.2e1, -t72 * pkin(15) * t63, -pkin(15) * t86, 0, t48, t30 ^ 2 / 0.2e1, t30 * t33, t30 * t52, t33 ^ 2 / 0.2e1, t33 * t52, t52 ^ 2 / 0.2e1, -t33 * t41 + t52 * t76, t30 * t41 + t52 * t78, (-t30 * t62 - t33 * t69) * t83, t39 + (t69 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1) * t85, t18 ^ 2 / 0.2e1, -t18 * t14, t18 * t47, t14 ^ 2 / 0.2e1, -t14 * t47, t47 ^ 2 / 0.2e1, t14 * t21 + t22 * t47, t18 * t21 - t23 * t47, -t14 * t23 - t18 * t22, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t12, t8 ^ 2 / 0.2e1, -t8 * t12, t12 ^ 2 / 0.2e1, t1 * t12 + t19 * t8, t10 * t19 - t12 * t2, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t59 ^ 2 * t91, t59 * t87, t59 * t80, t67 ^ 2 * t91, t67 * t80, qJD(6) ^ 2 / 0.2e1, -pkin(14) * t87, t72 * pkin(14) * t59, 0, pkin(14) ^ 2 * t91, t32 ^ 2 / 0.2e1, t32 * t28, t32 * t51, t28 ^ 2 / 0.2e1, t28 * t51, t51 ^ 2 / 0.2e1, -t28 * t41 + t51 * t75, t32 * t41 - t51 * t77, (t28 * t58 - t32 * t66) * t83, t39 + (t58 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1) * t85, t31 ^ 2 / 0.2e1, t31 * t26, t31 * t50, t26 ^ 2 / 0.2e1, t26 * t50, t49 / 0.2e1, t26 * t82, -t31 * t82, 0, t48, t16 ^ 2 / 0.2e1, -t16 * t13, -t16 * t46, t13 ^ 2 / 0.2e1, t13 * t46, t46 ^ 2 / 0.2e1, -t13 * t24 - t64 * t79, -t16 * t24 + t56 * t79, (-t13 * t56 - t16 * t64) * t90, t24 ^ 2 / 0.2e1 + (t56 ^ 2 / 0.2e1 + t64 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t49, t5 ^ 2 / 0.2e1, t4 * t5, t45 * t5, t4 ^ 2 / 0.2e1, t4 * t45, t45 ^ 2 / 0.2e1, -t11 * t4 + t45 * t6, t11 * t5 - t45 * t7, t4 * t7 - t5 * t6, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t9;
