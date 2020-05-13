% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh1m2TE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energykin_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:47:04
% EndTime: 2020-05-01 20:47:05
% DurationCPUTime: 0.46s
% Computational Cost: add. (892->78), mult. (1828->223), div. (0->0), fcn. (2061->22), ass. (0->96)
t65 = cos(pkin(20));
t70 = sin(pkin(18));
t75 = cos(pkin(18));
t96 = sin(pkin(20));
t33 = t70 * t65 - t75 * t96;
t36 = t75 * t65 + t70 * t96;
t59 = pkin(22) + pkin(21);
t53 = sin(t59);
t54 = cos(t59);
t15 = (-t33 * t54 + t36 * t53) * qJD(1);
t13 = t15 ^ 2;
t109 = t13 / 0.2e1;
t78 = qJD(1) ^ 2;
t108 = t78 / 0.2e1;
t67 = sin(qJ(4));
t107 = t15 * t67;
t72 = cos(qJ(4));
t106 = t15 * t72;
t71 = sin(pkin(17));
t76 = cos(pkin(17));
t37 = t70 * t76 - t75 * t71;
t38 = t71 * t70 + t76 * t75;
t69 = sin(qJ(2));
t74 = cos(qJ(2));
t26 = t69 * t37 + t38 * t74;
t105 = t26 * t78;
t51 = pkin(1) * t69 - pkin(15);
t104 = t51 * t78;
t62 = qJ(3) + qJ(2);
t56 = cos(t62);
t103 = t56 * t78;
t102 = t74 * t78;
t77 = qJD(2) ^ 2;
t101 = t77 * pkin(1) ^ 2;
t68 = sin(qJ(3));
t73 = cos(qJ(3));
t100 = (t68 ^ 2 + t73 ^ 2) * t101 / 0.2e1;
t99 = pkin(1) * qJD(2);
t19 = t33 * t73 - t68 * t36;
t22 = t68 * t33 + t36 * t73;
t8 = t19 * t74 - t22 * t69;
t9 = t19 * t69 + t22 * t74;
t80 = -t53 * t8 - t9 * t54;
t60 = qJD(2) + qJD(3);
t88 = t68 * t99;
t81 = t60 * pkin(5) + t88;
t82 = -t53 * t9 + t8 * t54;
t89 = t73 * t99;
t3 = -t80 * t81 + t82 * t89;
t98 = t3 ^ 2 / 0.2e1;
t97 = cos(pkin(22));
t95 = sin(pkin(22));
t94 = qJD(1) * pkin(15);
t64 = sin(pkin(19));
t66 = cos(pkin(19));
t31 = t73 * t64 + t68 * t66;
t34 = -t68 * t64 + t73 * t66;
t18 = (-t31 * t69 + t34 * t74) * qJD(1);
t14 = -t18 * pkin(2) - t94;
t93 = qJD(1) * t14;
t92 = qJD(1) * t60;
t41 = t51 * qJD(1);
t91 = qJD(1) * qJD(2);
t90 = t69 * t102;
t86 = pkin(1) * t91;
t30 = (-t68 * t69 + t73 * t74) * qJD(1);
t27 = -t30 * pkin(5) + t41;
t85 = t69 * t91;
t84 = pkin(2) * t60;
t83 = t60 * t89;
t16 = (-t33 * t53 - t36 * t54) * qJD(1);
t35 = -t70 * t97 + t75 * t95;
t61 = t77 / 0.2e1;
t57 = pkin(15) ^ 2 * t108;
t55 = sin(t62);
t52 = t60 ^ 2 / 0.2e1;
t50 = t74 * t91;
t49 = t74 ^ 2 * t108;
t48 = t69 ^ 2 * t108;
t40 = t60 * t88;
t39 = t51 ^ 2 * t108;
t32 = -t70 * t95 - t97 * t75;
t29 = (t68 * t74 + t69 * t73) * qJD(1);
t25 = -t37 * t74 + t69 * t38;
t24 = t34 * t84;
t23 = t31 * t84;
t21 = t69 * t32 - t35 * t74;
t20 = -t32 * t74 - t35 * t69;
t17 = (t31 * t74 + t34 * t69) * qJD(1);
t11 = -qJD(4) + t16;
t10 = (-cos(pkin(21)) * t32 - t35 * sin(pkin(21))) * pkin(4) + t51;
t7 = -t16 * pkin(9) - t15 * pkin(11) + t27;
t5 = -t80 * t89 - t82 * t81;
t2 = t72 * t5 + t67 * t7;
t1 = -t67 * t5 + t72 * t7;
t4 = [0, 0, 0, 0, 0, t108, 0, 0, 0, 0, t49, -t90, t50, t48, -t85, t61, -t78 * pkin(15) * t69, -pkin(15) * t102, 0, t57, t29 ^ 2 / 0.2e1, t29 * t30, t60 * t29, t30 ^ 2 / 0.2e1, t60 * t30, t52, -t30 * t41 + t40, t29 * t41 + t83, (-t29 * t68 - t30 * t73) * t99, t39 + t100, t109, t16 * t15, 0, t16 ^ 2 / 0.2e1, 0, 0, -t27 * t16, t27 * t15, t3 * t15 + t5 * t16, t5 ^ 2 / 0.2e1 + t98 + t27 ^ 2 / 0.2e1, t72 ^ 2 * t109, -t67 * t13 * t72, -t11 * t106, t67 ^ 2 * t109, t11 * t107, t11 ^ 2 / 0.2e1, -t1 * t11 + t3 * t107, t3 * t106 + t2 * t11, (-t1 * t72 - t2 * t67) * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t98, t26 ^ 2 * t108, -t25 * t105, t26 * t91, t25 ^ 2 * t108, -t25 * t91, t61, pkin(14) * t78 * t25, pkin(14) * t105, 0, pkin(14) ^ 2 * t108, t35 ^ 2 * t108, t35 * t78 * t32, 0, t32 ^ 2 * t108, 0, 0, -t32 * t104, t35 * t104, (t20 * t32 + t21 * t35) * t86, t39 + (t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1) * t101, t17 ^ 2 / 0.2e1, t18 * t17, t60 * t17, t18 ^ 2 / 0.2e1, t18 * t60, t52, t18 * t94, -t17 * t94, 0, t57, t49, -t90, t50, t48, -t85, t61, t23 * qJD(2) + t69 * t93, -t24 * qJD(2) + t74 * t93, (-t23 * t74 - t24 * t69) * qJD(1), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t55 ^ 2 * t108, t55 * t103, t55 * t92, t56 ^ 2 * t108, t56 * t92, t52, -t10 * t103 + t40, t10 * t78 * t55 + t83, (-t55 * t68 - t56 * t73) * t86, t10 ^ 2 * t108 + t100;];
T_reg = t4;
