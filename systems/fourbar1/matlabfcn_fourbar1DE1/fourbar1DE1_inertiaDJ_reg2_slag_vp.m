% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% MMD_reg [((1+1)*1/2)x(1*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = fourbar1DE1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_inertiaDJ_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_inertiaDJ_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_inertiaDJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:21
% EndTime: 2020-04-24 19:57:23
% DurationCPUTime: 0.40s
% Computational Cost: add. (831->79), mult. (1368->157), div. (36->9), fcn. (246->4), ass. (0->77)
t54 = pkin(2) ^ 2;
t57 = pkin(1) ^ 2;
t45 = cos(qJ(1));
t98 = pkin(2) * t45;
t84 = pkin(1) * t98;
t89 = -0.2e1 * t84 + t57;
t27 = t54 + t89;
t96 = 0.1e1 / t27 ^ 2 * t54;
t102 = pkin(2) * t96;
t46 = pkin(3) + pkin(4);
t95 = (pkin(2) + t46) * (pkin(2) - t46);
t21 = t89 + t95;
t47 = pkin(3) - pkin(4);
t94 = (pkin(2) + t47) * (pkin(2) - t47);
t22 = t89 + t94;
t97 = t21 * t22;
t58 = sqrt(-t97);
t12 = 0.1e1 / t58;
t101 = t12 ^ 2;
t44 = sin(qJ(1));
t100 = pkin(1) * t44;
t99 = pkin(1) * t45;
t93 = t44 ^ 2 * t57;
t92 = t47 ^ 2 * t46 ^ 2;
t91 = t44 * t58;
t90 = t46 * t47;
t49 = pkin(4) ^ 2;
t88 = -t49 / 0.2e1 + t54;
t50 = pkin(3) ^ 2;
t87 = t49 - t50;
t86 = qJD(1) * t12;
t85 = pkin(2) * t100;
t83 = -0.2e1 * t96;
t82 = pkin(2) * t93;
t81 = t44 * t45 * t57;
t80 = t58 * t90;
t29 = t50 / 0.2e1 + t88;
t79 = t54 - t87;
t78 = qJD(1) * t100;
t33 = t49 + t50;
t14 = t27 * t33 - t92;
t31 = -pkin(2) + t99;
t77 = 0.2e1 * t14 * pkin(1) * t31;
t8 = (-t21 - t22) * pkin(2) * t78;
t76 = 0.2e1 * t8 * t96;
t75 = 0.1e1 / t27 * t100 * t102;
t39 = t45 ^ 2;
t59 = pkin(1) * t57;
t74 = pkin(2) * (-t50 / 0.2e1 + t88) * t39 * t59;
t73 = t31 * t80;
t72 = 0.4e1 * t75;
t71 = 0.2e1 * t78 * t102;
t70 = qJD(1) * t72;
t69 = 0.1e1 / t22 ^ 2 * t71;
t19 = 0.1e1 / t22;
t68 = t19 / t21 ^ 2 * t71;
t55 = t57 ^ 2;
t52 = t54 ^ 2;
t51 = 0.1e1 / pkin(3);
t48 = 0.2e1 * t54;
t28 = t57 + t29;
t24 = t27 + t87;
t23 = t79 + t89;
t17 = 0.1e1 / t21;
t16 = (t54 - t50 / 0.6e1 - t49 / 0.6e1) * t57 + t52 + (-0.5e1 / 0.6e1 * t50 - 0.5e1 / 0.6e1 * t49) * t54 + t92 / 0.6e1;
t15 = t55 + (t48 - 0.2e1 / 0.3e1 * t50 - 0.2e1 / 0.3e1 * t49) * t57 + t52 + (-0.4e1 / 0.3e1 * t50 - 0.4e1 / 0.3e1 * t49) * t54 + t92 / 0.3e1;
t13 = t12 / t97;
t10 = t31 * t58;
t9 = -0.4e1 * t28 * t84 + t55 + t54 * t79 + (0.4e1 * t29 * t39 + t48 + t87) * t57;
t7 = t100 * t24 + t10;
t6 = -t100 * t23 + t10;
t5 = t7 ^ 2;
t4 = t6 ^ 2;
t3 = t31 * t12 * t8;
t2 = t44 * t77 + t58 * t9;
t1 = -0.4e1 * t45 * t74 + t59 ^ 2 / 0.2e1 + (0.6e1 * t16 * t39 + 0.3e1 / 0.2e1 * t52 - t92 / 0.2e1) * t57 + (0.3e1 / 0.2e1 * t55 - t33 * t57 + t94 * t95 / 0.2e1) * t54 + (-0.3e1 * t15 * t98 + t44 * t73) * pkin(1);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t68 + (t5 * t69 + (t7 * (t3 + (0.2e1 * t82 + (t24 * t45 - t91) * pkin(1)) * qJD(1)) * t83 + t5 * t70) * t19) * t17, ((-t9 * t101 / 0.2e1 - t2 * t13 / 0.2e1) * t76 + (-((0.4e1 * t28 * t85 - 0.8e1 * t29 * t81) * t58 + t45 * t77 + 0.4e1 * t33 * t31 * t82 - 0.2e1 * t14 * t93) * t96 + t2 * t72) * t86) * t51, ((-t100 * t101 * t31 * t90 - t1 * t13) * t76 + ((0.3e1 * t15 * t85 - 0.12e2 * t16 * t81 + 0.12e2 * t44 * t74 + t73 * t99 - t80 * t93) * t83 + 0.8e1 * t1 * t75) * t86) * t51, 0, 0, 0, 0, 0, 0, 0, t4 * t68 + (t4 * t69 + (t6 * (t3 + (-0.2e1 * t82 + (-t23 * t45 - t91) * pkin(1)) * qJD(1)) * t83 + t4 * t70) * t19) * t17, 0, 0, 0, 0;];
MMD_reg = t11;
