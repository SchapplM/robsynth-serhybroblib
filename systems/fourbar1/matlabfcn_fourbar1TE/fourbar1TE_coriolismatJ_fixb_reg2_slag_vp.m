% Calculate inertial parameters regressor of coriolis matrix for
% fourbar1TE
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
% cmat_reg [(1*1)x(1*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1TE_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:41
% EndTime: 2020-04-24 19:49:42
% DurationCPUTime: 0.29s
% Computational Cost: add. (831->78), mult. (1331->146), div. (36->9), fcn. (246->4), ass. (0->70)
t54 = pkin(2) ^ 2;
t57 = pkin(1) ^ 2;
t45 = cos(qJ(1));
t99 = pkin(1) * t45;
t82 = pkin(2) * t99;
t86 = -0.2e1 * t82 + t57;
t27 = t54 + t86;
t95 = 0.1e1 / t27 ^ 2 * t54;
t49 = pkin(4) ^ 2;
t50 = pkin(3) ^ 2;
t33 = t49 + t50;
t46 = pkin(3) + pkin(4);
t47 = pkin(3) - pkin(4);
t90 = (t47 ^ 2 * t46 ^ 2);
t102 = 0.2e1 * t27 * t33 - (2 * t90);
t44 = sin(qJ(1));
t100 = pkin(1) * t44;
t101 = pkin(2) * t100;
t93 = (pkin(2) + t46) * (pkin(2) - t46);
t21 = t86 + t93;
t92 = (pkin(2) + t47) * (pkin(2) - t47);
t22 = t86 + t92;
t96 = t21 * t22;
t58 = sqrt(-t96);
t12 = 0.1e1 / t58;
t79 = t12 * t95;
t76 = t101 * t95;
t9 = (-t21 - t22) * t101;
t98 = t12 * t9;
t48 = 0.2e1 * t54;
t52 = t54 ^ 2;
t55 = t57 ^ 2;
t97 = (t55 + (t48 - 0.2e1 / 0.3e1 * t50 - 0.2e1 / 0.3e1 * t49) * t57 + t52 + (-0.4e1 / 0.3e1 * t50 - 0.4e1 / 0.3e1 * t49) * t54 + t90 / 0.3e1) * pkin(2);
t31 = -pkin(2) + t99;
t94 = t31 * t44;
t91 = t44 ^ 2 * t57;
t89 = t44 * t58;
t88 = t45 * t57;
t87 = t46 * t47;
t85 = -t49 / 0.2e1 + t54;
t84 = t49 - t50;
t83 = qJD(1) / pkin(3);
t3 = t31 * t98;
t80 = pkin(2) * t91;
t78 = t58 * t87;
t29 = t50 / 0.2e1 + t85;
t77 = t54 - t84;
t75 = 0.1e1 / t27 * t76;
t74 = 0.1e1 / t96 * t9 * t79;
t39 = t45 ^ 2;
t59 = pkin(1) * t57;
t73 = t59 * (-t50 / 0.2e1 + t85) * pkin(2) * t39;
t72 = 0.2e1 * t75;
t71 = t12 * t75;
t70 = 0.1e1 / t22 ^ 2 * t76;
t19 = 0.1e1 / t22;
t69 = 0.1e1 / t21 ^ 2 * t19 * t76;
t68 = t31 * t99 - t91;
t28 = t57 + t29;
t24 = t27 + t84;
t23 = t77 + t86;
t17 = 0.1e1 / t21;
t16 = (t54 - t50 / 0.6e1 - t49 / 0.6e1) * t57 + t52 + (-0.5e1 / 0.6e1 * t50 - 0.5e1 / 0.6e1 * t49) * t54 + t90 / 0.6e1;
t10 = t31 * t58;
t8 = -0.4e1 * t28 * t82 + t55 + t54 * t77 + (0.4e1 * t29 * t39 + t48 + t84) * t57;
t7 = t24 * t100 + t10;
t6 = -t23 * t100 + t10;
t5 = t7 ^ 2;
t4 = t6 ^ 2;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t5 * t69 + (t5 * t70 + (-t7 * (0.2e1 * t80 + t3 + (t24 * t45 - t89) * pkin(1)) * t95 + t5 * t72) * t19) * t17) * qJD(1), (-(t8 * t98 + 0.4e1 * t31 * t33 * t80 + (0.4e1 * pkin(1) * pkin(2) * t28 - 0.8e1 * t29 * t88) * t89 + t68 * t102) * t79 / 0.2e1 + (-t74 / 0.2e1 + 0.2e1 * t71) * (pkin(1) * t94 * t102 + t58 * t8)) * t83, (-(t68 * t78 + (0.12e2 * t73 - 0.12e2 * t16 * t88 + (t3 * t87 + 0.3e1 * t97) * pkin(1)) * t44) * t79 + (0.4e1 * t71 - t74) * (-0.4e1 * t45 * t73 + t59 ^ 2 / 0.2e1 + (0.6e1 * t16 * t39 + 0.3e1 / 0.2e1 * t52 - t90 / 0.2e1) * t57 + (0.3e1 / 0.2e1 * t55 - t33 * t57 + t92 * t93 / 0.2e1) * t54 + (-0.3e1 * t45 * t97 + t78 * t94) * pkin(1))) * t83, 0, 0, 0, 0, 0, 0, 0, (t4 * t69 + (t4 * t70 + (-t6 * (-0.2e1 * t80 + t3 + (-t23 * t45 - t89) * pkin(1)) * t95 + t4 * t72) * t19) * t17) * qJD(1), 0, 0, 0, 0;];
cmat_reg = t1;
