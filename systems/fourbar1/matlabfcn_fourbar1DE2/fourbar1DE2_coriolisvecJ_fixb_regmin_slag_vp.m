% Calculate minimal parameter regressor of coriolis joint torque vector for
% fourbar1DE2
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
% tauc_reg [1x9]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:38
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1DE2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:38:29
% EndTime: 2020-06-26 17:38:31
% DurationCPUTime: 0.52s
% Computational Cost: add. (1339->91), mult. (2202->182), div. (50->9), fcn. (396->4), ass. (0->88)
t48 = pkin(3) + pkin(4);
t49 = pkin(3) - pkin(4);
t100 = (t49 ^ 2 * t48 ^ 2);
t57 = pkin(2) ^ 2;
t60 = pkin(1) ^ 2;
t47 = cos(qJ(1));
t107 = pkin(2) * t47;
t92 = pkin(1) * t107;
t95 = -0.2e1 * t92 + t60;
t29 = t57 + t95;
t52 = pkin(4) ^ 2;
t53 = pkin(3) ^ 2;
t35 = t52 + t53;
t113 = 0.2e1 * t29 * t35 - (2 * t100);
t46 = sin(qJ(1));
t109 = pkin(1) * t46;
t112 = pkin(2) * t109;
t102 = (pkin(2) + t48) * (pkin(2) - t48);
t23 = t95 + t102;
t101 = (pkin(2) + t49) * (pkin(2) - t49);
t24 = t95 + t101;
t104 = t23 * t24;
t61 = sqrt(-t104);
t14 = 0.1e1 / t61;
t27 = 0.1e1 / t29 ^ 2;
t50 = 0.2e1 * t57;
t55 = t57 ^ 2;
t58 = t60 ^ 2;
t17 = t58 + (t50 - 0.2e1 / 0.3e1 * t53 - 0.2e1 / 0.3e1 * t52) * t60 + t55 + (-0.4e1 / 0.3e1 * t53 - 0.4e1 / 0.3e1 * t52) * t57 + t100 / 0.3e1;
t111 = 0.3e1 * t17;
t110 = -2 * qJD(1);
t108 = pkin(1) * t47;
t11 = (-t23 - t24) * t112;
t106 = t11 * t14;
t33 = -pkin(2) + t108;
t105 = t14 * t33;
t103 = t27 * t57;
t99 = t46 * t61;
t98 = t47 * t60;
t97 = t48 * t49;
t96 = t60 * t46 ^ 2;
t94 = -t52 / 0.2e1 + t57;
t93 = t52 - t53;
t91 = t57 * t112;
t90 = pkin(2) * t96;
t41 = t47 ^ 2;
t62 = pkin(1) * t60;
t89 = (-t53 / 0.2e1 + t94) * t41 * t62;
t88 = t33 * t97;
t31 = t53 / 0.2e1 + t94;
t87 = t57 - t93;
t86 = qJD(1) * t103;
t51 = qJD(1) ^ 2;
t85 = t51 * t91;
t18 = (t57 - t53 / 0.6e1 - t52 / 0.6e1) * t60 + t55 + (-0.5e1 / 0.6e1 * t53 - 0.5e1 / 0.6e1 * t52) * t57 + t100 / 0.6e1;
t84 = -0.12e2 * t18 * t98;
t83 = 0.12e2 * t89;
t82 = t46 * t88;
t54 = 0.1e1 / pkin(3);
t81 = t54 * t86;
t80 = t27 * t85;
t79 = t51 * t54 * t103 / 0.2e1;
t28 = 0.1e1 / t29 * t27;
t78 = 0.2e1 * t28 * t85;
t77 = 0.1e1 / t24 ^ 2 * t80;
t76 = t33 * t108 - t96;
t21 = 0.1e1 / t24;
t75 = 0.1e1 / t23 ^ 2 * t21 * t80;
t74 = t76 * t61 * t97;
t25 = t87 + t95;
t73 = -0.2e1 * t90 + (-t25 * t47 - t99) * pkin(1);
t26 = t29 + t93;
t72 = 0.2e1 * t90 + (t26 * t47 - t99) * pkin(1);
t30 = t60 + t31;
t71 = 0.4e1 * t33 * t35 * t90 + (0.4e1 * pkin(1) * pkin(2) * t30 - 0.8e1 * t31 * t98) * t99 + t76 * t113;
t19 = 0.1e1 / t23;
t15 = t14 / t104;
t12 = t33 * t61;
t10 = -0.4e1 * t30 * t92 + t58 + t57 * t87 + (0.4e1 * t31 * t41 + t50 + t93) * t60;
t9 = qJD(1) * t11;
t8 = t26 * t109 + t12;
t7 = -t25 * t109 + t12;
t6 = t8 ^ 2;
t5 = t7 ^ 2;
t4 = t11 * t105;
t3 = t9 * t105;
t1 = -0.4e1 * t89 * t107 + t62 ^ 2 / 0.2e1 + (0.6e1 * t18 * t41 + 0.3e1 / 0.2e1 * t55 - t100 / 0.2e1) * t60 + (0.3e1 / 0.2e1 * t58 - t35 * t60 + t101 * t102 / 0.2e1) * t57 + (-0.3e1 * t17 * t107 + t61 * t82) * pkin(1);
t2 = [0, 0, 0, t6 * t75 + (t6 * t77 + (t6 * t78 + ((t72 * qJD(1) + t3) * t110 + t51 * (t4 + t72)) * t8 * t103) * t21) * t19, ((t10 * t106 + t71) * t79 + (-t10 * t14 * t9 - t71 * qJD(1)) * t81) * t14 + (t14 * t54 * t78 + (t11 * t79 - t81 * t9) * t15) * (t33 * t109 * t113 + t10 * t61), ((t14 * (t74 + (pkin(2) * t83 + t84 + (pkin(2) * t111 + t88 * t106) * pkin(1)) * t46) * t103 + (t15 * t11 * t103 + 0.4e1 * t14 * t28 * t91) * t1) * t51 + 0.2e1 * (-t14 * (t74 + (t84 + (pkin(1) * t111 + t83) * pkin(2)) * t46) * qJD(1) + (-pkin(1) * t14 ^ 2 * t82 - t1 * t15) * t9) * t86) * t54, t5 * t75 + (t5 * t77 + (t5 * t78 + ((t73 * qJD(1) + t3) * t110 + t51 * (t4 + t73)) * t7 * t103) * t21) * t19, 0, 0;];
tauc_reg = t2;
