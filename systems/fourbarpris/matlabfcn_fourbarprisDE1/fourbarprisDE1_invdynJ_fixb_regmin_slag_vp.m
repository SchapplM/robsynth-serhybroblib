% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% tau_reg [1x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbarprisDE1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_invdynJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:53
% EndTime: 2020-06-27 17:14:57
% DurationCPUTime: 0.46s
% Computational Cost: add. (1244->71), mult. (676->152), div. (184->15), fcn. (38->2), ass. (0->79)
t39 = (qJ(1) + pkin(3));
t79 = -pkin(2) - t39;
t31 = pkin(1) + t79;
t94 = 1 / t31;
t104 = 2 * t94;
t28 = pkin(1) - t79;
t78 = -pkin(2) + t39;
t29 = pkin(1) - t78;
t30 = pkin(1) + t78;
t87 = t30 * t31;
t74 = t29 * t87;
t65 = t28 * t74;
t48 = sqrt(-t65);
t10 = 0.1e1 / t48;
t103 = t10 * (-t74 + (t87 + (t30 - t31) * t29) * t28);
t102 = 2 * g(1);
t101 = 0.1e1 / t65 * t103;
t18 = 0.1e1 / t28;
t22 = 0.1e1 / t30;
t100 = t18 * t22;
t35 = (t39 ^ 2);
t36 = 1 / t39;
t42 = (qJ(1) ^ 2);
t43 = (pkin(3) ^ 2);
t89 = (pkin(3) * qJ(1));
t72 = -2 * t89 - t42 - t43;
t44 = pkin(2) ^ 2;
t46 = (pkin(1) ^ 2);
t82 = (t46 - t44);
t15 = t72 + t82;
t96 = t15 ^ 2;
t99 = t36 / t35 * t96;
t95 = 0.1e1 / t29;
t98 = t94 * t95;
t80 = qJDD(1) * t96;
t97 = t96 * t98;
t21 = 0.1e1 / t29 ^ 2;
t25 = 0.1e1 / t31 ^ 2;
t47 = 0.1e1 / pkin(1);
t93 = -t47 / 0.2e1;
t92 = t47 / 0.2e1;
t14 = -t72 + t82;
t71 = t103 / 0.2e1;
t3 = g(2) * t71;
t57 = ((t39 * t102) + t3) * t92;
t37 = 1 / t39 ^ 2;
t70 = t37 * t93;
t8 = t48 * g(2);
t91 = t36 * t57 + ((g(1) * t14) + t8) * t70;
t90 = (t43 * g(1)) + t8;
t41 = (qJD(1) ^ 2);
t88 = t15 * t41;
t86 = t35 * t41;
t85 = t37 * t41;
t84 = t39 * t41;
t83 = (-t44 - t46);
t17 = t39 * qJD(1);
t81 = qJD(1) * t17;
t77 = t96 * t85;
t76 = t41 * t99;
t75 = t94 * t86;
t73 = 4 * t81;
t69 = t94 * t77;
t68 = t25 * t77;
t66 = t22 * t75;
t64 = t22 * t69;
t23 = 0.1e1 / t30 ^ 2;
t63 = t23 * t69;
t62 = t42 * t98 * t100;
t61 = t95 * t64;
t60 = t21 * t64;
t59 = t37 * t62;
t58 = t41 * t62;
t2 = g(1) * t71 - 0.2e1 * g(2) * t39;
t19 = 0.1e1 / t28 ^ 2;
t56 = t19 * t61;
t45 = 0.1e1 / pkin(2);
t7 = g(1) * t48 - g(2) * t14;
t1 = [t56 / 0.2e1 + (-t60 / 0.2e1 + (t63 / 0.2e1 + (-t68 / 0.2e1 + ((t76 + (-t80 + 2 * (2 * t81 - t84) * t15) * t37) * t94)) * t22) * t95) * t18, t91, (-t2 * t36 / 0.2e1 + t7 * t37 / 0.2e1) * t47, (-t10 * t88 + t7 * t93) * t37 + (-t88 * t101 / 0.2e1 + t2 * t92 + 0.2e1 * (qJDD(1) * t15 + t84) * t10 + (t15 * qJD(1) * t101 - 0.4e1 * t17 * t10) * qJD(1)) * t36, -t18 * t61 + (t56 + (-t60 + (t63 + (-t68 + ((t76 + (-t80 + (t73 - 2 * t84) * t15) * t37) * t104)) * t22) * t95) * t18) * qJ(1) + t91, -t59 * t80 + t58 * t99 + qJDD(1) + (qJ(1) * t3 + ((3 * t42 + t83 + 4 * t89) * g(1)) + t90) * t36 * t92 + (((t42 - t46) * pkin(3) * t102) + (t90 + ((t42 + t83) * g(1))) * qJ(1)) * t70 + (-0.2e1 * t36 * t58 + t59 * t73) * t15 + ((t23 * t18 + t22 * t19) * t42 * t97 / 0.2e1 + (-qJ(1) * t97 - (t94 * t21 + t25 * t95) * t96 * t42 / 0.2e1) * t100) * t85, 0.2e1 * t95 * t19 * t66 + 0.2e1 * (-t21 * t66 + (t23 * t75 + (-t25 * t86 + ((-qJDD(1) * t35 - t84) * t104)) * t22) * t95) * t18, t45 * t57, t2 * t45 * t93;];
tau_reg = t1;
