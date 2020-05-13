% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [1x(1*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbarprisDE1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisDE1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_invdynJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:10:10
% EndTime: 2020-05-07 09:10:15
% DurationCPUTime: 0.45s
% Computational Cost: add. (1389->71), mult. (770->152), div. (224->15), fcn. (38->2), ass. (0->80)
t40 = (qJ(1) + pkin(3));
t80 = -pkin(2) - t40;
t32 = pkin(1) + t80;
t95 = 1 / t32;
t105 = 2 * t95;
t29 = pkin(1) - t80;
t79 = -pkin(2) + t40;
t30 = pkin(1) - t79;
t31 = pkin(1) + t79;
t88 = t31 * t32;
t75 = t30 * t88;
t66 = t29 * t75;
t49 = sqrt(-t66);
t11 = 0.1e1 / t49;
t104 = t11 * (-t75 + (t88 + (t31 - t32) * t30) * t29);
t103 = 2 * g(1);
t102 = 0.1e1 / t66 * t104;
t19 = 0.1e1 / t29;
t23 = 0.1e1 / t31;
t101 = t19 * t23;
t36 = (t40 ^ 2);
t37 = 1 / t40;
t43 = (qJ(1) ^ 2);
t44 = (pkin(3) ^ 2);
t90 = (pkin(3) * qJ(1));
t73 = -2 * t90 - t43 - t44;
t45 = pkin(2) ^ 2;
t47 = (pkin(1) ^ 2);
t83 = (t47 - t45);
t16 = t73 + t83;
t97 = t16 ^ 2;
t100 = t37 / t36 * t97;
t96 = 0.1e1 / t30;
t99 = t95 * t96;
t81 = qJDD(1) * t97;
t98 = t97 * t99;
t22 = 0.1e1 / t30 ^ 2;
t26 = 0.1e1 / t32 ^ 2;
t48 = 0.1e1 / pkin(1);
t94 = -t48 / 0.2e1;
t93 = t48 / 0.2e1;
t15 = -t73 + t83;
t72 = t104 / 0.2e1;
t4 = g(2) * t72;
t58 = ((t103 * t40) + t4) * t93;
t38 = 1 / t40 ^ 2;
t71 = t38 * t94;
t9 = t49 * g(2);
t92 = t37 * t58 + ((g(1) * t15) + t9) * t71;
t91 = (t44 * g(1)) + t9;
t42 = (qJD(1) ^ 2);
t89 = t16 * t42;
t87 = t36 * t42;
t86 = t38 * t42;
t85 = t40 * t42;
t84 = (-t45 - t47);
t18 = t40 * qJD(1);
t82 = qJD(1) * t18;
t78 = t97 * t86;
t77 = t42 * t100;
t76 = t95 * t87;
t74 = 4 * t82;
t70 = t95 * t78;
t69 = t26 * t78;
t67 = t23 * t76;
t65 = t23 * t70;
t24 = 0.1e1 / t31 ^ 2;
t64 = t24 * t70;
t63 = t43 * t99 * t101;
t62 = t96 * t65;
t61 = t22 * t65;
t60 = t38 * t63;
t59 = t42 * t63;
t3 = g(1) * t72 - 0.2e1 * g(2) * t40;
t20 = 0.1e1 / t29 ^ 2;
t57 = t20 * t62;
t46 = 0.1e1 / pkin(2);
t8 = g(1) * t49 - g(2) * t15;
t1 = t57 / 0.2e1 + (-t61 / 0.2e1 + (t64 / 0.2e1 + (-t69 / 0.2e1 + ((t77 + (-t81 + 2 * (2 * t82 - t85) * t16) * t38) * t95)) * t23) * t96) * t19;
t2 = [0, 0, 0, 0, 0, t1, t92, (-t3 * t37 / 0.2e1 + t8 * t38 / 0.2e1) * t48, 0, 0, 0, 0, 0, t1, 0, 0, (-t11 * t89 + t8 * t94) * t38 + (-t89 * t102 / 0.2e1 + t3 * t93 + 0.2e1 * (qJDD(1) * t16 + t85) * t11 + (qJD(1) * t102 * t16 - 0.4e1 * t18 * t11) * qJD(1)) * t37, 0, -t19 * t62 + (t57 + (-t61 + (t64 + (-t69 + ((t77 + (-t81 + (t74 - 2 * t85) * t16) * t38) * t105)) * t23) * t96) * t19) * qJ(1) + t92, -t60 * t81 + t59 * t100 + qJDD(1) + (qJ(1) * t4 + ((3 * t43 + t84 + 4 * t90) * g(1)) + t91) * t37 * t93 + (((t43 - t47) * pkin(3) * t103) + (t91 + ((t43 + t84) * g(1))) * qJ(1)) * t71 + (-0.2e1 * t37 * t59 + t60 * t74) * t16 + ((t24 * t19 + t23 * t20) * t43 * t98 / 0.2e1 + (-qJ(1) * t98 - (t95 * t22 + t26 * t96) * t97 * t43 / 0.2e1) * t101) * t86, 0, 0, 0, 0, 0, 0.2e1 * t96 * t20 * t67 + 0.2e1 * (-t22 * t67 + (t24 * t76 + (-t26 * t87 + ((-qJDD(1) * t36 - t85) * t105)) * t23) * t96) * t19, t46 * t58, t3 * t46 * t94, 0, 0;];
tau_reg = t2;
