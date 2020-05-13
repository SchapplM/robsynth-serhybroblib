% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = fourbar1turnOL_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:05
% EndTime: 2020-04-12 19:41:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (409->109), mult. (989->156), div. (0->0), fcn. (718->8), ass. (0->75)
t62 = cos(qJ(3));
t51 = qJDD(2) + qJDD(3);
t63 = cos(qJ(2));
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t76 = t59 * t58;
t27 = (-t63 * t62 + t76) * qJD(1);
t28 = (-t63 * t58 - t59 * t62) * qJD(1);
t82 = t27 * t28;
t9 = t51 + t82;
t85 = t62 * t9;
t46 = t63 * qJDD(1);
t72 = qJD(1) * qJD(2);
t70 = t59 * t72;
t36 = t46 - t70;
t84 = t36 - t70;
t25 = t27 ^ 2;
t26 = t28 ^ 2;
t52 = qJD(2) + qJD(3);
t50 = t52 ^ 2;
t83 = cos(qJ(1));
t81 = t52 * t27;
t80 = t52 * t28;
t56 = t63 ^ 2;
t66 = qJD(1) ^ 2;
t79 = t56 * t66;
t10 = -t82 + t51;
t78 = t58 * t10;
t77 = t58 * t27;
t60 = sin(qJ(1));
t37 = t60 * g(1) - t83 * g(2);
t14 = t84 * pkin(2) + t37;
t75 = t62 * t14;
t74 = t62 * t28;
t41 = t63 * t66 * t59;
t73 = qJDD(2) + t41;
t71 = qJD(1) * qJD(4);
t69 = t63 * t72;
t38 = t83 * g(1) + t60 * g(2);
t23 = t63 * g(3) - t59 * t38;
t15 = t73 * pkin(2) - t23;
t24 = -t59 * g(3) - t63 * t38;
t65 = qJD(2) ^ 2;
t16 = (-t65 - t79) * pkin(2) + t24;
t5 = -t62 * t15 + t58 * t16;
t44 = t59 * qJDD(1);
t34 = t44 + t69;
t68 = t58 * t34 - t62 * t36;
t6 = t58 * t15 + t62 * t16;
t67 = t62 * t5 - t58 * t6;
t57 = sin(qJ(4));
t43 = t57 * qJDD(1);
t61 = cos(qJ(4));
t33 = 0.2e1 * t61 * t71 + t43;
t45 = t61 * qJDD(1);
t35 = -0.2e1 * t57 * t71 + t45;
t8 = t27 * qJD(3) - t62 * t34 - t58 * t36;
t7 = -t28 * qJD(3) + t68;
t4 = t8 + t81;
t64 = qJD(4) ^ 2;
t54 = t59 ^ 2;
t48 = t61 ^ 2 * t66;
t47 = t57 ^ 2 * t66;
t39 = t61 * t66 * t57;
t32 = -t66 * pkin(1) - t38;
t31 = qJDD(1) * pkin(1) + t37;
t20 = -t57 * g(3) + t61 * t32;
t19 = t61 * g(3) + t57 * t32;
t18 = -t26 + t50;
t17 = t25 - t50;
t11 = t26 - t25;
t3 = -t8 + t81;
t2 = -t80 + t7;
t1 = t80 + t7;
t12 = [0, 0, 0, 0, 0, qJDD(1), t37, t38, 0, 0, (t34 + t69) * t59, t59 * (t46 - 0.2e1 * t70) + t63 * (t44 + 0.2e1 * t69), t59 * t73 + t63 * (-t54 * t66 + t65), t84 * t63, t59 * (-t65 + t79) + t63 * (qJDD(2) - t41), 0, t63 * t37, -t59 * t37, t59 * t23 + t63 * t24, 0, t59 * (t58 * t80 - t62 * t8) + t63 * (-t52 * t74 - t58 * t8), t59 * (-t62 * t2 + t58 * t4) + t63 * (-t58 * t2 - t62 * t4), t59 * (t58 * t18 - t85) + t63 * (-t62 * t18 - t58 * t9), t59 * (t58 * t7 + t62 * t81) + t63 * (t52 * t77 - t62 * t7), t59 * (-t62 * t17 + t78) + t63 * (-t62 * t10 - t58 * t17), (t59 * (-t27 * t62 - t28 * t58) + t63 * (t74 - t77)) * t52, t14 * t76 + t63 * (pkin(2) * t2 - t75), t59 * t75 + t63 * (-pkin(2) * t4 + t58 * t14), t59 * t67 + t63 * (t58 * t5 + t62 * t6 - pkin(2) * (-t25 - t26)), t63 * pkin(2) * t14, t33 * t57, t61 * t33 + t57 * t35, t57 * (qJDD(4) + t39) + t61 * (-t47 + t64), t35 * t61, t57 * (t48 - t64) + t61 * (qJDD(4) - t39), 0, pkin(1) * t35 + t61 * t31, -pkin(1) * t33 - t57 * t31, t57 * t19 + t61 * t20 + pkin(1) * (t47 + t48), pkin(1) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, (t54 - t56) * t66, t44, t41, t46, qJDD(2), -t23, -t24, 0, 0, -t82, t11, -t3, t82, t1, t51, pkin(2) * (-t58 * (-t50 - t25) - t85) + t5, pkin(2) * (t78 - t62 * (-t26 - t50)) + t6, pkin(2) * (-t58 * (-(qJD(3) - t52) * t28 + t68) - t62 * t3), -pkin(2) * t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t11, -t3, t82, t1, t51, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t47 - t48, t43, t39, t45, qJDD(4), -t19, -t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauJ_reg = t12;
