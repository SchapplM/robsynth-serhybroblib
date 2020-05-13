% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% fourbar1turnIC
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
% tau_reg [2x(5*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnIC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:34
% EndTime: 2020-05-07 11:33:35
% DurationCPUTime: 0.50s
% Computational Cost: add. (365->101), mult. (860->194), div. (34->3), fcn. (653->12), ass. (0->79)
t89 = -qJ(4) + qJ(2);
t28 = sin(qJ(3) + t89);
t25 = 0.1e1 / t28;
t80 = (-pkin(3) * t28 + pkin(2) * sin(t89)) / pkin(3) * t25;
t69 = 0.1e1 + t80;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t44 = cos(qJ(3));
t45 = cos(qJ(2));
t13 = t40 * t45 + t44 * t41;
t10 = qJD(1) * t13;
t93 = qJD(1) * t45;
t94 = qJD(1) * t41;
t9 = t40 * t94 - t44 * t93;
t62 = t69 * t9;
t106 = t10 * t62;
t82 = qJD(2) * qJD(1);
t73 = t45 * t82;
t81 = qJD(3) * qJD(1);
t87 = t41 * qJDD(1);
t105 = t45 * t81 + t73 + t87;
t75 = t41 * t82;
t85 = t45 * qJDD(1);
t104 = t41 * t81 + t75 - t85;
t2 = t104 * t44 + t105 * t40;
t33 = qJD(2) + qJD(3);
t103 = -t10 * t33 + t2;
t42 = sin(qJ(1));
t46 = cos(qJ(1));
t14 = g(1) * t46 + g(2) * t42;
t100 = pkin(2) * t40;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t98 = t39 * t43;
t34 = t39 ^ 2;
t36 = t43 ^ 2;
t97 = t34 - t36;
t35 = t41 ^ 2;
t37 = t45 ^ 2;
t96 = t35 - t37;
t95 = pkin(2) * qJD(2);
t92 = qJD(2) * t41;
t91 = qJD(3) * t33;
t90 = qJDD(1) * pkin(1);
t88 = t37 * qJDD(1);
t86 = t44 * qJDD(2);
t83 = qJD(1) * qJD(4);
t49 = qJD(1) ^ 2;
t79 = t41 * t49 * t45;
t32 = qJDD(2) + qJDD(3);
t78 = pkin(2) * t93;
t77 = t44 * t95;
t76 = qJD(3) * t95;
t71 = t25 / pkin(4) * t100;
t70 = -0.2e1 * pkin(1) * t83;
t68 = t83 * t98;
t67 = t41 * t73;
t66 = t49 * t71;
t65 = g(1) * t42 - g(2) * t46;
t64 = qJDD(1) * t71;
t12 = t40 * t41 - t44 * t45;
t38 = qJ(2) + qJ(3);
t30 = sin(t38);
t63 = -g(3) * t30 + qJDD(2) * t100 + t44 * t76 + t9 * t78;
t61 = t104 * t40;
t60 = t66 * t98;
t59 = pkin(1) * t49 + t14;
t58 = t65 + 0.2e1 * t90;
t57 = -t9 * t33 + t61;
t56 = -g(3) * t45 + t14 * t41;
t31 = cos(t38);
t55 = g(3) * t31 - t10 * t78 - t14 * t30 + t40 * t76;
t52 = pkin(2) ^ 2;
t48 = qJD(2) ^ 2;
t47 = qJD(4) ^ 2;
t4 = t33 * t13;
t3 = t33 * t12;
t1 = (-t33 * t93 - t87) * t44 + t61;
t5 = [0, 0, 0, 0, 0, qJDD(1), t65, t14, 0, 0, t35 * qJDD(1) + 0.2e1 * t67, 0.2e1 * t41 * t85 - 0.2e1 * t96 * t82, qJDD(2) * t41 + t48 * t45, -0.2e1 * t67 + t88, qJDD(2) * t45 - t48 * t41, 0, t65 * t45, -t65 * t41, -t14, 0, -t1 * t13 - t10 * t3, t1 * t12 - t10 * t4 - t13 * t2 + t3 * t9, -t13 * t32 + t3 * t33, t2 * t12 + t9 * t4, t12 * t32 + t4 * t33, 0, -t65 * t31 + ((-qJD(1) * t12 - t9) * t92 + (qJD(1) * t4 + qJDD(1) * t12 + t2) * t45) * pkin(2), t65 * t30 + (-0.2e1 * t10 * t92 + (-qJD(1) * t3 + qJDD(1) * t13 - t1) * t45) * pkin(2), ((-t12 * t40 - t13 * t44) * qJDD(2) + (t3 * t44 - t4 * t40 + (-t12 * t44 + t13 * t40) * qJD(3)) * qJD(2)) * pkin(2) - t14, t52 * t88 + (t65 * pkin(2) - 0.2e1 * t52 * t75) * t45, t34 * qJDD(1) + 0.2e1 * t68, 0.2e1 * qJDD(1) * t98 - 0.2e1 * t97 * t83, qJDD(4) * t39 + t47 * t43, t36 * qJDD(1) - 0.2e1 * t68, qJDD(4) * t43 - t47 * t39, 0, t39 * t70 + t58 * t43, -t58 * t39 + t43 * t70, -t14, (t65 + t90) * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t96 * t49, t87, t79, t85, qJDD(2), t56, g(3) * t41 + t14 * t45, 0, 0, t106, t69 * (t10 ^ 2 - t9 ^ 2), -t69 * t44 * t105 + t57 * t80 + t57, -t106, t69 * t103, t69 * t32, t55 * t80 + (t40 * t91 - t44 * t32 - t86 + t9 * t94 + (-qJD(2) * t33 * t40 - t86) * t80) * pkin(2) + t55, t32 * t100 + (-t33 * t77 + t63) * t80 + t63 + (t10 * t94 + t44 * t91) * pkin(2) - t69 * t31 * t14, -t62 * t77 + ((t1 + (qJD(2) * t80 - qJD(3)) * t9) * t44 - t103 * t40) * pkin(2), (t79 + (t40 ^ 2 + t44 ^ 2) * qJDD(2)) * t52 + t56 * pkin(2), -t60, t97 * t66, t39 * t64, t60, t43 * t64, qJDD(4) * t71, (-g(3) * t43 + t59 * t39) * t71, (g(3) * t39 + t59 * t43) * t71, 0, 0;];
tau_reg = t5;
