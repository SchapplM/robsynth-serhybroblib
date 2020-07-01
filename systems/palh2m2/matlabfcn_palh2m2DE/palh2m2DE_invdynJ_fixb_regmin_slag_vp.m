% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tau_reg [4x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh2m2DE_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.31s
% Computational Cost: add. (261->78), mult. (552->142), div. (0->0), fcn. (339->8), ass. (0->65)
t34 = sin(qJ(2));
t63 = qJD(1) * qJD(2);
t86 = 0.2e1 * t34 * t63;
t35 = sin(qJ(1));
t39 = cos(qJ(1));
t16 = g(1) * t39 + g(2) * t35;
t33 = sin(qJ(3));
t37 = cos(qJ(3));
t40 = qJD(3) ^ 2;
t51 = qJDD(3) * t33 + t37 * t40;
t38 = cos(qJ(2));
t41 = qJD(2) ^ 2;
t52 = qJDD(2) * t34 + t38 * t41;
t85 = pkin(4) * t52 + pkin(5) * t51 + t16;
t15 = g(1) * t35 - g(2) * t39;
t27 = qJD(1) + qJD(4);
t84 = t27 ^ 2;
t83 = pkin(4) * t86 - t15;
t82 = qJD(2) - qJD(3);
t66 = t38 * pkin(4) + pkin(1);
t17 = pkin(2) + t66;
t81 = -0.2e1 * t17;
t80 = 2 * qJDD(1);
t79 = pkin(4) * t34;
t78 = pkin(5) * t33;
t36 = cos(qJ(4));
t76 = t85 * t36;
t42 = qJD(1) ^ 2;
t75 = t34 * t42;
t14 = -qJD(2) * t79 - qJD(3) * t78;
t74 = t36 * t14;
t73 = t37 * t42;
t28 = t33 ^ 2;
t72 = -t37 ^ 2 + t28;
t29 = t34 ^ 2;
t71 = -t38 ^ 2 + t29;
t56 = t37 * pkin(5) + t17;
t9 = pkin(3) + t56;
t70 = qJD(1) * t9;
t69 = t14 * qJD(1);
t68 = t9 * qJDD(1);
t67 = -qJD(4) + t27;
t65 = t37 * qJDD(1);
t64 = t38 * qJDD(1);
t62 = qJD(1) * qJD(3);
t32 = sin(qJ(4));
t61 = t32 * t70;
t59 = t62 * t81;
t58 = -qJD(4) * t74 + (t68 + t69) * t36 + t85 * t32;
t57 = -0.2e1 * pkin(1) * t63;
t55 = qJD(4) * t9 * (-qJD(1) - t27);
t54 = t33 * t38 - t34 * t37;
t53 = t33 * t34 + t37 * t38;
t50 = (-qJD(1) - t67) * t14 - t15;
t48 = t17 * t42 + t16;
t47 = pkin(1) * t42 + t16;
t46 = pkin(1) * t80 + t15;
t26 = qJDD(1) + qJDD(4);
t45 = t26 * t32 + t36 * t84;
t44 = t26 * t36 - t32 * t84;
t43 = qJDD(1) * t81 + t83;
t20 = pkin(4) * t75;
t2 = t82 * t54;
t1 = t82 * t53;
t3 = [qJDD(1), t15, t16, qJDD(1) * t29 + t38 * t86, 0.2e1 * t34 * t64 - 0.2e1 * t63 * t71, t52, qJDD(2) * t38 - t34 * t41, 0, t34 * t57 + t38 * t46, -t34 * t46 + t38 * t57, t66 * t80 - t83, 0.2e1 * t33 * t37 * t62 + qJDD(1) * t28, 0.2e1 * t33 * t65 - 0.2e1 * t62 * t72, t51, qJDD(3) * t37 - t33 * t40, 0, t33 * t59 - t37 * t43, t33 * t43 + t37 * t59, t56 * t80 + t15 + 0.2e1 * t69, t26, t32 * t55 + (t27 * t14 + t26 * t9 + t15) * t36 + t58, t36 * t55 + ((-qJDD(1) - t26) * t9 + t50) * t32 + t76; 0, 0, 0, -t38 * t75, t71 * t42, t34 * qJDD(1), t64, qJDD(2), -t38 * g(3) + t34 * t47, g(3) * t34 + t38 * t47, t20, 0, 0, 0, 0, 0, (t34 * t73 + qJDD(3) * t53 + (-qJD(2) * t54 + t2) * qJD(3)) * pkin(4), (-t33 * t75 - qJDD(3) * t54 + (-qJD(2) * t53 + t1) * qJD(3)) * pkin(4), t20, 0, t45 * t79, t44 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t73, t72 * t42, t33 * qJDD(1), t65, qJDD(3), -t37 * g(3) + t48 * t33 + (qJDD(2) * t53 + (qJD(3) * t54 + t2) * qJD(2)) * pkin(4), g(3) * t33 + t48 * t37 + (-qJDD(2) * t54 + (qJD(3) * t53 + t1) * qJD(2)) * pkin(4), t42 * t78, 0, t45 * t78, t44 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -qJD(4) * t61 - t27 * (-t61 - t74) + t15 * t36 + t58, t67 * t36 * t70 + (t50 - t68) * t32 + t76;];
tau_reg = t3;
