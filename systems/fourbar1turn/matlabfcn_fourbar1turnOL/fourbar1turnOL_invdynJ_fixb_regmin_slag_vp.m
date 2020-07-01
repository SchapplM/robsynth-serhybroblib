% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnOL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:29
% EndTime: 2020-06-27 16:56:31
% DurationCPUTime: 0.31s
% Computational Cost: add. (220->67), mult. (578->136), div. (0->0), fcn. (438->10), ass. (0->56)
t37 = cos(qJ(2));
t70 = 0.2e1 * t37;
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t36 = cos(qJ(3));
t46 = t32 * t37 + t36 * t33;
t11 = t46 * qJD(1);
t12 = t32 * t33 - t36 * t37;
t25 = qJD(2) + qJD(3);
t7 = t25 * t46;
t4 = qJD(1) * t7 + qJDD(1) * t12;
t64 = qJD(1) * t33;
t52 = t32 * t64;
t63 = qJD(1) * t37;
t10 = -t36 * t63 + t52;
t69 = t11 * t10;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t68 = t31 * t35;
t26 = t31 ^ 2;
t67 = -t35 ^ 2 + t26;
t27 = t33 ^ 2;
t66 = -t37 ^ 2 + t27;
t62 = qJD(2) * t25;
t61 = qJD(2) * t33;
t60 = qJD(3) * t25;
t58 = t33 * qJDD(1);
t57 = t35 * qJDD(1);
t56 = t37 * qJDD(1);
t55 = qJD(1) * qJD(2);
t54 = qJD(1) * qJD(4);
t53 = pkin(2) * t63;
t51 = pkin(2) * qJD(2) * qJD(3);
t50 = t33 * t55;
t49 = -0.2e1 * pkin(1) * t54;
t34 = sin(qJ(1));
t38 = cos(qJ(1));
t48 = g(1) * t38 + g(2) * t34;
t47 = g(1) * t34 - g(2) * t38;
t41 = qJD(1) ^ 2;
t45 = pkin(1) * t41 + t48;
t44 = 0.2e1 * pkin(1) * qJDD(1) + t47;
t30 = qJ(2) + qJ(3);
t22 = sin(t30);
t23 = cos(t30);
t43 = g(3) * t23 - t11 * t53 - t48 * t22 + t32 * t51;
t42 = t32 * qJDD(2) * pkin(2) - g(3) * t22 + t10 * t53 - t48 * t23 + t36 * t51;
t3 = qJD(3) * t52 + (-t25 * t63 - t58) * t36 + (-t56 + t50) * t32;
t40 = qJD(2) ^ 2;
t39 = qJD(4) ^ 2;
t24 = qJDD(2) + qJDD(3);
t6 = t25 * t12;
t5 = -t10 ^ 2 + t11 ^ 2;
t2 = -t11 * t25 + t4;
t1 = -t10 * t25 + t3;
t8 = [qJDD(1), t47, t48, t27 * qJDD(1) + t50 * t70, 0.2e1 * t33 * t56 - 0.2e1 * t66 * t55, qJDD(2) * t33 + t40 * t37, qJDD(2) * t37 - t40 * t33, 0, t47 * t37, -t47 * t33, -t11 * t6 - t3 * t46, t6 * t10 - t11 * t7 + t3 * t12 - t4 * t46, -t24 * t46 + t6 * t25, t12 * t24 + t7 * t25, 0, -t47 * t23 + ((-qJD(1) * t12 - t10) * t61 + t4 * t70) * pkin(2), t47 * t22 + (-0.2e1 * t11 * t61 + (-qJD(1) * t6 + qJDD(1) * t46 - t3) * t37) * pkin(2), t26 * qJDD(1) + 0.2e1 * t54 * t68, 0.2e1 * t31 * t57 - 0.2e1 * t67 * t54, qJDD(4) * t31 + t39 * t35, qJDD(4) * t35 - t39 * t31, 0, t31 * t49 + t44 * t35, -t44 * t31 + t35 * t49; 0, 0, 0, -t33 * t41 * t37, t66 * t41, t58, t56, qJDD(2), -g(3) * t37 + t48 * t33, g(3) * t33 + t48 * t37, t69, t5, t1, t2, t24, (t10 * t64 + t32 * t60 + (-qJDD(2) - t24) * t36) * pkin(2) + t43, (t11 * t64 + t24 * t32 + t36 * t60) * pkin(2) + t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t5, t1, t2, t24, (-qJDD(2) * t36 - t32 * t62) * pkin(2) + t43, -t36 * pkin(2) * t62 + t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t68, t67 * t41, t31 * qJDD(1), t57, qJDD(4), -g(3) * t35 + t45 * t31, g(3) * t31 + t45 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t8;
