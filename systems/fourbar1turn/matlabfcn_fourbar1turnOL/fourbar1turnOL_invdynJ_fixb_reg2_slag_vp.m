% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnOL_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:10
% EndTime: 2020-04-12 19:41:11
% DurationCPUTime: 0.39s
% Computational Cost: add. (290->84), mult. (817->166), div. (0->0), fcn. (610->10), ass. (0->64)
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t36 = cos(qJ(3));
t37 = cos(qJ(2));
t14 = t32 * t37 + t36 * t33;
t12 = qJD(1) * t14;
t25 = qJD(2) + qJD(3);
t13 = t32 * t33 - t36 * t37;
t7 = t25 * t14;
t4 = qJD(1) * t7 + qJDD(1) * t13;
t2 = -t12 * t25 + t4;
t6 = t25 * t13;
t3 = qJD(1) * t6 - qJDD(1) * t14;
t68 = qJD(1) * t37;
t69 = qJD(1) * t33;
t11 = t32 * t69 - t36 * t68;
t77 = t11 * t12;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t75 = t31 * t35;
t26 = t31 ^ 2;
t28 = t35 ^ 2;
t74 = t26 - t28;
t27 = t33 ^ 2;
t29 = t37 ^ 2;
t73 = t27 - t29;
t72 = pkin(2) * qJD(2);
t67 = qJD(2) * t33;
t66 = qJD(3) * t25;
t65 = qJDD(1) * pkin(1);
t62 = t29 * qJDD(1);
t61 = t35 * qJDD(1);
t60 = t37 * qJDD(1);
t59 = qJD(1) * qJD(2);
t58 = qJD(1) * qJD(4);
t41 = qJD(1) ^ 2;
t57 = t41 * t75;
t56 = t33 * t41 * t37;
t55 = pkin(2) * t68;
t54 = t36 * t72;
t53 = qJD(3) * t72;
t52 = t33 * t59;
t51 = -0.2e1 * pkin(1) * t58;
t50 = t58 * t75;
t49 = t37 * t52;
t34 = sin(qJ(1));
t38 = cos(qJ(1));
t15 = g(1) * t38 + g(2) * t34;
t48 = g(1) * t34 - g(2) * t38;
t47 = pkin(1) * t41 + t15;
t46 = t48 + 0.2e1 * t65;
t45 = -g(3) * t37 + t15 * t33;
t30 = qJ(2) + qJ(3);
t22 = sin(t30);
t23 = cos(t30);
t44 = g(3) * t23 - t12 * t55 - t15 * t22 + t32 * t53;
t43 = t32 * qJDD(2) * pkin(2) - g(3) * t22 + t11 * t55 - t15 * t23 + t36 * t53;
t42 = pkin(2) ^ 2;
t40 = qJD(2) ^ 2;
t39 = qJD(4) ^ 2;
t24 = qJDD(2) + qJDD(3);
t5 = -t11 ^ 2 + t12 ^ 2;
t1 = -t11 * t25 + t3;
t8 = [0, 0, 0, 0, 0, qJDD(1), t48, t15, 0, 0, t27 * qJDD(1) + 0.2e1 * t49, 0.2e1 * t33 * t60 - 0.2e1 * t73 * t59, qJDD(2) * t33 + t40 * t37, -0.2e1 * t49 + t62, qJDD(2) * t37 - t40 * t33, 0, t48 * t37, -t48 * t33, -t15, 0, -t12 * t6 - t3 * t14, t6 * t11 - t12 * t7 + t3 * t13 - t14 * t4, -t14 * t24 + t6 * t25, t11 * t7 + t4 * t13, t13 * t24 + t7 * t25, 0, -t48 * t23 + ((-qJD(1) * t13 - t11) * t67 + 0.2e1 * t37 * t4) * pkin(2), t48 * t22 + 0.2e1 * (-t12 * t67 - t3 * t37) * pkin(2), ((-t13 * t32 - t14 * t36) * qJDD(2) + (-t32 * t7 + t36 * t6 + (-t13 * t36 + t14 * t32) * qJD(3)) * qJD(2)) * pkin(2) - t15, t42 * t62 + (t48 * pkin(2) - 0.2e1 * t42 * t52) * t37, t26 * qJDD(1) + 0.2e1 * t50, 0.2e1 * t31 * t61 - 0.2e1 * t74 * t58, qJDD(4) * t31 + t39 * t35, t28 * qJDD(1) - 0.2e1 * t50, qJDD(4) * t35 - t39 * t31, 0, t31 * t51 + t46 * t35, -t46 * t31 + t35 * t51, -t15, (t48 + t65) * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t73 * t41, t33 * qJDD(1), t56, t60, qJDD(2), t45, g(3) * t33 + t15 * t37, 0, 0, t77, t5, t1, -t77, t2, t24, (t11 * t69 + t32 * t66 + (-qJDD(2) - t24) * t36) * pkin(2) + t44, (t12 * t69 + t24 * t32 + t36 * t66) * pkin(2) + t43, -t11 * t54 + ((-qJD(3) * t11 + t3) * t36 - t2 * t32) * pkin(2), (t56 + (t32 ^ 2 + t36 ^ 2) * qJDD(2)) * t42 + t45 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t5, t1, -t77, t2, t24, (-qJD(2) * t25 * t32 - qJDD(2) * t36) * pkin(2) + t44, -t25 * t54 + t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t74 * t41, t31 * qJDD(1), t57, t61, qJDD(4), -g(3) * t35 + t47 * t31, g(3) * t31 + t47 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t8;
