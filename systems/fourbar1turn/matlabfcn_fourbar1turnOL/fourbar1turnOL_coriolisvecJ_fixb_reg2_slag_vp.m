% Calculate inertial parameters regressor of coriolis joint torque vector for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1turnOL_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:05
% EndTime: 2020-04-12 19:41:07
% DurationCPUTime: 0.26s
% Computational Cost: add. (168->40), mult. (597->108), div. (0->0), fcn. (422->6), ass. (0->47)
t42 = (qJD(1) * qJD(4));
t57 = -2 * t42;
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t27 = cos(qJ(3));
t28 = cos(qJ(2));
t14 = t24 * t28 + t27 * t25;
t12 = t14 * qJD(1);
t18 = qJD(2) + qJD(3);
t7 = t18 * t14;
t5 = qJD(1) * t7;
t2 = -t12 * t18 + t5;
t46 = qJD(1) * t28;
t47 = qJD(1) * t25;
t11 = t24 * t47 - t27 * t46;
t48 = pkin(2) * qJD(2);
t37 = qJD(3) * t48;
t39 = pkin(2) * t46;
t56 = t11 * t39 + t27 * t37;
t55 = t11 * t12;
t53 = t25 * t28;
t26 = cos(qJ(4));
t31 = qJD(1) ^ 2;
t52 = t26 * t31;
t51 = -t12 * t39 + t24 * t37;
t23 = sin(qJ(4));
t50 = t23 ^ 2 - t26 ^ 2;
t49 = t25 ^ 2 - t28 ^ 2;
t45 = qJD(2) * t25;
t44 = qJD(3) * t18;
t43 = qJD(1) * qJD(2);
t41 = t23 * t52;
t40 = t31 * t53;
t38 = t27 * t48;
t36 = pkin(1) * t57;
t35 = t23 * t26 * t42;
t34 = t43 * t53;
t33 = -0.2e1 * t34;
t13 = t24 * t25 - t27 * t28;
t6 = t18 * t13;
t4 = qJD(1) * t6;
t32 = pkin(2) ^ 2;
t30 = qJD(2) ^ 2;
t29 = qJD(4) ^ 2;
t3 = -t11 ^ 2 + t12 ^ 2;
t1 = -t11 * t18 + t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, -0.2e1 * t49 * t43, t30 * t28, t33, -t30 * t25, 0, 0, 0, 0, 0, -t12 * t6 - t4 * t14, t6 * t11 - t12 * t7 + t4 * t13 - t14 * t5, t6 * t18, t11 * t7 + t5 * t13, t7 * t18, 0, (-t11 * t45 + t28 * t5 + (-t13 * t45 + t28 * t7) * qJD(1)) * pkin(2), (-t12 * t45 - t28 * t4 + (-t14 * t45 - t28 * t6) * qJD(1)) * pkin(2), (-t24 * t7 + t27 * t6 + (-t13 * t27 + t14 * t24) * qJD(3)) * t48, t32 * t33, 0.2e1 * t35, t50 * t57, t29 * t26, -0.2e1 * t35, -t29 * t23, 0, t23 * t36, t26 * t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t49 * t31, 0, t40, 0, 0, 0, 0, 0, 0, t55, t3, t1, -t55, t2, 0, (t11 * t47 + t24 * t44) * pkin(2) + t51, (t12 * t47 + t27 * t44) * pkin(2) + t56, -t11 * t38 + ((-qJD(3) * t11 + t4) * t27 - t2 * t24) * pkin(2), t32 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t3, t1, -t55, t2, 0, -t24 * t18 * t48 + t51, -t18 * t38 + t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t50 * t31, 0, t41, 0, 0, t31 * pkin(1) * t23, pkin(1) * t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t8;
