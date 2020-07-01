% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1turnOL_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:29
% EndTime: 2020-06-27 16:56:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (124->34), mult. (422->92), div. (0->0), fcn. (294->6), ass. (0->43)
t19 = qJD(2) + qJD(3);
t25 = sin(qJ(3));
t26 = sin(qJ(2));
t28 = cos(qJ(3));
t29 = cos(qJ(2));
t33 = t25 * t29 + t28 * t26;
t11 = t33 * qJD(1);
t48 = pkin(2) * qJD(2);
t36 = qJD(3) * t48;
t46 = qJD(1) * t29;
t41 = pkin(2) * t46;
t54 = -t11 * t41 + t25 * t36;
t38 = t28 * t46;
t47 = qJD(1) * t26;
t39 = t25 * t47;
t10 = -t38 + t39;
t53 = t10 * t41 + t28 * t36;
t52 = t11 * t10;
t27 = cos(qJ(4));
t32 = qJD(1) ^ 2;
t51 = t27 * t32;
t24 = sin(qJ(4));
t50 = t24 ^ 2 - t27 ^ 2;
t49 = t26 ^ 2 - t29 ^ 2;
t45 = qJD(2) * t26;
t44 = qJD(3) * t19;
t43 = qJD(1) * qJD(2);
t42 = qJD(1) * qJD(4);
t40 = t19 * t48;
t37 = 0.2e1 * t42;
t35 = t26 * t43;
t34 = -0.2e1 * pkin(1) * t42;
t12 = t25 * t26 - t28 * t29;
t4 = qJD(3) * t39 - t19 * t38 + t25 * t35;
t7 = t19 * t33;
t5 = qJD(1) * t7;
t31 = qJD(2) ^ 2;
t30 = qJD(4) ^ 2;
t6 = t19 * t12;
t3 = -t10 ^ 2 + t11 ^ 2;
t2 = -t11 * t19 + t5;
t1 = -t10 * t19 + t4;
t8 = [0, 0, 0, 0.2e1 * t29 * t35, -0.2e1 * t49 * t43, t31 * t29, -t31 * t26, 0, 0, 0, -t11 * t6 - t33 * t4, t6 * t10 - t11 * t7 + t4 * t12 - t33 * t5, t6 * t19, t7 * t19, 0, (-t10 * t45 + t29 * t5 + (-t12 * t45 + t29 * t7) * qJD(1)) * pkin(2), (-t11 * t45 - t29 * t4 + (-t29 * t6 - t33 * t45) * qJD(1)) * pkin(2), t24 * t27 * t37, -t50 * t37, t30 * t27, -t30 * t24, 0, t24 * t34, t27 * t34; 0, 0, 0, -t26 * t32 * t29, t49 * t32, 0, 0, 0, 0, 0, t52, t3, t1, t2, 0, (t10 * t47 + t25 * t44) * pkin(2) + t54, (t11 * t47 + t28 * t44) * pkin(2) + t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t3, t1, t2, 0, -t25 * t40 + t54, -t28 * t40 + t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t51, t50 * t32, 0, 0, 0, t32 * pkin(1) * t24, pkin(1) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg = t8;
