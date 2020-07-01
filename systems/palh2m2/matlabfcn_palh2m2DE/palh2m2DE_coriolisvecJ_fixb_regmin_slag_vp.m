% Calculate minimal parameter regressor of coriolis joint torque vector for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% tauc_reg [4x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = palh2m2DE_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.21s
% Computational Cost: add. (141->44), mult. (377->105), div. (0->0), fcn. (198->6), ass. (0->53)
t17 = (qJD(1) + qJD(4));
t64 = t17 ^ 2;
t41 = (qJD(1) * qJD(2));
t63 = -2 * t41;
t62 = qJD(2) - qJD(3);
t61 = 2 * qJD(1);
t24 = sin(qJ(2));
t60 = pkin(4) * t24;
t23 = sin(qJ(3));
t59 = pkin(5) * t23;
t25 = cos(qJ(4));
t27 = cos(qJ(2));
t12 = t27 * pkin(4) + pkin(1) + pkin(2);
t26 = cos(qJ(3));
t7 = pkin(5) * t26 + pkin(3) + t12;
t58 = t25 * t7;
t48 = pkin(4) * qJD(2);
t40 = t24 * t48;
t46 = qJD(3) * t23;
t10 = -pkin(5) * t46 - t40;
t22 = sin(qJ(4));
t57 = t22 * t10;
t30 = qJD(1) ^ 2;
t56 = t23 * t30;
t55 = t24 * t30;
t54 = t25 * t10;
t53 = t26 * t30;
t28 = qJD(3) ^ 2;
t52 = t28 * t26;
t29 = qJD(2) ^ 2;
t51 = t29 * t27;
t50 = t23 ^ 2 - t26 ^ 2;
t49 = t24 ^ 2 - t27 ^ 2;
t47 = pkin(1) * t30;
t45 = t10 * qJD(1);
t44 = -qJD(1) - t17;
t43 = qJD(4) - t17;
t11 = pkin(4) * t51 + pkin(5) * t52;
t42 = -qJD(4) * t54 + t22 * t11 + t25 * t45;
t39 = t22 * t7 * qJD(1);
t38 = qJD(3) * t61;
t37 = t24 * t41;
t36 = t44 * t22;
t35 = pkin(1) * t63;
t34 = t64 * t22;
t33 = t64 * t25;
t32 = t23 * t27 - t24 * t26;
t31 = t23 * t24 + t26 * t27;
t13 = pkin(4) * t55;
t6 = t25 * t11;
t2 = t62 * t32;
t1 = t62 * t31;
t3 = [0, 0, 0, 0.2e1 * t27 * t37, t49 * t63, t51, -t29 * t24, 0, t24 * t35, t27 * t35, -0.2e1 * pkin(4) * t37, t23 * t26 * t38, -t50 * t38, t52, -t28 * t23, 0, (-t12 * t46 - t26 * t40) * t61, (-qJD(3) * t12 * t26 + t23 * t40) * t61, 0.2e1 * t45, 0, qJD(4) * t36 * t7 + t17 * t54 + t42, t6 + t10 * t36 + (t44 * t58 + t57) * qJD(4); 0, 0, 0, -t27 * t55, t49 * t30, 0, 0, 0, t24 * t47, t27 * t47, t13, 0, 0, 0, 0, 0, (t24 * t53 + (-qJD(2) * t32 + t2) * qJD(3)) * pkin(4), (-t23 * t55 + (-qJD(2) * t31 + t1) * qJD(3)) * pkin(4), t13, 0, t33 * t60, -t34 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23 * t53, t50 * t30, 0, 0, 0, t12 * t56 + (qJD(3) * t32 + t2) * t48, t12 * t53 + (qJD(3) * t31 + t1) * t48, pkin(5) * t56, 0, t33 * t59, -t34 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t39 - t17 * (-t39 - t54) + t42, t6 + t43 * t57 + (-t43 * t58 - t57) * qJD(1);];
tauc_reg = t3;
