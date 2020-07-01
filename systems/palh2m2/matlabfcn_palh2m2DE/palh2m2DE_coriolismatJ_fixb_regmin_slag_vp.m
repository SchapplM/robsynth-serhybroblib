% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(4*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = palh2m2DE_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:48
% EndTime: 2020-06-30 17:56:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (97->46), mult. (224->79), div. (0->0), fcn. (155->6), ass. (0->54)
t26 = cos(qJ(3));
t57 = 0.2e1 * t26 ^ 2 - 0.1e1;
t27 = cos(qJ(2));
t56 = 0.2e1 * t27 ^ 2 - 0.1e1;
t22 = sin(qJ(4));
t12 = t27 * pkin(4) + pkin(1) + pkin(2);
t6 = t26 * pkin(5) + pkin(3) + t12;
t55 = t22 * t6;
t25 = cos(qJ(4));
t24 = sin(qJ(2));
t45 = t24 * qJD(2);
t17 = pkin(4) * t45;
t23 = sin(qJ(3));
t46 = t23 * qJD(3);
t9 = -pkin(5) * t46 - t17;
t54 = t25 * t9;
t53 = pkin(4) * qJD(2);
t52 = pkin(4) * qJD(3);
t51 = qJD(4) * t6;
t50 = qJD(1) * t23;
t49 = qJD(1) * t24;
t48 = qJD(1) * t25;
t47 = qJD(1) * t26;
t44 = t26 * qJD(3);
t43 = t27 * qJD(2);
t42 = pkin(1) * qJD(1);
t41 = pkin(1) * qJD(2);
t19 = qJD(1) + qJD(4);
t40 = pkin(4) * t19 * t24;
t39 = pkin(5) * t19 * t23;
t38 = pkin(5) * t50;
t15 = pkin(4) * t49;
t37 = t12 * t50;
t36 = t23 * t47;
t35 = t27 * t49;
t34 = t12 * t47;
t33 = t24 * t42;
t32 = t27 * t42;
t31 = t22 * t40;
t30 = t22 * t39;
t29 = t23 * t15;
t28 = t26 * t15;
t14 = t56 * qJD(1);
t13 = t57 * qJD(1);
t11 = t25 * t40;
t10 = t25 * t39;
t8 = -t27 * t23 + t26 * t24;
t7 = -t24 * t23 - t27 * t26;
t5 = t22 * t9;
t4 = pkin(4) * (t22 * t43 - t24 * t48);
t3 = pkin(4) * (t22 * t49 + t25 * t43);
t2 = pkin(5) * (t22 * t44 - t23 * t48);
t1 = pkin(5) * (t22 * t50 + t25 * t44);
t16 = [0, 0, 0, t24 * t43, t56 * qJD(2), 0, 0, 0, -t24 * t41, -t27 * t41, -t17, t23 * t44, t57 * qJD(3), 0, 0, 0, -t12 * t46 - t26 * t17, -t12 * t44 + t23 * t17, t9, 0, -t22 * t51 + t54, -t25 * t51 - t5; 0, 0, 0, t35, t14, t43, -t45, 0, -t33, -t32, -t15, 0, 0, 0, 0, 0, -t28, t29, -t15, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t13, t44, -t46, 0, -t37, -t34, -t38, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t55, -t25 * t6 * t19; 0, 0, 0, -t35, -t14, 0, 0, 0, t33, t32, t15, 0, 0, 0, 0, 0, t28, -t29, t15, 0, t11, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t52, t7 * t52, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t13, 0, 0, 0, t37, t34, t38, 0, t10, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t53, -t7 * t53, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t55 + t54, t6 * t48 - t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
