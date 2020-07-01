% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = palh2m2DE_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (49->20), mult. (159->60), div. (0->0), fcn. (97->6), ass. (0->25)
t27 = pkin(4) * (qJD(2) - qJD(3));
t11 = sin(qJ(4));
t26 = qJD(4) * t11;
t14 = cos(qJ(4));
t25 = qJD(4) * t14;
t12 = sin(qJ(3));
t24 = t12 * qJD(3);
t13 = sin(qJ(2));
t23 = t13 * qJD(2);
t15 = cos(qJ(3));
t22 = t15 * qJD(3);
t16 = cos(qJ(2));
t21 = t16 * qJD(2);
t20 = pkin(4) * t23;
t19 = -0.2e1 * pkin(1) * qJD(2);
t9 = t16 * pkin(4) + pkin(1) + pkin(2);
t8 = -pkin(5) * t24 - t20;
t7 = pkin(5) * t15 + pkin(3) + t9;
t6 = (t11 * t21 + t13 * t25) * pkin(4);
t5 = (t11 * t22 + t12 * t25) * pkin(5);
t4 = (-t12 * t26 + t14 * t22) * pkin(5);
t3 = (-t13 * t26 + t14 * t21) * pkin(4);
t2 = -t11 * t8 - t25 * t7;
t1 = t14 * t8 - t26 * t7;
t10 = [0, 0, 0, 0.2e1 * t13 * t21, 0.2e1 * (-t13 ^ 2 + t16 ^ 2) * qJD(2), 0, 0, 0, t13 * t19, t16 * t19, -0.2e1 * t20, 0.2e1 * t12 * t22, 0.2e1 * (-t12 ^ 2 + t15 ^ 2) * qJD(3), 0, 0, 0, -0.2e1 * t15 * t20 - 0.2e1 * t24 * t9, 0.2e1 * t12 * t20 - 0.2e1 * t22 * t9, 0.2e1 * t8, 0, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, t21, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t24, 0, 0, 0, 0, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t12 * t16 - t13 * t15) * t27, (t12 * t13 + t15 * t16) * t27, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
