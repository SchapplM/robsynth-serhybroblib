% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = fourbar1turnOL_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:32
% EndTime: 2020-06-27 16:56:33
% DurationCPUTime: 0.27s
% Computational Cost: add. (98->33), mult. (260->57), div. (0->0), fcn. (284->6), ass. (0->43)
t16 = sin(qJ(3));
t17 = sin(qJ(2));
t19 = cos(qJ(3));
t20 = cos(qJ(2));
t10 = t16 * t20 + t19 * t17;
t9 = t16 * t17 - t19 * t20;
t44 = t10 * t9;
t43 = pkin(2) * qJD(2);
t42 = pkin(2) * qJD(3);
t3 = -t10 ^ 2 + t9 ^ 2;
t41 = t3 * qJD(1);
t4 = (t10 * t20 - t17 * t9) * pkin(2);
t40 = t4 * qJD(1);
t5 = (-t10 * t17 - t20 * t9) * pkin(2);
t39 = t5 * qJD(1);
t18 = cos(qJ(4));
t38 = qJD(1) * t18;
t37 = qJD(1) * t20;
t15 = sin(qJ(4));
t11 = -t15 ^ 2 + t18 ^ 2;
t36 = t11 * qJD(1);
t12 = -t17 ^ 2 + t20 ^ 2;
t35 = t12 * qJD(1);
t34 = t15 * qJD(4);
t33 = t18 * qJD(4);
t32 = t20 * qJD(2);
t31 = qJD(2) + qJD(3);
t30 = pkin(1) * t15 * qJD(1);
t29 = pkin(1) * t38;
t28 = pkin(2) * t37;
t27 = t20 * t42;
t26 = t16 * t43;
t25 = t19 * t43;
t24 = t15 * t38;
t23 = t17 * t37;
t22 = t9 * t28;
t21 = t10 * t28;
t14 = t19 * t42;
t13 = t16 * t42;
t8 = qJD(1) * t44;
t7 = t31 * t10;
t6 = t31 * t9;
t1 = [0, 0, 0, t17 * t32, t12 * qJD(2), 0, 0, 0, 0, 0, -t31 * t44, t31 * t3, 0, 0, 0, t4 * qJD(2) + t10 * t27, t5 * qJD(2) - t9 * t27, t15 * t33, t11 * qJD(4), 0, 0, 0, -pkin(1) * t34, -pkin(1) * t33; 0, 0, 0, t23, t35, t32, -t17 * qJD(2), 0, 0, 0, -t8, t41, t6, t7, 0, t40, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t41, t6, t7, 0, t21, -t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t36, t33, -t34, 0, -t30, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t23, -t35, 0, 0, 0, 0, 0, t8, -t41, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 + t26, t14 + t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t41, 0, 0, 0, -t21, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t36, 0, 0, 0, t30, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
