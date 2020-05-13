% Calculate inertial parameters regressor of fixed base kinetic energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% T_reg [1x(10*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = palh3m2OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energykin_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_energykin_fixb_reg2_slag_vp: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energykin_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:40:29
% EndTime: 2020-05-07 04:40:29
% DurationCPUTime: 0.28s
% Computational Cost: add. (453->67), mult. (1160->200), div. (0->0), fcn. (886->16), ass. (0->64)
t58 = qJD(1) ^ 2;
t73 = t58 / 0.2e1;
t41 = qJD(2) + qJD(7);
t72 = pkin(3) * t41;
t71 = cos(qJ(5));
t70 = sin(qJ(8));
t53 = cos(qJ(6));
t69 = t53 * t58;
t56 = cos(qJ(2));
t68 = t56 * t58;
t57 = qJD(2) ^ 2;
t67 = t57 * pkin(1) ^ 2;
t66 = pkin(1) * qJD(2);
t65 = qJD(1) * qJD(2);
t64 = qJD(1) * qJD(6);
t42 = qJD(2) + qJD(3);
t49 = sin(qJ(3));
t63 = t49 * t66;
t45 = sin(qJ(7));
t62 = t45 * t66;
t52 = cos(qJ(7));
t61 = t52 * t66;
t55 = cos(qJ(3));
t60 = t55 * t66;
t50 = sin(qJ(2));
t24 = (t49 * t50 - t55 * t56) * qJD(1);
t26 = (-t49 * t56 - t50 * t55) * qJD(1);
t48 = sin(qJ(4));
t54 = cos(qJ(4));
t13 = -t54 * t24 + t48 * t26;
t35 = (-pkin(1) * t56 - pkin(12)) * qJD(1);
t34 = t42 * pkin(4) - t60;
t20 = t48 * t34 - t54 * t63;
t19 = t54 * t34 + t48 * t63;
t18 = -t24 * pkin(4) + t35;
t51 = cos(qJ(8));
t47 = sin(qJ(5));
t46 = sin(qJ(6));
t44 = cos(pkin(15));
t43 = sin(pkin(15));
t40 = qJD(4) + t42;
t39 = qJD(8) + t41;
t33 = t35 ^ 2 / 0.2e1;
t32 = -t43 * t70 - t44 * t51;
t31 = t43 * t51 - t44 * t70;
t30 = t44 * t72 + t61;
t29 = t43 * t72 + t62;
t25 = (t45 * t56 + t50 * t52) * qJD(1);
t22 = (t45 * t50 - t52 * t56) * qJD(1);
t17 = t40 * pkin(10) + t20;
t16 = -t40 * pkin(8) - t19;
t15 = t48 * t24 + t54 * t26;
t12 = qJD(5) + t13;
t11 = t35 + (t22 * t44 - t25 * t43) * pkin(3);
t10 = t71 * t15 + t47 * t40;
t8 = t47 * t15 - t71 * t40;
t7 = t32 * t29 + t31 * t30;
t6 = -t31 * t29 + t32 * t30;
t5 = -t31 * t22 + t32 * t25;
t4 = -t32 * t22 - t31 * t25;
t3 = t13 * pkin(8) - t15 * pkin(10) + t18;
t2 = t71 * t17 + t47 * t3;
t1 = -t47 * t17 + t71 * t3;
t9 = [0, 0, 0, 0, 0, t73, 0, 0, 0, 0, t50 ^ 2 * t73, t50 * t68, t50 * t65, t56 ^ 2 * t73, t56 * t65, t57 / 0.2e1, pkin(12) * t68, -t58 * pkin(12) * t50, 0, pkin(12) ^ 2 * t73, t26 ^ 2 / 0.2e1, t26 * t24, t26 * t42, t24 ^ 2 / 0.2e1, t24 * t42, t42 ^ 2 / 0.2e1, -t35 * t24 - t42 * t60, t35 * t26 + t42 * t63, (-t24 * t49 + t26 * t55) * t66, t33 + (t49 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1) * t67, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t40, t13 ^ 2 / 0.2e1, -t13 * t40, t40 ^ 2 / 0.2e1, t18 * t13 + t19 * t40, t18 * t15 - t20 * t40, -t20 * t13 - t19 * t15, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t12, t8 ^ 2 / 0.2e1, -t8 * t12, t12 ^ 2 / 0.2e1, t1 * t12 + t16 * t8, t16 * t10 - t2 * t12, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t46 ^ 2 * t73, t46 * t69, t46 * t64, t53 ^ 2 * t73, t53 * t64, qJD(6) ^ 2 / 0.2e1, -pkin(6) * t69, t58 * pkin(6) * t46, 0, pkin(6) ^ 2 * t73, t25 ^ 2 / 0.2e1, -t25 * t22, t25 * t41, t22 ^ 2 / 0.2e1, -t22 * t41, t41 ^ 2 / 0.2e1, t35 * t22 + t41 * t61, t35 * t25 - t41 * t62, (-t22 * t45 - t25 * t52) * t66, t33 + (t45 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * t67, t5 ^ 2 / 0.2e1, t5 * t4, t5 * t39, t4 ^ 2 / 0.2e1, t4 * t39, t39 ^ 2 / 0.2e1, -t11 * t4 + t6 * t39, t11 * t5 - t7 * t39, t7 * t4 - t6 * t5, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t9;
