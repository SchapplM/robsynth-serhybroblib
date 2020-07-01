% Calculate minimal parameter regressor of fixed base kinetic energy for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnTE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:29
% EndTime: 2020-06-27 16:22:30
% DurationCPUTime: 0.37s
% Computational Cost: add. (2794->46), mult. (3959->141), div. (160->12), fcn. (1070->4), ass. (0->59)
t72 = pkin(4) ^ 2;
t39 = pkin(2) ^ 2;
t40 = pkin(1) ^ 2;
t32 = cos(qJ(2));
t62 = pkin(2) * t32;
t56 = -0.2e1 * pkin(1) * t62 + t40;
t26 = t39 + t56;
t24 = 0.1e1 / t26;
t71 = qJD(2) * t24;
t70 = -0.2e1 * pkin(2);
t55 = pkin(3) ^ 2 - t72;
t22 = t26 + t55;
t28 = pkin(1) * t32 - pkin(2);
t31 = sin(qJ(2));
t65 = -pkin(3) - pkin(4);
t20 = (pkin(2) - t65) * (pkin(2) + t65) + t56;
t64 = -pkin(3) + pkin(4);
t21 = (pkin(2) - t64) * (pkin(2) + t64) + t56;
t41 = sqrt(-t20 * t21);
t58 = t31 * t41;
t13 = -pkin(1) * t58 - t28 * t22;
t69 = -t13 / 0.2e1;
t68 = -t24 / 0.2e1;
t67 = -t32 / 0.2e1;
t33 = qJD(1) ^ 2;
t66 = t33 / 0.2e1;
t63 = pkin(1) * t31;
t16 = t22 * t63 - t28 * t41;
t61 = t16 * t31;
t30 = t31 ^ 2;
t60 = t24 * t30;
t23 = t26 - t55;
t59 = t31 * t23;
t57 = t32 * t41;
t54 = qJD(1) * t24;
t53 = qJD(1) * t32;
t25 = 0.1e1 / t26 ^ 2;
t52 = t25 * t70;
t51 = qJD(1) * qJD(2);
t50 = (-t20 - t21) * pkin(2) * t63 / t41 * t71;
t49 = t25 * t33 / t72;
t38 = 0.1e1 / pkin(3);
t48 = t38 * t54;
t35 = 0.1e1 / pkin(4);
t27 = pkin(1) - t62;
t14 = -pkin(2) * t58 + t27 * t23;
t11 = 0.1e1 / t14 ^ 2;
t15 = pkin(2) * t59 + t27 * t41;
t12 = t15 ^ 2;
t43 = t50 / 0.2e1;
t5 = 0.2e1 * (-(t27 * t43 + (t39 * pkin(1) * t60 + ((t23 * t32 + t58) * t24 / 0.2e1 - t15 * t25 * t63) * pkin(2)) * qJD(2)) / t14 - (t31 * t43 + ((-t57 + t59) * t68 + (t14 * t25 - t24 * t27) * t63) * qJD(2)) * pkin(2) * t15 * t11) * pkin(4) * t26 * t35 / (t12 * t11 + 0.1e1);
t45 = t35 * t5 * t54;
t10 = 0.1e1 / t13 ^ 2;
t4 = qJD(2) + ((-t28 * t50 + (0.2e1 * t40 * pkin(2) * t60 + ((t22 * t32 + t58) * t24 + t52 * t61) * pkin(1)) * qJD(2)) / t13 - (-t31 * t50 + (-t24 * t57 + ((t28 * t70 + t22) * t24 + t13 * t52) * t31) * qJD(2)) * pkin(1) * t16 * t10) * pkin(3) * t26 * t38 / (t16 ^ 2 * t10 + 0.1e1);
t44 = t38 * t4 * t71;
t42 = pkin(1) * t33 * t35 * t68;
t7 = (t61 / 0.2e1 + t13 * t67) * t48;
t6 = (t16 * t67 + t31 * t69) * t48;
t1 = [t66, 0, 0, t30 * t66, t31 * t33 * t32, t31 * t51, t32 * t51, qJD(2) ^ 2 / 0.2e1, 0, 0, t6 ^ 2 / 0.2e1, t6 * t7, t6 * t4, t7 ^ 2 / 0.2e1, t7 * t4, t4 ^ 2 / 0.2e1, (t44 * t69 + t7 * t53) * pkin(2), (-t6 * t53 + t16 * t44 / 0.2e1) * pkin(2), t12 * t49 / 0.8e1, -t15 * t14 * t49 / 0.4e1, t15 * t45 / 0.2e1, -t14 * t45 / 0.2e1, t5 ^ 2 / 0.2e1, t14 * t42, t15 * t42;];
T_reg = t1;
