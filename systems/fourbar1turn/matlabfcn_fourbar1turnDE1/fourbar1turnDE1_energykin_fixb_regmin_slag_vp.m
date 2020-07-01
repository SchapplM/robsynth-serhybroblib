% Calculate minimal parameter regressor of fixed base kinetic energy for
% fourbar1turnDE1
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
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnDE1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energykin_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_energykin_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:22
% EndTime: 2020-06-27 16:35:22
% DurationCPUTime: 0.40s
% Computational Cost: add. (3830->52), mult. (5468->141), div. (258->14), fcn. (1498->8), ass. (0->62)
t79 = pkin(4) ^ 2;
t78 = pkin(3) ^ 2;
t77 = -2 * pkin(2);
t38 = qJD(1) ^ 2;
t76 = t38 / 0.2e1;
t75 = (-pkin(3) - pkin(4));
t74 = (-pkin(3) + pkin(4));
t36 = sin(qJ(2));
t73 = pkin(1) * t36;
t37 = cos(qJ(2));
t72 = pkin(2) * t37;
t45 = pkin(2) ^ 2;
t46 = pkin(1) ^ 2;
t65 = -0.2e1 * pkin(1) * t72 + t46;
t31 = t45 + t65;
t64 = t78 - t79;
t27 = t31 + t64;
t33 = pkin(1) * t37 - pkin(2);
t25 = ((pkin(2) - t75) * (pkin(2) + t75)) + t65;
t26 = ((pkin(2) - t74) * (pkin(2) + t74)) + t65;
t47 = sqrt(-t25 * t26);
t21 = t27 * t73 - t33 * t47;
t71 = t21 * t36;
t29 = 0.1e1 / t31;
t35 = t36 ^ 2;
t70 = t29 * t35;
t30 = 0.1e1 / t31 ^ 2;
t69 = t30 / t79;
t28 = t31 - t64;
t68 = t36 * t28;
t67 = t36 * t47;
t66 = t37 * t47;
t63 = qJD(1) * t37;
t62 = t30 * t77;
t61 = qJD(1) * qJD(2);
t32 = pkin(1) - t72;
t20 = pkin(2) * t68 + t32 * t47;
t16 = t20 ^ 2;
t19 = -pkin(2) * t67 + t32 * t28;
t49 = t19 ^ 2;
t11 = (t16 + t49) * t69;
t40 = 0.1e1 / pkin(4);
t60 = t29 * t40 * t11 ^ (-0.1e1 / 0.2e1);
t17 = t21 ^ 2;
t43 = 0.1e1 / pkin(3);
t18 = -pkin(1) * t67 - t33 * t27;
t48 = t18 ^ 2;
t59 = t29 * t43 * ((t17 + t48) / t78 * t30) ^ (-0.1e1 / 0.2e1);
t58 = (-t25 - t26) * qJD(2) * pkin(2) * t73 / t47 * t29;
t55 = 0.1e1 / t11 * t38 * t69;
t54 = qJD(1) * t59;
t53 = t58 / 0.2e1;
t52 = pkin(1) * t38 * t60;
t15 = 0.1e1 / t49;
t5 = 0.2e1 * (-(t32 * t53 + (t45 * pkin(1) * t70 + ((t28 * t37 + t67) * t29 / 0.2e1 - t20 * t30 * t73) * pkin(2)) * qJD(2)) / t19 - (t36 * t53 + (-(-t66 + t68) * t29 / 0.2e1 + (t19 * t30 - t29 * t32) * t73) * qJD(2)) * pkin(2) * t20 * t15) * pkin(4) / (t16 * t15 + 0.1e1) * t31 * t40;
t51 = qJD(1) * t5 * t60;
t14 = 0.1e1 / t48;
t4 = qJD(2) + ((-t33 * t58 + (0.2e1 * t46 * pkin(2) * t70 + ((t27 * t37 + t67) * t29 + t62 * t71) * pkin(1)) * qJD(2)) / t18 - (-t36 * t58 + (-t29 * t66 + ((t33 * t77 + t27) * t29 + t18 * t62) * t36) * qJD(2)) * pkin(1) * t21 * t14) * pkin(3) / (t17 * t14 + 0.1e1) * t31 * t43;
t50 = qJD(2) * t4 * t59;
t7 = (-t18 * t37 + t71) * t54;
t6 = (-t18 * t36 - t21 * t37) * t54;
t1 = [t76, 0, 0, t35 * t76, t36 * t38 * t37, t36 * t61, t37 * t61, qJD(2) ^ 2 / 0.2e1, 0, 0, t6 ^ 2 / 0.2e1, t6 * t7, t6 * t4, t7 ^ 2 / 0.2e1, t7 * t4, t4 ^ 2 / 0.2e1, (-t18 * t50 + t7 * t63) * pkin(2), (t21 * t50 - t6 * t63) * pkin(2), t16 * t55 / 0.2e1, -t20 * t19 * t55, t20 * t51, -t19 * t51, t5 ^ 2 / 0.2e1, -t19 * t52, -t20 * t52;];
T_reg = t1;
