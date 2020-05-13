% Calculate inertial parameters regressor of fixed base kinetic energy for
% fourbar1turnDE2
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
% T_reg [1x(2*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnDE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:34:54
% EndTime: 2020-04-12 19:34:55
% DurationCPUTime: 0.53s
% Computational Cost: add. (4411->55), mult. (6340->156), div. (315->15), fcn. (1732->8), ass. (0->68)
t85 = pkin(4) ^ 2;
t84 = pkin(3) ^ 2;
t83 = -2 * pkin(2);
t42 = qJD(1) ^ 2;
t82 = t42 / 0.2e1;
t81 = (-pkin(3) - pkin(4));
t80 = (-pkin(3) + pkin(4));
t39 = sin(qJ(2));
t79 = pkin(1) * t39;
t40 = cos(qJ(2));
t78 = pkin(2) * t40;
t49 = pkin(2) ^ 2;
t50 = pkin(1) ^ 2;
t70 = -0.2e1 * pkin(1) * t78 + t50;
t33 = t49 + t70;
t69 = t84 - t85;
t29 = t33 + t69;
t35 = pkin(1) * t40 - pkin(2);
t27 = ((pkin(2) - t81) * (pkin(2) + t81)) + t70;
t28 = ((pkin(2) - t80) * (pkin(2) + t80)) + t70;
t51 = sqrt(-t27 * t28);
t23 = t29 * t79 - t35 * t51;
t77 = t23 * t39;
t31 = 0.1e1 / t33;
t37 = t39 ^ 2;
t76 = t31 * t37;
t32 = 0.1e1 / t33 ^ 2;
t75 = t32 / t85;
t74 = t32 / t84;
t30 = t33 - t69;
t73 = t39 * t30;
t72 = t39 * t51;
t71 = t40 * t51;
t68 = qJD(1) * t40;
t67 = t32 * t83;
t66 = qJD(1) * qJD(2);
t34 = pkin(1) - t78;
t21 = -pkin(2) * t72 + t30 * t34;
t16 = t21 ^ 2;
t22 = pkin(2) * t73 + t34 * t51;
t18 = t22 ^ 2;
t10 = (t16 + t18) * t75;
t44 = 0.1e1 / pkin(4);
t65 = t31 * t44 * t10 ^ (-0.1e1 / 0.2e1);
t20 = -pkin(1) * t72 - t29 * t35;
t14 = t20 ^ 2;
t19 = t23 ^ 2;
t11 = (t14 + t19) * t74;
t47 = 0.1e1 / pkin(3);
t64 = t31 * t47 * t11 ^ (-0.1e1 / 0.2e1);
t63 = (-t27 - t28) * qJD(2) * pkin(2) * t79 / t51 * t31;
t62 = t40 ^ 2 * t82;
t59 = t42 / t10 * t75;
t58 = qJD(1) * t64;
t57 = qJD(2) * t64;
t56 = t63 / 0.2e1;
t55 = pkin(1) * t42 * t65;
t17 = 0.1e1 / t21 ^ 2;
t2 = 0.2e1 * (-(t34 * t56 + (t49 * pkin(1) * t76 + ((t30 * t40 + t72) * t31 / 0.2e1 - t22 * t32 * t79) * pkin(2)) * qJD(2)) / t21 - (t39 * t56 + (-(-t71 + t73) * t31 / 0.2e1 + (t21 * t32 - t31 * t34) * t79) * qJD(2)) * pkin(2) * t22 * t17) * pkin(4) / (t17 * t18 + 0.1e1) * t33 * t44;
t54 = qJD(1) * t2 * t65;
t15 = 0.1e1 / t20 ^ 2;
t1 = qJD(2) + ((-t35 * t63 + (0.2e1 * t50 * pkin(2) * t76 + ((t29 * t40 + t72) * t31 + t67 * t77) * pkin(1)) * qJD(2)) / t20 - (-t39 * t63 + (-t31 * t71 + ((t35 * t83 + t29) * t31 + t20 * t67) * t39) * qJD(2)) * pkin(1) * t23 * t15) * pkin(3) / (t15 * t19 + 0.1e1) * t33 * t47;
t53 = t1 * t57;
t52 = t59 / 0.2e1;
t41 = qJD(2) ^ 2;
t5 = (-t20 * t40 + t77) * t58;
t3 = (-t20 * t39 - t23 * t40) * t58;
t4 = [0, 0, 0, 0, 0, t82, 0, 0, 0, 0, t37 * t82, t39 * t42 * t40, t39 * t66, t62, t40 * t66, t41 / 0.2e1, 0, 0, 0, 0, t3 ^ 2 / 0.2e1, t3 * t5, t3 * t1, t5 ^ 2 / 0.2e1, t5 * t1, t1 ^ 2 / 0.2e1, (-t20 * t53 + t5 * t68) * pkin(2), (t23 * t53 - t3 * t68) * pkin(2), (t20 * t3 - t23 * t5) * pkin(2) * t57, (t62 + (t19 / 0.2e1 + t14 / 0.2e1) / t11 * t41 * t74) * t49, t18 * t52, -t22 * t21 * t59, t22 * t54, t16 * t52, -t21 * t54, t2 ^ 2 / 0.2e1, -t21 * t55, -t22 * t55, 0, t50 * t82;];
T_reg = t4;
