% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(2*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1turnTE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:19:50
% EndTime: 2020-04-12 19:19:50
% DurationCPUTime: 0.44s
% Computational Cost: add. (2988->51), mult. (4267->155), div. (178->13), fcn. (1148->4), ass. (0->66)
t82 = pkin(4) ^ 2;
t81 = pkin(3) ^ 2;
t43 = pkin(2) ^ 2;
t44 = pkin(1) ^ 2;
t34 = cos(qJ(2));
t72 = pkin(2) * t34;
t66 = -0.2e1 * pkin(1) * t72 + t44;
t27 = t43 + t66;
t25 = 0.1e1 / t27;
t80 = qJD(2) * t25;
t79 = -0.2e1 * pkin(2);
t65 = t81 - t82;
t23 = t27 + t65;
t29 = pkin(1) * t34 - pkin(2);
t33 = sin(qJ(2));
t75 = -pkin(3) - pkin(4);
t21 = (pkin(2) - t75) * (pkin(2) + t75) + t66;
t74 = -pkin(3) + pkin(4);
t22 = (pkin(2) - t74) * (pkin(2) + t74) + t66;
t45 = sqrt(-t21 * t22);
t68 = t33 * t45;
t14 = -pkin(1) * t68 - t29 * t23;
t78 = -t14 / 0.2e1;
t77 = -t25 / 0.2e1;
t36 = qJD(1) ^ 2;
t76 = t36 / 0.2e1;
t73 = pkin(1) * t33;
t17 = t23 * t73 - t29 * t45;
t71 = t17 * t33;
t31 = t33 ^ 2;
t70 = t25 * t31;
t24 = t27 - t65;
t69 = t33 * t24;
t67 = t34 * t45;
t64 = qJD(1) * t25;
t63 = qJD(1) * t34;
t26 = 0.1e1 / t27 ^ 2;
t62 = t26 * t79;
t61 = qJD(1) * qJD(2);
t60 = (-t21 - t22) * pkin(2) * t73 / t45 * t80;
t59 = t26 * t36 / t82;
t41 = 0.1e1 / pkin(3);
t58 = t41 * t64;
t57 = t41 * t80;
t56 = t34 ^ 2 * t76;
t28 = pkin(1) - t72;
t15 = -pkin(2) * t68 + t28 * t24;
t47 = t15 ^ 2;
t11 = 0.1e1 / t47;
t16 = pkin(2) * t69 + t28 * t45;
t12 = t16 ^ 2;
t38 = 0.1e1 / pkin(4);
t51 = t60 / 0.2e1;
t2 = 0.2e1 * (-(t28 * t51 + (t43 * pkin(1) * t70 + ((t24 * t34 + t68) * t25 / 0.2e1 - t16 * t26 * t73) * pkin(2)) * qJD(2)) / t15 - (t33 * t51 + ((-t67 + t69) * t77 + (t15 * t26 - t25 * t28) * t73) * qJD(2)) * pkin(2) * t16 * t11) * pkin(4) * t27 * t38 / (t12 * t11 + 0.1e1);
t53 = t2 * t38 * t64;
t46 = t14 ^ 2;
t10 = 0.1e1 / t46;
t13 = t17 ^ 2;
t1 = qJD(2) + ((-t29 * t60 + (0.2e1 * t44 * pkin(2) * t70 + ((t23 * t34 + t68) * t25 + t62 * t71) * pkin(1)) * qJD(2)) / t14 - (-t33 * t60 + (-t25 * t67 + ((t29 * t79 + t23) * t25 + t14 * t62) * t33) * qJD(2)) * pkin(1) * t17 * t10) * pkin(3) * t27 * t41 / (t13 * t10 + 0.1e1);
t52 = t1 * t57;
t50 = t59 / 0.8e1;
t48 = pkin(1) * t36 * t38 * t77;
t35 = qJD(2) ^ 2;
t5 = (t71 / 0.2e1 + t34 * t78) * t58;
t3 = (t14 * t33 + t17 * t34) * t58 / 0.2e1;
t4 = [0, 0, 0, 0, 0, t76, 0, 0, 0, 0, t31 * t76, t33 * t36 * t34, t33 * t61, t56, t34 * t61, t35 / 0.2e1, 0, 0, 0, 0, t3 ^ 2 / 0.2e1, -t3 * t5, -t3 * t1, t5 ^ 2 / 0.2e1, t5 * t1, t1 ^ 2 / 0.2e1, (t5 * t63 + t52 * t78) * pkin(2), (t3 * t63 + t17 * t52 / 0.2e1) * pkin(2), (-t17 * t5 / 0.2e1 + t3 * t78) * pkin(2) * t57, (t56 + (t13 / 0.8e1 + t46 / 0.8e1) / t81 * t35 * t26) * t43, t12 * t50, -t16 * t15 * t59 / 0.4e1, t16 * t53 / 0.2e1, t47 * t50, -t15 * t53 / 0.2e1, t2 ^ 2 / 0.2e1, t15 * t48, t16 * t48, 0, t44 * t76;];
T_reg = t4;
