% Calculate inertial parameters regressor of fixed base kinetic energy for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% T_reg [1x(1*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1DE1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energykin_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_energykin_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:16
% EndTime: 2020-04-24 19:57:17
% DurationCPUTime: 0.18s
% Computational Cost: add. (188->53), mult. (328->86), div. (10->3), fcn. (48->4), ass. (0->39)
t36 = pkin(1) ^ 2;
t25 = cos(qJ(1));
t57 = pkin(2) * t25;
t48 = pkin(1) * t57;
t51 = -0.2e1 * t48 + t36;
t27 = pkin(3) - pkin(4);
t55 = (pkin(2) + t27) * (pkin(2) - t27);
t26 = pkin(3) + pkin(4);
t56 = (pkin(2) + t26) * (pkin(2) - t26);
t59 = (t51 + t56) * (t51 + t55);
t24 = sin(qJ(1));
t58 = pkin(1) * t24;
t12 = pkin(1) * t25 - pkin(2);
t37 = sqrt(-t59);
t3 = t12 * t37;
t54 = t27 ^ 2 * t26 ^ 2;
t29 = qJD(1) ^ 2;
t34 = pkin(2) ^ 2;
t53 = t29 * t34;
t19 = t25 ^ 2;
t52 = t36 * t19;
t30 = pkin(4) ^ 2;
t50 = -t30 / 0.2e1 + t34;
t31 = pkin(3) ^ 2;
t49 = t30 - t31;
t11 = t34 + t51;
t47 = 0.1e1 / t11 ^ 2 * t53;
t46 = t31 / 0.2e1 + t50;
t45 = t34 - t49;
t44 = 0.1e1 / pkin(3) / t37 * t47;
t43 = -t47 / t59 / 0.2e1;
t38 = pkin(1) * t36;
t35 = t36 ^ 2;
t33 = t34 ^ 2;
t28 = 0.2e1 * t34;
t14 = t30 + t31;
t2 = t3 + (t11 + t49) * t58;
t1 = t3 - (t45 + t51) * t58;
t4 = [0, 0, 0, 0, 0, t29 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 * t43, -((0.4e1 * t46 * t52 - 0.4e1 * (t36 + t46) * t48 + t35 + (t28 + t49) * t36 + t34 * t45) * t37 + 0.2e1 * (t11 * t14 - t54) * t12 * t58) * t44 / 0.2e1, -(0.6e1 * ((t34 - t31 / 0.6e1 - t30 / 0.6e1) * t36 + t33 + (-0.5e1 / 0.6e1 * t31 - 0.5e1 / 0.6e1 * t30) * t34 + t54 / 0.6e1) * t52 + 0.3e1 / 0.2e1 * t35 * t34 + (0.3e1 / 0.2e1 * t33 - t14 * t34 - t54 / 0.2e1) * t36 + t34 * t55 * t56 / 0.2e1 + (t24 * t27 * t26 * t3 - 0.3e1 * (t35 + (t28 - 0.2e1 / 0.3e1 * t31 - 0.2e1 / 0.3e1 * t30) * t36 + t33 + (-0.4e1 / 0.3e1 * t31 - 0.4e1 / 0.3e1 * t30) * t34 + t54 / 0.3e1) * t57) * pkin(1) + (-0.4e1 * (-t31 / 0.2e1 + t50) * t19 * t57 + t38 / 0.2e1) * t38) * t44, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t1 ^ 2 * t43, 0, 0, 0, 0;];
T_reg = t4;
