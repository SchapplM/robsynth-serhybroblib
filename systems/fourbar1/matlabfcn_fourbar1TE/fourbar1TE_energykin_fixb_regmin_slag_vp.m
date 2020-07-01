% Calculate minimal parameter regressor of fixed base kinetic energy for
% fourbar1TE
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
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbar1TE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_energykin_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_energykin_fixb_regmin_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:20:58
% EndTime: 2020-06-26 17:20:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (188->53), mult. (324->85), div. (10->3), fcn. (48->4), ass. (0->38)
t36 = pkin(1) ^ 2;
t25 = cos(qJ(1));
t56 = pkin(2) * t25;
t48 = pkin(1) * t56;
t51 = -0.2e1 * t48 + t36;
t27 = pkin(3) - pkin(4);
t54 = (pkin(2) + t27) * (pkin(2) - t27);
t26 = pkin(3) + pkin(4);
t55 = (pkin(2) + t26) * (pkin(2) - t26);
t58 = (t51 + t55) * (t51 + t54);
t24 = sin(qJ(1));
t57 = pkin(1) * t24;
t12 = pkin(1) * t25 - pkin(2);
t37 = sqrt(-t58);
t3 = t12 * t37;
t53 = t27 ^ 2 * t26 ^ 2;
t19 = t25 ^ 2;
t52 = t36 * t19;
t30 = pkin(4) ^ 2;
t34 = pkin(2) ^ 2;
t50 = -t30 / 0.2e1 + t34;
t31 = pkin(3) ^ 2;
t49 = t30 - t31;
t11 = t34 + t51;
t29 = qJD(1) ^ 2;
t47 = 0.1e1 / t11 ^ 2 * t29 * t34;
t46 = t31 / 0.2e1 + t50;
t45 = t34 - t49;
t44 = 0.1e1 / pkin(3) / t37 * t47;
t43 = -t47 / t58 / 0.2e1;
t38 = pkin(1) * t36;
t35 = t36 ^ 2;
t33 = t34 ^ 2;
t28 = 0.2e1 * t34;
t14 = t30 + t31;
t2 = t3 + (t11 + t49) * t57;
t1 = t3 - (t45 + t51) * t57;
t4 = [t29 / 0.2e1, 0, 0, t2 ^ 2 * t43, -((0.4e1 * t46 * t52 - 0.4e1 * (t36 + t46) * t48 + t35 + (t28 + t49) * t36 + t34 * t45) * t37 + 0.2e1 * (t11 * t14 - t53) * t12 * t57) * t44 / 0.2e1, -(0.6e1 * ((t34 - t31 / 0.6e1 - t30 / 0.6e1) * t36 + t33 + (-0.5e1 / 0.6e1 * t31 - 0.5e1 / 0.6e1 * t30) * t34 + t53 / 0.6e1) * t52 + 0.3e1 / 0.2e1 * t35 * t34 + (0.3e1 / 0.2e1 * t33 - t14 * t34 - t53 / 0.2e1) * t36 + t34 * t54 * t55 / 0.2e1 + (t24 * t27 * t26 * t3 - 0.3e1 * (t35 + (t28 - 0.2e1 / 0.3e1 * t31 - 0.2e1 / 0.3e1 * t30) * t36 + t33 + (-0.4e1 / 0.3e1 * t31 - 0.4e1 / 0.3e1 * t30) * t34 + t53 / 0.3e1) * t56) * pkin(1) + (-0.4e1 * (-t31 / 0.2e1 + t50) * t19 * t56 + t38 / 0.2e1) * t38) * t44, t1 ^ 2 * t43, 0, 0;];
T_reg = t4;
