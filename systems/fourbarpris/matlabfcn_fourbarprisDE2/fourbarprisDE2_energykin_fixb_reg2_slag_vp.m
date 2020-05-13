% Calculate inertial parameters regressor of fixed base kinetic energy for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_reg [1x(1*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fourbarprisDE2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energykin_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:16
% EndTime: 2020-05-07 09:45:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (101->14), mult. (65->24), div. (25->3), fcn. (2->2), ass. (0->14)
t14 = (qJ(1) + pkin(3));
t20 = -pkin(2) + t14;
t21 = -pkin(2) - t14;
t27 = ((pkin(1) - t20) * (pkin(1) + t21) * (pkin(1) - t21) * (pkin(1) + t20));
t26 = 1 / t27;
t17 = (t14 ^ 2);
t16 = (qJ(1) ^ 2);
t3 = pkin(1) ^ 2 - pkin(2) ^ 2 - t16 - (2 * qJ(1) + pkin(3)) * pkin(3);
t24 = 1 / t17 * t3 ^ 2;
t15 = (qJD(1) ^ 2);
t19 = (t15 * t26);
t18 = t19 * t24;
t1 = -t18 / 0.2e1;
t2 = [0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, t15 * t3 / t14 * (-t27) ^ (-0.1e1 / 0.2e1), 0, -qJ(1) * t18, (-(t16 * t24 * t26) / 0.2e1 + 0.1e1 / 0.2e1) * t15, 0, 0, 0, 0, 0, -2 * t17 * t19, 0, 0, 0, 0;];
T_reg = t2;
