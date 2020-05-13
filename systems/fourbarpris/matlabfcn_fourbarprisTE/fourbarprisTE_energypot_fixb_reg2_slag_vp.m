% Calculate inertial parameters regressor of potential energy for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% U_reg [1x(1*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbarprisTE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energypot_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:38
% EndTime: 2020-05-07 09:01:38
% DurationCPUTime: 0.05s
% Computational Cost: add. (135->24), mult. (103->31), div. (14->3), fcn. (14->2), ass. (0->19)
t35 = qJ(1) ^ 2;
t51 = (-pkin(2) ^ 2 + pkin(3) ^ 2 + t35);
t48 = qJ(1) + pkin(3);
t45 = -pkin(2) + t48;
t46 = -pkin(2) - t48;
t31 = sqrt(-(pkin(1) + t46) * (pkin(1) + t45) * (pkin(1) - t45) * (pkin(1) - t46));
t50 = t31 * g(1);
t40 = 0.1e1 / pkin(1);
t49 = 0.1e1 / t48 * t40;
t44 = t49 / 0.2e1;
t43 = 0.1e1 / pkin(2) * t40 / 0.2e1;
t42 = -2 * pkin(3) * qJ(1) - t51;
t39 = pkin(1) ^ 2;
t33 = t39 - t42;
t32 = t39 + t42;
t30 = t31 * g(2);
t29 = -g(2) * t33 + t50;
t28 = (g(1) * t33 + t30) * t44;
t1 = [0, 0, 0, 0, 0, 0, t28, -t29 * t49 / 0.2e1, -g(3), -pkin(1) * g(1), 0, 0, 0, 0, 0, 0, t29 * t44, g(3), t28, (0.2e1 * (t35 - t39) * pkin(3) * g(1) + (t30 + (-t39 + t51) * g(1)) * qJ(1)) * t44, 0, 0, 0, 0, 0, 0, (-g(1) * t32 + t30) * t43, (-g(2) * t32 - t50) * t43, -g(3), 0;];
U_reg = t1;
