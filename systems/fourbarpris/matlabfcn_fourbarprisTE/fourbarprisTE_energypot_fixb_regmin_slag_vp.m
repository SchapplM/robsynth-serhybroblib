% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:07
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbarprisTE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energypot_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:07:02
% EndTime: 2020-06-27 17:07:02
% DurationCPUTime: 0.05s
% Computational Cost: add. (132->21), mult. (102->30), div. (14->3), fcn. (14->2), ass. (0->19)
t34 = qJ(1) ^ 2;
t50 = (-pkin(2) ^ 2 + pkin(3) ^ 2 + t34);
t47 = qJ(1) + pkin(3);
t44 = -pkin(2) + t47;
t45 = -pkin(2) - t47;
t30 = sqrt(-(pkin(1) + t45) * (pkin(1) + t44) * (pkin(1) - t44) * (pkin(1) - t45));
t49 = t30 * g(1);
t39 = 0.1e1 / pkin(1);
t48 = 0.1e1 / t47 * t39;
t43 = t48 / 0.2e1;
t42 = 0.1e1 / pkin(2) * t39 / 0.2e1;
t41 = -2 * pkin(3) * qJ(1) - t50;
t38 = pkin(1) ^ 2;
t32 = t38 - t41;
t31 = t38 + t41;
t29 = t30 * g(2);
t28 = -g(2) * t32 + t49;
t27 = (g(1) * t32 + t29) * t43;
t1 = [0, t27, -t28 * t48 / 0.2e1, t28 * t43, t27, (0.2e1 * (t34 - t38) * pkin(3) * g(1) + (t29 + (-t38 + t50) * g(1)) * qJ(1)) * t43, 0, (-g(1) * t31 + t29) * t42, (-g(2) * t31 - t49) * t42;];
U_reg = t1;
